import numpy as np
import awkward as ak
import dask_awkward as dak
from copy import deepcopy
from coffea import processor
from coffea.nanoevents import PFNanoAODSchema
from coffea.lumi_tools import LumiData, LumiList
from coffea.analysis_tools import Weights, PackedSelection
from analysis.configs import load_config
from analysis.corrections.pileup import add_pileup_weight
from analysis.corrections.electron import ElectronWeights, ElectronSS
from analysis.histograms import HistBuilder, fill_histogram
from analysis.selections import (
    object_selector,
    get_lumi_mask,
    get_trigger_mask,
    get_trigger_match_mask,
)


PFNanoAODSchema.warn_missing_crossrefs = False


class ZToEEProcessor(processor.ProcessorABC):
    def __init__(self, year: str):

        self.year = year
        self.processor_config = load_config(
            config_type="processor", config_name="ztoee", year=year
        )
        self.histogram_config = self.processor_config.histogram_config
        self.histograms = HistBuilder(self.histogram_config).build_histogram()

    def process(self, events):
        # get goldenjson and hlt paths
        goldenjson = self.processor_config.goldenjson
        hlt_paths = self.processor_config.hlt_paths
        # check if dataset is MC or Data
        is_mc = hasattr(events, "genWeight")
        # initialize output dictionary
        output = {}
        # initialize metadata info
        nevents = ak.num(events, axis=0)
        output["metadata"] = {}
        output["metadata"].update({"raw_initial_nevents": nevents})

        # -------------------------------------------------------------
        # Object corrections
        # -------------------------------------------------------------
        # electron scale and smearing corrections
        electron_ss = ElectronSS(
            events=events,
            year=self.year,
            variation="nominal",
        )
        if is_mc:
            # energies in MC are smeared
            electron_ss.apply_smearing()
        else:
            # energies in data are scaled
            electron_ss.apply_scale()
        # --------------------------------------------------------------
        # Weights
        # --------------------------------------------------------------
        # initialize weights container
        weights_container = Weights(None, storeIndividual=True)
        if is_mc:
            # add genweights
            weights_container.add("genweight", events.genWeight)
            # add pileup weights
            add_pileup_weight(
                events=events,
                year=self.year,
                variation="nominal",
                weights_container=weights_container,
            )
            # add electron id, reco and trigger weights
            electron_weights = ElectronWeights(
                events=events,
                year=self.year,
                weights=weights_container,
                variation="nominal",
                id_wp=self.processor_config.object_selection["electrons"]["cuts"][
                    "electron_id"
                ],
                hlt_paths=hlt_paths,
            )
            electron_weights.add_id_weights()
            electron_weights.add_hlt_weights()
            electron_weights.add_reco_weights("RecoBelow20")
            electron_weights.add_reco_weights("Reco20to75")
            electron_weights.add_reco_weights("RecoAbove75")
        else:
            weights_container.add("genweight", ak.ones_like(events.PV.npvsGood))

        # save nevents (sum of weights) before selections
        sumw = ak.sum(weights_container.weight())
        output["metadata"].update({"sumw": sumw})

        # --------------------------------------------------------------
        # Object selection
        # --------------------------------------------------------------
        objects = object_selector(events, self.processor_config.object_selection)

        # --------------------------------------------------------------
        # Event selection
        # --------------------------------------------------------------
        selections = {}
        for selection, str_mask in self.processor_config.event_selection.items():
            selections[selection] = eval(str_mask)

        event_selection = PackedSelection()
        event_selection.add_multiple(selections)
        region_selection = event_selection.all(*selections.keys())

        # save cutflow
        output["metadata"].update({"cutflow": {"initial": sumw}})
        current_selection = []
        for cut_name in selections.keys():
            current_selection.append(cut_name)
            output["metadata"]["cutflow"][cut_name] = ak.sum(
                weights_container.weight()[event_selection.all(*current_selection)]
            )
        # save raw and weighted number of events after selection
        final_nevents = dak.sum(region_selection)
        weighted_final_nevents = ak.sum(weights_container.weight()[region_selection])
        output["metadata"].update(
            {
                "weighted_final_nevents": weighted_final_nevents,
                "raw_final_nevents": final_nevents,
            }
        )
        if not is_mc:
            # compute integrated luminosity (/pb)
            lumi_mask = selections["lumimask"]
            lumi_data = LumiData(self.processor_config.lumidata)
            lumi_list = LumiList(
                events[lumi_mask].run, events[lumi_mask].luminosityBlock
            )
            lumi = lumi_data.get_lumi(lumi_list)
            # save luminosity to metadata
            output["metadata"].update({"lumi": lumi})

        # --------------------------------------------------------------
        # Histogram filling
        # --------------------------------------------------------------
        if final_nevents > 0:
            histograms = deepcopy(self.histograms)
            # get analysis features
            feature_map = {}
            for feature, axis_info in self.histogram_config.axes.items():
                feature_map[feature] = eval(axis_info["expression"])[region_selection]
            # fill histograms
            if is_mc:
                # get event weight systematic variations for MC samples
                variations = ["nominal"] + list(weights_container.variations)
                for variation in variations:
                    if variation == "nominal":
                        region_weight = weights_container.weight()[region_selection]
                    else:
                        region_weight = weights_container.weight(modifier=variation)[
                            region_selection
                        ]

                    fill_histogram(
                        histograms=histograms,
                        histogram_config=self.histogram_config,
                        feature_map=feature_map,
                        weights=region_weight,
                        variation=variation,
                        flow=True,
                    )
            else:
                region_weight = weights_container.weight()[region_selection]
                fill_histogram(
                    histograms=histograms,
                    histogram_config=self.histogram_config,
                    feature_map=feature_map,
                    weights=region_weight,
                    variation="nominal",
                    flow=True,
                )
        # add histograms to output dictionary
        output["histograms"] = histograms
        return output

    def postprocess(self, accumulator):
        pass
