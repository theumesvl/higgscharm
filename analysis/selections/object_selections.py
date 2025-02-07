import inspect
import numpy as np
import awkward as ak
from coffea.nanoevents.methods import candidate
from coffea.nanoevents.methods.vector import LorentzVector
from analysis.working_points import working_points
from analysis.selections import (
    delta_r_higher,
    delta_r_lower,
    select_dileptons,
    transverse_mass
)


class ObjectSelector:

    def __init__(self, object_selection_config, year):
        self.object_selection_config = object_selection_config
        self.year = year

    def select_objects(self, events):
        self.objects = {}
        self.events = events
        for obj_name, obj_config in self.object_selection_config.items():
            # check if object field is read from events or from user defined function
            if "events" in obj_config["field"]:
                self.objects[obj_name] = eval(obj_config["field"])
            else:
                selection_function = getattr(self, obj_config["field"])
                parameters = inspect.signature(selection_function).parameters.keys()
                if "cuts" in parameters:
                    selection_function(obj_config["cuts"])
                    break
                else:
                    selection_function()
            if "cuts" in obj_config:
                selection_mask = self.get_selection_mask(
                    events=events, obj_name=obj_name, cuts=obj_config["cuts"]
                )
                self.objects[obj_name] = self.objects[obj_name][selection_mask]
        return self.objects

    def get_selection_mask(self, events, obj_name, cuts):
        # bring 'objects' and to local scope
        objects = self.objects
        # initialize selection mask
        selection_mask = ak.ones_like(self.objects[obj_name].pt, dtype=bool)
        # iterate over all cuts
        for selection, str_mask in cuts.items():
            # check if 'str_mask' contains 'events' or 'objects' and evaluate string expression
            if "events" in str_mask or "objects" in str_mask:
                mask = eval(str_mask)
            # read the mask from the working points function
            else:
                signature = inspect.signature(getattr(working_points, selection))
                parameters = signature.parameters.keys()
                if "year" in parameters:
                    mask = getattr(working_points, selection)(
                        self.events, str_mask, self.year
                    )
                else:
                    mask = getattr(working_points, selection)(self.events, str_mask)
            # update selection mask
            selection_mask = np.logical_and(selection_mask, mask)
        return selection_mask

    # --------------------------------------------------------------------------------
    # WWto2L
    # --------------------------------------------------------------------------------
    # 2L = two muons
    def select_dimuons(self):
        if "muons" not in self.objects:
            raise ValueError(f"'muons' object has not been defined!")
        self.objects["dimuons"] = select_dileptons(self.objects, "muons")
    # 2L = two electrons
    def select_dielectrons(self):
        if "electrons" not in self.objects:
            raise ValueError(f"'electrons' object has not been defined!")
        self.objects["dielectrons"] = select_dileptons(self.objects, "electrons")
    # # 2L = muon and electron 
    # def select_muonElectron(self):
    #     if "muons" not in self.objects or "electrons" not in self.objects:
    #         raise ValueError(f"'muons' or 'electrons' object has not been defined!")
    #     self.objects["muonElectron"] = select_dileptons(
    #         self.objects, "muonElectron"
    #     )
    # --------------------------------------------------------------------------------
    # WWTo2L
    # --------------------------------------------------------------------------------
#     # add muons and electrons together in one large pool of leptons
#     def select_WWto2L_leptons(self):
#         leptons = ak.cartesian(
#             [self.objects["muons"], self.objects["electrons"]], axis=1
#         )
#         leptons = leptons[ak.argsort(leptons.pt, axis=1)]
#         self.objects["leptons"] = ak.zip(
#             {
#                 "pt": leptons.pt,
#                 "eta": leptons.eta,
#                 "phi": leptons.phi,
#                 "mass": leptons.mass,
#                 "charge": leptons.charge,
#                 "pdgId": leptons.pdgId,
#             },
#             with_name="PtEtaPhiMCandidate",
#             behavior=candidate.behavior,
#         )
    
#     # select two leptons 
#     def select_WWto2L_ll_pairs(self):
#         self.objects["ll_pairs"] = ak.combinations(
#             self.objects["leptons"], 2, fields=["l1", "l2"]
#         )
#         self.objects["ll_pairs"].pt = (
#             self.objects["ll_pairs"].l1.pt + self.objects["ll_pairs"].l2.pt
#         )
    # --------------------------------------------------------------------------------
    # WWto2L 
    # --------------------------------------------------------------------------------
    def select_hww_leptons(self):
        # set 'leptons' by concatenating electrons and muons
        leptons = ak.concatenate(
            [self.objects["muons"], self.objects["electrons"]], axis=1
        )
        leptons = leptons[ak.argsort(leptons.pt, axis=1)]
        self.objects["leptons"] = ak.zip(
            {
                "pt": leptons.pt,
                "eta": leptons.eta,
                "phi": leptons.phi,
                "mass": leptons.mass,
                "charge": leptons.charge,
                "pdgId": leptons.pdgId,
            },
            with_name="PtEtaPhiMCandidate",
            behavior=candidate.behavior,
        )
        
    def select_hww_pairs(self):
        self.objects["ll_pairs"] = ak.combinations(self.objects["leptons"], 2, fields=["l1", "l2"])
        self.objects["ll_pairs"].pt = (
            self.objects["ll_pairs"].l1.pt + self.objects["ll_pairs"].l2.pt
        )
        
        
    # def select_hww_pTll(self):
    #     self.objects["ll_pair"] = self.objects["ll_pairs"].l1 + self.objects["ll_pairs"].l2
    #     self.objects["pTll"] = self.objects["ll_pair"].pt
        
    def select_hww_mTll(self):
        self.objects["mTll"] = transverse_mass(self.objects["ll_pairs"].l1+self.objects["ll_pairs"].l2, self.objects["met"])
        
    def select_hww_mTl1(self):
        self.objects["mTl1"] = transverse_mass(self.objects["ll_pairs"].l1, self.objects["met"])
        
    def select_hww_mTl2(self):
        self.objects["mTl2"] = transverse_mass(self.objects["ll_pairs"].l2, self.objects["met"])
        
    def select_candidate_cjet(self):
        self.objects["candidate_cjet"] = self.objects["cjets"][ak.argmax(self.objects["cjets"].btagDeepFlavCvL, axis=1) == ak.local_index(self.objects["cjets"], axis=1)]

        #not to be used as we cannot reconstruct the W mass due to the missing energy
#     def select_WWto2L_WWpairs(self):
#         # sort lepton pairs by proximity to W mass
#         wmass = 80.3602
#         dist_from_w_all_pairs = np.abs(
#             (self.objects["ll_pairs"].l1 + self.objects["ll_pairs"].l2).mass - wmass
#         )
#         sorted_ll_pairs = self.objects["ll_pairs"][
#             ak.argsort(dist_from_w_all_pairs, axis=1)
#         ]
#         # WW candidates
#         ww_pairs = select_WWtoMuEle_WW_candidates(sorted_ll_pairs)
#         # mass-sorted alternative pairing W candidates
#         alt_sorted_ll_pairs = self.objects["ll_pairs"][
#             ak.argsort(
#                 -(self.objects["ll_pairs"].l1 + self.objects["ll_pairs"].l2).mass,
#                 axis=1,
#             )
#         ]
#         alt_ww_pairs = select_WWtoMuEle_WW_candidates(alt_sorted_ll_pairs)
#         # 'smart cut': require NOT(|mZa - mZ| < |mw1 − mZ| AND mZb < 12)
#         # This cut discards 4µ and 4e candidates where the alternative pairing looks like an on-shell Z + low-mass l+l−
#         #smart_cut = ~(
#           #  (
#            #     np.abs(alt_ww_pairs.w1.p4.mass - wmass)
#           #      < np.abs(ww_pairs.w1.p4.mass - wmass)
#           #  )
#           #  & (alt_ww_pairs.w2.p4.mass < 12)
#         #)
#         #self.objects["ww_pairs"] = ww_pairs[smart_cut]
#         self.objects["ww_pairs"].pt = (
#             self.objects["ww_pairs"].w1.p4.pt + self.objects["ww_pairs"].w2.p4.pt
#         )

#     def select_WWtoMuEle_WWcandidate(self):
#         """
#         selects best ww candidate as the one with w1 closest in mass to nominal W boson mass
#         and w2 from the candidates whose lepton give higher pT sum
#         """
#         # get mask of w1's closest to W
#         wmass = 80.3602
#         w1_dist_to_w = np.abs(self.objects["ww_pairs"].w1.p4.mass - wmass)
#         min_w1_dist_to_w = ak.min(w1_dist_to_w, axis=1)
#         closest_w1_mask = w1_dist_to_w == min_w1_dist_to_w
#         # get mask of w2's with higher pT sum
#         w2_pt_sum = (
#             self.objects["ww_pairs"].w2.l1.pt + self.objects["ww_pairs"].w2.l2.pt
#         )
#         max_w2_pt_sum = ak.max(w2_pt_sum[closest_w1_mask], axis=1)
#         best_candidate_mask = (w2_pt_sum == max_w2_pt_sum) & closest_w1_mask
#         # select best candidate from ww_pairs
#         self.objects["ww_candidate"] = self.objects["ww_pairs"][best_candidate_mask]
