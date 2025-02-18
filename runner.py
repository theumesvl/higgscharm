import os
import argparse


ERAS = {"2022preEE": ["C", "D"], "2022postEE": ["E", "F", "G"]}
PRIMARY_DATASETS = ["Muon", "MuonEG", "EGamma"]
DATA_SAMPLES = {}
for year, eras in ERAS.items():
    DATA_SAMPLES[year] = {}
    for primary_dataset in PRIMARY_DATASETS:
        DATA_SAMPLES[year][primary_dataset] = []
        for era in eras:
            DATA_SAMPLES[year][primary_dataset].append(f"{primary_dataset}{era}")

            
MC_DATASETS = {
    "ttbar": [
        "TTto2L2Nu", 
        "TTto4Q", 
        "TTtoLNu2Q"
    ],
    "singletop": [
        "TbarWplusto2L2Nu",
        "TWminusto2L2Nu",
        "TWminustoLNu2Q",
        "TbarWplusto4Q",
        "TbarWplustoLNu2Q",
        "TWminusto4Q",
        "TbarBQ",
        "TBbarQ",
        "TbarQto2Q",
        "TbarQtoLNu",
        "TQbarto2Q",
        "TQbartoLNu",
        "TbarBtoLminusNuB",
        "TBbartoLplusNuBbar",
    ],
    "diboson": [
        #"WW",
        "ZZto2L2Nu",
        "ZZto2Nu2Q",
        "ZZto2L2Q",
        "ZZto4L",
        "WZto3LNu",
        "WZto2L2Q",
        "WZtoLNu2Q",
        "WZtoL3Nu-1Jets",
        "WWtoLNu2Q",
        "WWto4Q",
        "WWto2L2Nu",
        "WGtoLNuG-PTG-10to100",
        "WGtoLNuG-PTG-100to200",
        "WGtoLNuG-PTG-200to400",
        "WGtoLNuG-PTG-400to600",
        "WGtoLNuG-PTG-600",
        #"WZ", 
        #"ZZ"
    ],
    "DY+jets": [
        "DYto2L_2Jets_50", 
        "DYto2L_2Jets_10to50"
    ],
    "V+jets": [
        #"WtoLNu-2Jets_0J",
        #"WtoLNu-2Jets_1J",
        #"WtoLNu-2Jets_2J",
        "WtoLNu-2Jets",
    ],
    "higgs_bkg_WW": [
        "GluGluHto2Wto2L2Nu",
        "VBFHto2Wto2L2Nu",
    ],
    "higgs": [
        "bbH_Hto2Zto4L",
        "GluGluHtoZZto4L",
        "TTH_Hto2Z",
        "VBFHto2Zto4L",
        "WminusH_Hto2Zto4L",
        "WplusH_Hto2Zto4L",
        "ZHto2Zto4L",
    ],
    "ggtozz": [
        "GluGluToContinto2Zto2E2Mu",
        "GluGluToContinto2Zto2E2Tau",
        "GluGluToContinto2Zto2Mu2Tau",
        "GluGlutoContinto2Zto4E",
        "GluGlutoContinto2Zto4Mu",
        "GluGlutoContinto2Zto4Tau",
    ],
    "qqtozz": ["ZZto4L"],
}
DATASETS = {
    "WWtoMuEle": {
        "mc": ["ttbar", "singletop", "diboson","DY+jets","V+jets","higgs_bkg_WW"],
        "data": ["Muon", "MuonEG", "EGamma"],
    },
    "zzto4l": {
        "mc": ["higgs", "ggtozz", "qqtozz"],
        "data": ["Muon", "MuonEG", "EGamma"],
    },
    "ztoee": {
        "mc": ["ttbar", "singletop", "diboson"], 
        "data": ["EGamma"]
    },
    "ztomumu": {
        "mc": ["ttbar", "singletop", "diboson"], 
        "data": ["EGamma"]
    },
}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--processor",
        dest="processor",
        type=str,
        default="WWtoMuEle",
        help="processor to be used {ztomumu, ztoee, zzto4l, WWtoMuEle} (default WWtoMuEle)",
    )
    parser.add_argument(
        "--year",
        dest="year",
        type=str,
        help="dataset year {2022preEE, 2022postEE}",
    )
    parser.add_argument(
        "--nfiles",
        dest="nfiles",
        type=int,
        default=20,
        help="number of root files to include in each dataset partition (default 20)",
    )
    parser.add_argument(
        "--submit",
        action="store_true",
        help="Enable Condor job submission. If not provided, it just builds condor files",
    )
    parser.add_argument(
        "--eos",
        action="store_true",
        help="Enable saving outputs to /eos",
    )
    args = parser.parse_args()
    # get datasets for processor and year
    mc = [
        sample
        for dataset in DATASETS[args.processor]["mc"]
        for sample in MC_DATASETS[dataset]
    ]
    data = [
        sample
        for dataset in DATASETS[args.processor]["data"]
        for sample in DATA_SAMPLES[args.year][dataset]
    ]
    datasets = mc + data
    # submit job for each dataset
    for dataset in datasets:
        cmd = f"python3 submit_condor.py --processor {args.processor} --year {args.year} --dataset {dataset} --nfiles {args.nfiles}"
        if args.submit:
            cmd += " --submit"

