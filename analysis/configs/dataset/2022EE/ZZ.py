from analysis.configs.dataset_config import DatasetConfig

dataset_config = DatasetConfig(
    name="ZZ",
    process="Diboson",
    path=(
        "/pnfs/iihe/cms/store/user/daocampo/PFNano_Run3/"
        "mc_summer22EE_MINIAODv4/ZZ_TuneCP5_13p6TeV_pythia8/"
        "Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2_BTV_Run3_2022_Comm_MINIAODv4/240518_172329/0000/"
    ),
    key="Events",
    year="2022EE",
    era="MC",
    xsec=12.75,
    partitions=3,
    stepsize=50_000,
    filenames=(
        "MC_defaultAK4_29.root",
        "MC_defaultAK4_25.root",
        "MC_defaultAK4_27.root",
        "MC_defaultAK4_10.root",
        "MC_defaultAK4_21.root",
        "MC_defaultAK4_16.root",
        "MC_defaultAK4_28.root",
        "MC_defaultAK4_13.root",
        "MC_defaultAK4_12.root",
        "MC_defaultAK4_11.root",
        "MC_defaultAK4_22.root",
        "MC_defaultAK4_26.root",
        "MC_defaultAK4_15.root",
        "MC_defaultAK4_14.root",
        "MC_defaultAK4_24.root",
        "MC_defaultAK4_23.root",
        "MC_defaultAK4_17.root",
        "MC_defaultAK4_20.root",
        "MC_defaultAK4_1.root",
        "MC_defaultAK4_2.root",
        "MC_defaultAK4_3.root",
        "MC_defaultAK4_4.root",
        "MC_defaultAK4_5.root",
        "MC_defaultAK4_6.root",
        "MC_defaultAK4_7.root",
        "MC_defaultAK4_8.root",
        "MC_defaultAK4_9.root",
    ),
)