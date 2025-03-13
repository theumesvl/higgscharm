#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022postEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022postEE/MuonEGF --dataset MuonEGF_10 --partition_fileset '{"MuonEGF_10": ["root://ruhex-osgce.rutgers.edu//store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/e2d56b01-2ce6-48c8-9be8-dccfdc5b06f2.root", "root://eoscms.cern.ch//eos/cms/store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/ece3ac63-510d-4a86-94f5-a82c806d5e78.root", "root://eoscms.cern.ch//eos/cms/store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/ef058380-4eee-4088-9232-526fe0ebbc2b.root", "root://maite.iihe.ac.be:1095//store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/f188bf88-79f9-4b0b-a177-62aee090a59b.root", "root://maite.iihe.ac.be:1095//store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/f1eb623a-0c2a-47b7-b9c4-5e9b804359e7.root", "root://hactar01.crc.nd.edu//store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/f7105902-4b2e-42c3-a8ce-5e53f0bd3e45.root", "root://maite.iihe.ac.be:1095//store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/f8f90fb8-59ab-4e82-b75c-589083d63f04.root", "root://ruhex-osgce.rutgers.edu//store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/f919eb16-e7b5-42fe-bc78-102860f2220d.root", "root://ruhex-osgce.rutgers.edu//store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/fc20cea4-137e-42cc-954f-11e7fb49f92b.root"]}' --output_format root