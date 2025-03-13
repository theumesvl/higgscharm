#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022postEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022postEE/MuonEGF --dataset MuonEGF_6 --partition_fileset '{"MuonEGF_6": ["root://ruhex-osgce.rutgers.edu//store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/15a5f576-51d7-44ce-9ee4-982ad3c8a25f.root", "root://maite.iihe.ac.be:1095//store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/259ac7b2-f513-45eb-9f7b-2d5d920da551.root", "root://cmsio2.rc.ufl.edu:1094//store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/2ce0d3a3-84d6-47d5-b263-2a1a20cc3ac6.root", "root://cmsio2.rc.ufl.edu:1094//store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/367b5515-6325-4a22-a724-797d6a6a4e66.root", "root://eoscms.cern.ch//eos/cms/store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/3cb0f11a-68d2-4e59-93db-7ab5c719ea8a.root", "root://ruhex-osgce.rutgers.edu//store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/42a95a78-2513-404a-a0ae-3b7544d7ca3f.root", "root://ruhex-osgce.rutgers.edu//store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/498844ea-65a3-4ca5-931a-11ea6a19e79c.root", "root://hactar01.crc.nd.edu//store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/4b274c49-1fef-4222-aae7-d96e86014733.root", "root://hactar01.crc.nd.edu//store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/4b4c2ad7-5747-44d3-93eb-a2f4b695b473.root"]}' --output_format root