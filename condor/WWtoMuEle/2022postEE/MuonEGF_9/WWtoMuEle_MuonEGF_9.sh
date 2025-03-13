#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022postEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022postEE/MuonEGF --dataset MuonEGF_9 --partition_fileset '{"MuonEGF_9": ["root://eoscms.cern.ch//eos/cms/store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/bb3128ba-cf10-4420-b83b-5e3529056009.root", "root://eoscms.cern.ch//eos/cms/store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/be04cc0b-f711-44e6-a8cb-61405e4961ad.root", "root://maite.iihe.ac.be:1095//store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/bea537cc-9b90-4fb8-a8a2-94c28c13225a.root", "root://hactar01.crc.nd.edu//store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/c7a5b077-7fd0-4a8d-9f6e-2fb719580e69.root", "root://hactar01.crc.nd.edu//store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/c8dd1e07-1d1e-47ef-a0dd-9ca291c99c83.root", "root://cmsio2.rc.ufl.edu:1094//store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/d1f1107d-bef8-4fc3-a7a2-e6818fe6bb8a.root", "root://maite.iihe.ac.be:1095//store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/d3b42279-1c8f-4455-ada5-9f6e24a5cbca.root", "root://hactar01.crc.nd.edu//store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/d98ef7a2-9a1f-4b28-86b1-007fc4d82224.root", "root://hactar01.crc.nd.edu//store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/dd5ee85e-c027-47ce-bf5c-018872410e25.root"]}' --output_format root