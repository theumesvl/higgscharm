#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022postEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022postEE/MuonEGF --dataset MuonEGF_8 --partition_fileset '{"MuonEGF_8": ["root://maite.iihe.ac.be:1095//store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/7d61fed9-c008-4182-b708-fea93ea17c37.root", "root://maite.iihe.ac.be:1095//store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/808acd24-1cba-4794-96e3-e561fca19be4.root", "root://hactar01.crc.nd.edu//store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/86532c24-07a2-4838-8511-6c548477d075.root", "root://maite.iihe.ac.be:1095//store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/9337d1ce-bde9-496f-8545-cd2cb9c48f24.root", "root://eoscms.cern.ch//eos/cms/store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/94ab9cd3-8cb1-454f-ae0a-fde4b3e6ef37.root", "root://eoscms.cern.ch//eos/cms/store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/9959fdac-5491-4ed5-974a-9bf656ff8fb7.root", "root://eoscms.cern.ch//eos/cms/store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/a1667d31-5294-4ed8-a293-4d9105a8f443.root", "root://cmsio2.rc.ufl.edu:1094//store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/ac63ab54-af8e-4d0d-829f-2998dc999de4.root", "root://hactar01.crc.nd.edu//store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/b11afa8f-ba1e-4d78-a3ef-aaee3fa86f01.root"]}' --output_format root