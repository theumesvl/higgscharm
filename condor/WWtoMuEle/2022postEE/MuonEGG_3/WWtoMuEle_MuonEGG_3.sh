#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022postEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022postEE/MuonEGG --dataset MuonEGG_3 --partition_fileset '{"MuonEGG_3": ["root://eos01.grid.cyfronet.pl:1094//eos/cms/store/data/Run2022G/MuonEG/NANOAOD/22Sep2023-v1/2520000/baf18ae5-5d75-4c55-aade-6587b21c9f10.root", "root://cmsio2.rc.ufl.edu:1094//store/data/Run2022G/MuonEG/NANOAOD/22Sep2023-v1/2520000/be0de64d-8288-40ed-b8c1-aecbf25407ae.root", "root://hactar01.crc.nd.edu//store/data/Run2022G/MuonEG/NANOAOD/22Sep2023-v1/2520000/c1f12ec0-8c3a-4aa8-ba91-1ded3851f628.root", "root://maite.iihe.ac.be:1095//store/data/Run2022G/MuonEG/NANOAOD/22Sep2023-v1/2520000/cd404eb6-8218-4787-b5ed-af6cd9fe3750.root", "root://cmsio2.rc.ufl.edu:1094//store/data/Run2022G/MuonEG/NANOAOD/22Sep2023-v1/2520000/cf518da3-0dbb-4418-9181-1ee1271a5a9b.root", "root://ruhex-osgce.rutgers.edu//store/data/Run2022G/MuonEG/NANOAOD/22Sep2023-v1/2520000/d44742ec-959a-48be-a199-a23f2cf883ef.root", "root://ruhex-osgce.rutgers.edu//store/data/Run2022G/MuonEG/NANOAOD/22Sep2023-v1/2520000/e154069d-bb89-4e63-9461-5193afff31f5.root", "root://hactar01.crc.nd.edu//store/data/Run2022G/MuonEG/NANOAOD/22Sep2023-v1/2520000/f1d8ec12-59c3-4e2b-8b61-adcd3e68449c.root", "root://maite.iihe.ac.be:1095//store/data/Run2022G/MuonEG/NANOAOD/22Sep2023-v1/2520000/fb86661e-a9f1-47c3-ab58-0c2805c1ddbd.root"]}' --output_format root