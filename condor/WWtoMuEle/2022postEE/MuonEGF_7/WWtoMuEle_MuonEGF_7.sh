#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022postEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022postEE/MuonEGF --dataset MuonEGF_7 --partition_fileset '{"MuonEGF_7": ["root://eoscms.cern.ch//eos/cms/store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/4d76213a-ef14-411a-9558-559a6df3f978.root", "root://cmsio2.rc.ufl.edu:1094//store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/4e81fe81-8f4c-4436-8c53-8eda162f6923.root", "root://ruhex-osgce.rutgers.edu//store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/4fb72196-3b02-4499-8f6c-a54e15692b32.root", "root://cmsio2.rc.ufl.edu:1094//store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/56f478ca-3b6e-4ade-879d-4473b6503182.root", "root://hactar01.crc.nd.edu//store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/57ca02b2-a056-41de-9afb-7f9a2a55d1e7.root", "root://hactar01.crc.nd.edu//store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/5912311b-cdcd-42be-a1ff-03c0111b317d.root", "root://hactar01.crc.nd.edu//store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/5fc74264-b04b-4daf-9052-ca0b33b2e767.root", "root://eoscms.cern.ch//eos/cms/store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/695bf65e-2d27-4dbc-b8a5-7e82a80eb46a.root", "root://eoscms.cern.ch//eos/cms/store/data/Run2022F/MuonEG/NANOAOD/22Sep2023-v1/50000/7a9d9efd-464d-41c7-a09c-ae3ead628f3f.root"]}' --output_format root