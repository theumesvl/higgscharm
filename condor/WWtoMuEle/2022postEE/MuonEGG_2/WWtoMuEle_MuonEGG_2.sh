#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022postEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022postEE/MuonEGG --dataset MuonEGG_2 --partition_fileset '{"MuonEGG_2": ["root://eos01.grid.cyfronet.pl:1094//eos/cms/store/data/Run2022G/MuonEG/NANOAOD/22Sep2023-v1/2520000/59ff0fa7-935a-44ec-b939-41d38bad11f5.root", "root://hactar01.crc.nd.edu//store/data/Run2022G/MuonEG/NANOAOD/22Sep2023-v1/2520000/62dc133c-e2af-47fc-9d46-afed2f37e3e0.root", "root://ruhex-osgce.rutgers.edu//store/data/Run2022G/MuonEG/NANOAOD/22Sep2023-v1/2520000/654f214d-205a-4479-af44-73a1f2440f37.root", "root://cmsio2.rc.ufl.edu:1094//store/data/Run2022G/MuonEG/NANOAOD/22Sep2023-v1/2520000/68b813ce-362d-4ec7-8027-6d532727e29e.root", "root://hactar01.crc.nd.edu//store/data/Run2022G/MuonEG/NANOAOD/22Sep2023-v1/2520000/7a3c45c1-2740-4e6a-bb92-1aecda97f463.root", "root://hactar01.crc.nd.edu//store/data/Run2022G/MuonEG/NANOAOD/22Sep2023-v1/2520000/8c0a6cf3-0e73-4d9a-bc59-f0402faef283.root", "root://cmsio2.rc.ufl.edu:1094//store/data/Run2022G/MuonEG/NANOAOD/22Sep2023-v1/2520000/93e68f59-35ed-4c55-bc05-1fcbc75c249a.root", "root://k8s-redir.ultralight.org:1094//store/data/Run2022G/MuonEG/NANOAOD/22Sep2023-v1/2520000/9478336c-d059-4a62-baf6-5e6e0f1e8d52.root", "root://k8s-redir.ultralight.org:1094//store/data/Run2022G/MuonEG/NANOAOD/22Sep2023-v1/2520000/ba292b38-e8aa-44c8-9d6c-65f5a11f20b5.root"]}' --output_format root