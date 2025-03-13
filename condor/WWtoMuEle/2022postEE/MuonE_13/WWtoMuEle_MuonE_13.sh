#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022postEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022postEE/MuonE --dataset MuonE_13 --partition_fileset '{"MuonE_13": ["root://k8s-redir.ultralight.org:1094//store/data/Run2022E/Muon/NANOAOD/22Sep2023-v1/30000/361b56fd-35bd-4146-986e-d2123c70b34b.root", "root://hactar01.crc.nd.edu//store/data/Run2022E/Muon/NANOAOD/22Sep2023-v1/30000/38cd861e-21ec-4675-88b6-f455cfb723f9.root", "root://hactar01.crc.nd.edu//store/data/Run2022E/Muon/NANOAOD/22Sep2023-v1/30000/4a94c3da-488c-4574-a530-c48912188852.root", "root://cceos.ihep.ac.cn:1094//eos/ihep/cms/store/data/Run2022E/Muon/NANOAOD/22Sep2023-v1/30000/4b0c8ae0-6947-4cf3-a214-ae279c22d0f3.root", "root://eos.cms.rcac.purdue.edu//store/data/Run2022E/Muon/NANOAOD/22Sep2023-v1/30000/4d5edfe2-4b0d-4aa3-8252-b764a96dbdc8.root", "root://grid-cms-xrootd.physik.rwth-aachen.de:1094//store/data/Run2022E/Muon/NANOAOD/22Sep2023-v1/30000/580048cb-5685-4b55-b59d-3e4090d001bc.root", "root://ruhex-osgce.rutgers.edu//store/data/Run2022E/Muon/NANOAOD/22Sep2023-v1/30000/5873413e-fc09-4f9f-9bfa-bb16fa4d17e9.root", "root://eos.cms.rcac.purdue.edu//store/data/Run2022E/Muon/NANOAOD/22Sep2023-v1/30000/58f38b62-7c0d-478f-accb-b9cbd9a729b2.root", "root://hactar01.crc.nd.edu//store/data/Run2022E/Muon/NANOAOD/22Sep2023-v1/30000/6951968d-8a34-4e6d-9847-ea7025c92592.root"]}' --output_format root