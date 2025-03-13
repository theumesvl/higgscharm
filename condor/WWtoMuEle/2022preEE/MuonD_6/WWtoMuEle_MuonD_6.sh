#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022preEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022preEE/MuonD --dataset MuonD_6 --partition_fileset '{"MuonD_6": ["root://cceos.ihep.ac.cn:1094//eos/ihep/cms/store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/8c6b64bc-8249-4fa3-bcb5-375ddcecd94e.root", "root://xrootd-local.unl.edu:1094//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/8e06a5ea-9534-4d16-abcf-988426b31fd2.root", "root://xrootd-local.unl.edu:1094//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/8f08b2ae-946d-4b64-9bd5-ef7cc46002c6.root", "root://hactar01.crc.nd.edu//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/92ae04b7-552f-49ee-a98f-9c08fdd7307d.root", "root://eos.cms.rcac.purdue.edu//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/951a05c2-f3c0-4755-9fc5-b102b1d648dd.root", "root://xrootd-local.unl.edu:1094//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/956993a7-adc9-483f-a270-c3ed7d76160f.root", "root://hactar01.crc.nd.edu//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/9bdd38b2-e53d-4db3-9e7c-a822451b1621.root", "root://eos.cms.rcac.purdue.edu//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/9f8b3535-18f1-4ae2-9b9e-45bf5c2bf854.root", "root://cmsdcache-kit-disk.gridka.de:1094//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/a5519acd-97bd-42ef-9ca1-1fb4446aa75f.root"]}' --output_format root