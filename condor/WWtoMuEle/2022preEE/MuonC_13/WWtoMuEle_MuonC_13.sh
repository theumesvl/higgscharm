#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022preEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022preEE/MuonC --dataset MuonC_13 --partition_fileset '{"MuonC_13": ["root://hactar01.crc.nd.edu//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/50000/8c96e87f-466d-4a63-82da-7359aa2bd965.root", "root://gaexrdoor.ciemat.es:1094//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/50000/a14e16d1-3a38-4077-8941-193f614a9283.root", "root://gaexrdoor.ciemat.es:1094//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/50000/a1fb9be8-7b21-4ee4-9ce0-7e1a911d8b5a.root", "root://hactar01.crc.nd.edu//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/50000/b8bf1e1d-2a37-498c-962a-64b1fbb4ec0f.root", "root://dcache-cms-xrootd.desy.de:1094//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/50000/c1b07ced-0d21-4214-b15d-4f93c2295dec.root", "root://hactar01.crc.nd.edu//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/50000/d2a83da6-d422-4f4b-ad96-6f3ab6b913ab.root", "root://xrootd-local.unl.edu:1094//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/50000/e82012d4-ce19-4bff-98f2-a9a651d0efe6.root", "root://hactar01.crc.nd.edu//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/50000/edb51905-88b7-4eed-aecb-5b29ce7ca96d.root", "root://cceos.ihep.ac.cn:1094//eos/ihep/cms/store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/50000/ffad8310-8bee-4828-97ce-7938c2d7779d.root"]}' --output_format root