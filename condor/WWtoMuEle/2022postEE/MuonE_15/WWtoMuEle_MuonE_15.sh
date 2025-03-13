#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022postEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022postEE/MuonE --dataset MuonE_15 --partition_fileset '{"MuonE_15": ["root://hactar01.crc.nd.edu//store/data/Run2022E/Muon/NANOAOD/22Sep2023-v1/30000/c4e8fe81-c8e1-4d56-8a85-674b27c5fd44.root", "root://redirector.t2.ucsd.edu:1095//store/data/Run2022E/Muon/NANOAOD/22Sep2023-v1/30000/e10f5c92-287d-44ce-aa29-720692d7fc75.root", "root://redirector.t2.ucsd.edu:1095//store/data/Run2022E/Muon/NANOAOD/22Sep2023-v1/30000/e54e14e3-d9a3-4c63-b587-1b090a349d33.root", "root://k8s-redir.ultralight.org:1094//store/data/Run2022E/Muon/NANOAOD/22Sep2023-v1/30000/ed9bf0f9-40d7-4d00-b727-cf40a9cf36c8.root", "root://maite.iihe.ac.be:1095//store/data/Run2022E/Muon/NANOAOD/22Sep2023-v1/30000/f2e4991a-f4f6-4ea3-b061-6d2245708572.root", "root://hactar01.crc.nd.edu//store/data/Run2022E/Muon/NANOAOD/22Sep2023-v1/40000/18134f06-6928-4b90-8218-11f0853ef857.root", "root://ruhex-osgce.rutgers.edu//store/data/Run2022E/Muon/NANOAOD/22Sep2023-v1/40000/19ad937b-86d2-459b-9495-fd86792727f2.root", "root://gaexrdoor.ciemat.es:1094//store/data/Run2022E/Muon/NANOAOD/22Sep2023-v1/40000/1cd43e1a-4750-44bd-9d30-e5b90a3b576e.root", "root://k8s-redir.ultralight.org:1094//store/data/Run2022E/Muon/NANOAOD/22Sep2023-v1/40000/cb3afe1f-699f-4be6-a1ed-dd6fd0b2ca81.root"]}' --output_format root