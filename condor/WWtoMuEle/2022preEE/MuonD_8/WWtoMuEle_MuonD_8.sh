#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022preEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022preEE/MuonD --dataset MuonD_8 --partition_fileset '{"MuonD_8": ["root://xrootd-vanderbilt.sites.opensciencegrid.org:1094//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/d7495d9d-4a1e-4f18-b7eb-b56c06d0bbdd.root", "root://eos.grid.vbc.ac.at:1094//eos/vbc/experiments/cms/store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/d861a2c2-fd5c-4616-be42-16c1c8907a4f.root", "root://k8s-redir.ultralight.org:1094//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/de9c4285-a4e8-4fd9-af02-4e1735612ea1.root", "root://xrootd-vanderbilt.sites.opensciencegrid.org:1094//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/e511199a-f4b2-4e54-b0d2-66f660bd5ece.root", "root://k8s-redir.ultralight.org:1094//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/e52cefee-6876-4fcf-885c-4fc2507a9027.root", "root://redirector.t2.ucsd.edu:1095//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/e6969546-aa16-490b-8375-5064885dee1d.root", "root://eos.cms.rcac.purdue.edu//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/e779a65b-67f6-4048-8eb6-3029518e558c.root", "root://eos.grid.vbc.ac.at:1094//eos/vbc/experiments/cms/store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/e79123e9-148a-423f-8fe8-47fa351714cb.root", "root://maite.iihe.ac.be:1095//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/ea3daa90-8674-4909-9f17-f0375f711bb0.root"]}' --output_format root