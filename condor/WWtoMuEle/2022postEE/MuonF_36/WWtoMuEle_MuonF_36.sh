#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022postEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022postEE/MuonF --dataset MuonF_36 --partition_fileset '{"MuonF_36": ["root://ruhex-osgce.rutgers.edu//store/data/Run2022F/Muon/NANOAOD/22Sep2023-v2/70000/febd6929-9bda-43cd-8c28-78274f460ec3.root", "root://redirector.t2.ucsd.edu:1095//store/data/Run2022F/Muon/NANOAOD/22Sep2023-v2/80000/011ae799-69e3-47e0-a9c4-5d0b8e4b21e8.root", "root://grid-cms-xrootd.physik.rwth-aachen.de:1094//store/data/Run2022F/Muon/NANOAOD/22Sep2023-v2/80000/18ff64e8-c210-4448-b854-f118480c3bd2.root", "root://gaexrdoor.ciemat.es:1094//store/data/Run2022F/Muon/NANOAOD/22Sep2023-v2/80000/203fc2ec-7893-499a-ac3e-6233113bf499.root", "root://k8s-redir.ultralight.org:1094//store/data/Run2022F/Muon/NANOAOD/22Sep2023-v2/80000/24c1f525-c162-4fd9-8c5b-ce8dc1d4ff51.root", "root://cceos.ihep.ac.cn:1094//eos/ihep/cms/store/data/Run2022F/Muon/NANOAOD/22Sep2023-v2/80000/5480a042-f871-42db-85f2-b6303a18cc8a.root", "root://maite.iihe.ac.be:1095//store/data/Run2022F/Muon/NANOAOD/22Sep2023-v2/80000/5615d7ea-1b24-4370-80c0-30c6f907514c.root", "root://ruhex-osgce.rutgers.edu//store/data/Run2022F/Muon/NANOAOD/22Sep2023-v2/80000/e3168142-667a-4bd4-b8b6-f160368832a1.root", "root://grid-cms-xrootd.physik.rwth-aachen.de:1094//store/data/Run2022F/Muon/NANOAOD/22Sep2023-v2/80000/f28e42ee-9ce7-4554-883c-8dba709a2932.root"]}' --output_format root