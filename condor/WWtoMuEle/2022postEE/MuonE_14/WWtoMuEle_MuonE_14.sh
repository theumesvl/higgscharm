#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022postEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022postEE/MuonE --dataset MuonE_14 --partition_fileset '{"MuonE_14": ["root://k8s-redir.ultralight.org:1094//store/data/Run2022E/Muon/NANOAOD/22Sep2023-v1/30000/6e6647fa-291c-4a4e-960d-898860dfb90f.root", "root://k8s-redir.ultralight.org:1094//store/data/Run2022E/Muon/NANOAOD/22Sep2023-v1/30000/743aeff5-3954-449f-8e85-8057e8c436e1.root", "root://cceos.ihep.ac.cn:1094//eos/ihep/cms/store/data/Run2022E/Muon/NANOAOD/22Sep2023-v1/30000/7d4749dc-b8f2-4403-a5bc-f3b7807fc210.root", "root://cceos.ihep.ac.cn:1094//eos/ihep/cms/store/data/Run2022E/Muon/NANOAOD/22Sep2023-v1/30000/80134792-8121-4ea0-8976-afc779260384.root", "root://hactar01.crc.nd.edu//store/data/Run2022E/Muon/NANOAOD/22Sep2023-v1/30000/8050ac31-9f92-41e8-b53d-276791560806.root", "root://cmsio2.rc.ufl.edu:1094//store/data/Run2022E/Muon/NANOAOD/22Sep2023-v1/30000/89603db4-17ee-47cf-9ff1-a4c4c62904ad.root", "root://cceos.ihep.ac.cn:1094//eos/ihep/cms/store/data/Run2022E/Muon/NANOAOD/22Sep2023-v1/30000/9f310c3b-0f3d-4dc8-a34f-ae37da350ed8.root", "root://cceos.ihep.ac.cn:1094//eos/ihep/cms/store/data/Run2022E/Muon/NANOAOD/22Sep2023-v1/30000/a41677d2-36ae-422a-9d74-18a78694d538.root", "root://ruhex-osgce.rutgers.edu//store/data/Run2022E/Muon/NANOAOD/22Sep2023-v1/30000/aa6cb208-83b2-46bf-8e26-3dfd7714eebd.root"]}' --output_format root