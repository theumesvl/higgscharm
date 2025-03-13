#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022preEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022preEE/MuonD --dataset MuonD_9 --partition_fileset '{"MuonD_9": ["root://cceos.ihep.ac.cn:1094//eos/ihep/cms/store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/f07fe947-ef48-41ba-9c68-096be4dbf4fb.root", "root://xrootd-local.unl.edu:1094//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/f1f5bf4b-76e2-45b9-8fac-29df7ff82984.root", "root://cmsio2.rc.ufl.edu:1094//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/f421d262-2043-4600-a989-f105f5b64a36.root", "root://eos.cms.rcac.purdue.edu//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/f47cdfd4-89a5-44eb-9b29-6ddb3786d740.root", "root://cmsio2.rc.ufl.edu:1094//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/f4dcb9cc-dfee-4770-9484-6b3b49616cc2.root", "root://k8s-redir.ultralight.org:1094//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/f8f86b58-a285-46b9-86a9-3f8d4e7bfc12.root", "root://redirector.t2.ucsd.edu:1095//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/fa88ef2d-d755-4de5-b90f-64c5dacbc5d2.root", "root://gaexrdoor.ciemat.es:1094//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/fe6cc579-3e6c-4981-835e-483e7df686c9.root", "root://k8s-redir.ultralight.org:1094//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/fe7df7e9-614e-4c7e-9845-eb285df0498e.root"]}' --output_format root