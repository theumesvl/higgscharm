#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022preEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022preEE/MuonC --dataset MuonC_9 --partition_fileset '{"MuonC_9": ["root://hactar01.crc.nd.edu//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/30000/bc76b388-7339-426e-a07a-e4e185875ed4.root", "root://redirector.t2.ucsd.edu:1095//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/30000/be0088e7-303e-496e-98e7-1c98292b645c.root", "root://hactar01.crc.nd.edu//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/30000/be907269-d744-4dc4-b916-f8591a36978a.root", "root://cceos.ihep.ac.cn:1094//eos/ihep/cms/store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/30000/c0dbc6c1-b33e-491f-81b0-79dc7b72462e.root", "root://hactar01.crc.nd.edu//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/30000/c2478a88-177b-467b-947b-44bbd0ddc628.root", "root://xrootd-local.unl.edu:1094//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/30000/c6de94e3-4a64-41ee-88e2-b4de7bb77ce9.root", "root://eos.cms.rcac.purdue.edu//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/30000/c780eeaf-eb0c-4086-a269-c08479202354.root", "root://gaexrdoor.ciemat.es:1094//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/30000/c7ea9680-3526-4dd6-a54a-0950a8aa396d.root", "root://grid143.kfki.hu:1094//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/30000/c8f2766f-e6ba-4295-9477-642c35859bcf.root"]}' --output_format root