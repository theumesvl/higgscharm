#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022preEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022preEE/MuonC --dataset MuonC_8 --partition_fileset '{"MuonC_8": ["root://hactar01.crc.nd.edu//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/30000/9fdf3b62-4eda-4fc9-a2f6-ad0363672f62.root", "root://k8s-redir.ultralight.org:1094//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/30000/a1970ec7-4fcc-49b5-9369-a83889426a1f.root", "root://cceos.ihep.ac.cn:1094//eos/ihep/cms/store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/30000/a3241ae6-4c24-4d68-b148-55d8b7f190f0.root", "root://xrootd-cms.infn.it:1194//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/30000/aae17dfb-8d6a-4f79-ae2a-aec041064529.root", "root://cceos.ihep.ac.cn:1094//eos/ihep/cms/store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/30000/b0a27791-b3af-4347-ae90-4f88faecb833.root", "root://grid-cms-xrootd.physik.rwth-aachen.de:1094//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/30000/b4af60fc-c19f-482a-b99d-f20243810efe.root", "root://k8s-redir.ultralight.org:1094//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/30000/b599bd66-0aea-4bdb-96ba-2b0ace1f0f15.root", "root://dcache-cms-xrootd.desy.de:1094//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/30000/bae5ab41-ee30-4578-870b-33b77fb5f1f6.root", "root://grid143.kfki.hu:1094//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/30000/bbaac590-d287-477d-9a05-28c6b8e84e9c.root"]}' --output_format root