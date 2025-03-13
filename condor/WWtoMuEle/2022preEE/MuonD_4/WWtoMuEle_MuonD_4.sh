#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022preEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022preEE/MuonD --dataset MuonD_4 --partition_fileset '{"MuonD_4": ["root://eos.grid.vbc.ac.at:1094//eos/vbc/experiments/cms/store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/523d5db2-613a-45de-818c-bf8846cd93fb.root", "root://k8s-redir.ultralight.org:1094//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/5784d950-aba4-4922-b4e6-1ad6eed4045b.root", "root://gaexrdoor.ciemat.es:1094//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/5ab99ef1-6ce7-4cd1-abfd-2fa701dd3105.root", "root://hactar01.crc.nd.edu//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/5aba2972-601b-452f-a68e-3f27f5a7b1d0.root", "root://dcache-cms-xrootd.desy.de:1094//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/6133bf22-ddc0-4f61-b4b6-a41b5fb6f4d5.root", "root://eos.cms.rcac.purdue.edu//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/63a2727a-2594-4128-b426-64517693fb5c.root", "root://grid-cms-xrootd.physik.rwth-aachen.de:1094//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/659ceecb-e41a-4ec3-9899-f5620975edcd.root", "root://xrootd-local.unl.edu:1094//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/66151aa6-9eeb-4df5-850b-661127cb12f5.root", "root://hactar01.crc.nd.edu//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/6b687e02-a99d-422f-a18b-bd0c6e715181.root"]}' --output_format root