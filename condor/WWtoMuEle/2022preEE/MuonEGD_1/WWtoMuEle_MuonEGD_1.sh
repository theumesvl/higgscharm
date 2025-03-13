#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022preEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022preEE/MuonEGD --dataset MuonEGD_1 --partition_fileset '{"MuonEGD_1": ["root://cmsdcache-kit-disk.gridka.de:1094//store/data/Run2022D/MuonEG/NANOAOD/22Sep2023-v1/40000/27453cd2-36d5-4b34-9cf6-1303480b5dcf.root", "root://eos.cms.rcac.purdue.edu//store/data/Run2022D/MuonEG/NANOAOD/22Sep2023-v1/40000/3af87944-1b7d-408a-a809-806b00a7d561.root", "root://eos.cms.rcac.purdue.edu//store/data/Run2022D/MuonEG/NANOAOD/22Sep2023-v1/40000/3bf393d1-985b-45ab-863d-589b05b6e50c.root", "root://cmsdcache-kit-disk.gridka.de:1094//store/data/Run2022D/MuonEG/NANOAOD/22Sep2023-v1/40000/405c5f48-2581-4aa3-a547-ea2a690a56b4.root", "root://hactar01.crc.nd.edu//store/data/Run2022D/MuonEG/NANOAOD/22Sep2023-v1/40000/42803125-9b57-417b-bba1-4217bb983489.root", "root://eos.cms.rcac.purdue.edu//store/data/Run2022D/MuonEG/NANOAOD/22Sep2023-v1/40000/44183ea8-8843-455c-8cb3-3831a03106e8.root", "root://maite.iihe.ac.be:1095//store/data/Run2022D/MuonEG/NANOAOD/22Sep2023-v1/40000/685d399d-2a05-42c1-bb0f-119951424112.root", "root://grid143.kfki.hu:1094//store/data/Run2022D/MuonEG/NANOAOD/22Sep2023-v1/40000/6a93dde7-3b91-4bc4-93db-b59fc2c3f6b6.root"]}' --output_format root