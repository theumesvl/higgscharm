#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022postEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022postEE/MuonEGG --dataset MuonEGG_1 --partition_fileset '{"MuonEGG_1": ["root://ruhex-osgce.rutgers.edu//store/data/Run2022G/MuonEG/NANOAOD/22Sep2023-v1/2520000/1c9153d0-4fec-4e52-8aae-95ecb7909213.root", "root://hactar01.crc.nd.edu//store/data/Run2022G/MuonEG/NANOAOD/22Sep2023-v1/2520000/201190bd-157c-46f7-bfc1-7b7fa7cdd6fb.root", "root://k8s-redir.ultralight.org:1094//store/data/Run2022G/MuonEG/NANOAOD/22Sep2023-v1/2520000/2a297070-04a7-48f2-ab77-2baedf21c6a1.root", "root://ruhex-osgce.rutgers.edu//store/data/Run2022G/MuonEG/NANOAOD/22Sep2023-v1/2520000/30f0e075-1a2d-43d9-9dd8-7de2ff256a07.root", "root://eos01.grid.cyfronet.pl:1094//eos/cms/store/data/Run2022G/MuonEG/NANOAOD/22Sep2023-v1/2520000/3377ea06-2b74-48cd-a04c-e897dc7743b4.root", "root://k8s-redir.ultralight.org:1094//store/data/Run2022G/MuonEG/NANOAOD/22Sep2023-v1/2520000/405b00e8-13f0-4142-a58a-7011cc10d785.root", "root://maite.iihe.ac.be:1095//store/data/Run2022G/MuonEG/NANOAOD/22Sep2023-v1/2520000/40a3e92a-ec2f-4781-989e-c129775bcb78.root", "root://hactar01.crc.nd.edu//store/data/Run2022G/MuonEG/NANOAOD/22Sep2023-v1/2520000/4cef275b-6f48-432b-af1c-d38f4053f4ef.root", "root://maite.iihe.ac.be:1095//store/data/Run2022G/MuonEG/NANOAOD/22Sep2023-v1/2520000/4edad7fb-956b-4417-8493-5c5bd2335176.root"]}' --output_format root