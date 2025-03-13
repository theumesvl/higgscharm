#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022postEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022postEE/MuonG --dataset MuonG_8 --partition_fileset '{"MuonG_8": ["root://se.cis.gov.pl:1094//grid/cms/store/data/Run2022G/Muon/NANOAOD/22Sep2023-v1/30000/da658f76-629e-4b79-b322-9e17a7b8a72f.root", "root://eos.cms.rcac.purdue.edu//store/data/Run2022G/Muon/NANOAOD/22Sep2023-v1/30000/dc344ec4-47d9-4bec-8480-a93b346837ba.root", "root://grid-cms-xrootd.physik.rwth-aachen.de:1094//store/data/Run2022G/Muon/NANOAOD/22Sep2023-v1/30000/dcaf599f-7c85-4674-9ea2-1ce8a209dd1d.root", "root://k8s-redir.ultralight.org:1094//store/data/Run2022G/Muon/NANOAOD/22Sep2023-v1/30000/e501358e-9abf-4e5e-a22b-eb41ed8af32e.root", "root://gaexrdoor.ciemat.es:1094//store/data/Run2022G/Muon/NANOAOD/22Sep2023-v1/30000/e565e554-c7d0-4192-a045-19aefdf23b95.root", "root://eos.cms.rcac.purdue.edu//store/data/Run2022G/Muon/NANOAOD/22Sep2023-v1/30000/e81537d4-611d-48a1-bd56-8a180ee282b5.root", "root://maite.iihe.ac.be:1095//store/data/Run2022G/Muon/NANOAOD/22Sep2023-v1/30000/e8fc87af-b434-4019-ad2b-064e38f82393.root", "root://eos.cms.rcac.purdue.edu//store/data/Run2022G/Muon/NANOAOD/22Sep2023-v1/30000/ed07e2ee-5864-4cf8-a7d3-11ffbcc4d5fa.root", "root://hactar01.crc.nd.edu//store/data/Run2022G/Muon/NANOAOD/22Sep2023-v1/30000/ede81d67-0462-4f4d-91a1-865e3a4427a4.root"]}' --output_format root