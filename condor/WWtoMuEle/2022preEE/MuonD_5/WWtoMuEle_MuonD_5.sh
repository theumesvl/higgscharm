#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022preEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022preEE/MuonD --dataset MuonD_5 --partition_fileset '{"MuonD_5": ["root://gaexrdoor.ciemat.es:1094//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/6cbbe7a8-ec03-46f1-acbe-7754f0882892.root", "root://eos.cms.rcac.purdue.edu//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/7859ce73-1c25-4cef-8d4f-7610a00d2e55.root", "root://grid-cms-xrootd.physik.rwth-aachen.de:1094//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/78a50190-af2f-4d52-b4f0-9735eea6f374.root", "root://cceos.ihep.ac.cn:1094//eos/ihep/cms/store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/79bb47d4-fe42-4fa1-b632-5c8a6b4aaa24.root", "root://redirector.t2.ucsd.edu:1095//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/7dd31587-3cbb-4bac-8b2f-ff3441d19952.root", "root://grid-cms-xrootd.physik.rwth-aachen.de:1094//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/829cc24b-3c15-4ff5-af7a-92afe651b0a4.root", "root://cmsdcache-kit-disk.gridka.de:1094//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/8359be4a-25cb-4e66-89d5-13d8436e8ef1.root", "root://maite.iihe.ac.be:1095//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/8b042d18-fae3-4875-9eb9-42e33df950d5.root", "root://xrootd-local.unl.edu:1094//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/8c579dae-c362-4c09-b484-cf2b9b8138a3.root"]}' --output_format root