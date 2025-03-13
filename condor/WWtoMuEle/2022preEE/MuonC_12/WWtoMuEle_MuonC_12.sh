#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022preEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022preEE/MuonC --dataset MuonC_12 --partition_fileset '{"MuonC_12": ["root://cmsio2.rc.ufl.edu:1094//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/30000/fc884714-2086-4c1c-94a8-0ed100ecdd7e.root", "root://xrootd-cms.infn.it:1194//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/30000/fe3438fb-36fd-4ca0-a385-64cc73082d0c.root", "root://dcache-cms-xrootd.desy.de:1094//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/50000/0246d95b-cea7-4dcf-8c35-e132c1310b66.root", "root://cmsdcache-kit-disk.gridka.de:1094//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/50000/0a8a4c73-882a-468d-b29c-b2753289dbee.root", "root://eos.cms.rcac.purdue.edu//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/50000/12b9a2ba-c5b6-40fc-8455-d111d391bdd5.root", "root://grid143.kfki.hu:1094//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/50000/18f62f2a-c99e-452c-ac3f-348655a2e4d8.root", "root://grid143.kfki.hu:1094//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/50000/53995cbd-0258-4ec6-9873-e3ef184ac1c7.root", "root://k8s-redir.ultralight.org:1094//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/50000/6239896a-a85c-443a-88dd-bdd7c9b98349.root", "root://maite.iihe.ac.be:1095//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/50000/889f30d4-4412-4f3b-b0cc-c67b9a26168b.root"]}' --output_format root