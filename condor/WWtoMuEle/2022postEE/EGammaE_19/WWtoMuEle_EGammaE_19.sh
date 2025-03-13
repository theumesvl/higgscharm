#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022postEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022postEE/EGammaE --dataset EGammaE_19 --partition_fileset '{"EGammaE_19": ["root://gaexrdoor.ciemat.es:1094//store/data/Run2022E/EGamma/NANOAOD/22Sep2023-v1/2530000/ff8385b4-a4f8-4638-a57d-85f5ef39a66d.root", "root://gaexrdoor.ciemat.es:1094//store/data/Run2022E/EGamma/NANOAOD/22Sep2023-v1/40000/2c8e1cfa-d8f3-4383-bfea-1791fcc2857a.root", "root://cmsio2.rc.ufl.edu:1094//store/data/Run2022E/EGamma/NANOAOD/22Sep2023-v1/40000/35a67222-1476-4c9c-9d84-c31e580a5a7b.root", "root://ruhex-osgce.rutgers.edu//store/data/Run2022E/EGamma/NANOAOD/22Sep2023-v1/40000/4939a46d-6c24-4d4e-9f8c-c4c147ffc70a.root", "root://xrootd-cms.infn.it:1194//store/data/Run2022E/EGamma/NANOAOD/22Sep2023-v1/40000/53b0deda-473b-41cb-b8fd-c90ba95dd392.root", "root://hactar01.crc.nd.edu//store/data/Run2022E/EGamma/NANOAOD/22Sep2023-v1/40000/7c817ff3-b12a-43fe-8e92-c49c94e376f0.root", "root://gaexrdoor.ciemat.es:1094//store/data/Run2022E/EGamma/NANOAOD/22Sep2023-v1/40000/860df1c7-aaa4-42e2-80ff-7a19b5ca8fd3.root", "root://cmsxrootd.hep.wisc.edu:1094//store/data/Run2022E/EGamma/NANOAOD/22Sep2023-v1/50000/1ffc3233-81fd-4bc4-8525-40167aa5a6da.root", "root://storage01.lcg.cscs.ch:1096//pnfs/lcg.cscs.ch/cms/trivcat/store/data/Run2022E/EGamma/NANOAOD/22Sep2023-v1/50000/a2c2a0c9-b83e-4e96-b9e0-4ee52b136d86.root"]}' --output_format root