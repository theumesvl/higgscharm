#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022postEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022postEE/EGammaG --dataset EGammaG_8 --partition_fileset '{"EGammaG_8": ["root://xrootd.hep.kbfi.ee:1094//store/data/Run2022G/EGamma/NANOAOD/22Sep2023-v2/2560000/b5cff409-2e56-4607-8888-3ec767d38517.root", "root://cmsio2.rc.ufl.edu:1094//store/data/Run2022G/EGamma/NANOAOD/22Sep2023-v2/2560000/cc53410c-06d3-4713-aadf-fbec7bb144ac.root", "root://hactar01.crc.nd.edu//store/data/Run2022G/EGamma/NANOAOD/22Sep2023-v2/2560000/ce59f39b-fcf7-4561-8157-8541280debe4.root", "root://xrootd.hep.kbfi.ee:1094//store/data/Run2022G/EGamma/NANOAOD/22Sep2023-v2/2560000/cec93220-d45c-4e7d-bc84-43cc1092efc5.root", "root://storage01.lcg.cscs.ch:1096//pnfs/lcg.cscs.ch/cms/trivcat/store/data/Run2022G/EGamma/NANOAOD/22Sep2023-v2/2560000/dd135067-4f60-4909-bbd4-a783e548ab84.root", "root://xrootd.hep.kbfi.ee:1094//store/data/Run2022G/EGamma/NANOAOD/22Sep2023-v2/2560000/e0386f63-1045-4326-a595-daba9d76d7c5.root", "root://cmsio2.rc.ufl.edu:1094//store/data/Run2022G/EGamma/NANOAOD/22Sep2023-v2/2560000/e0b317f3-3ca9-4823-93c7-d0908af14c14.root", "root://hactar01.crc.nd.edu//store/data/Run2022G/EGamma/NANOAOD/22Sep2023-v2/2560000/e263c79a-1020-4133-9a0e-b12551917359.root", "root://xrootd.hep.kbfi.ee:1094//store/data/Run2022G/EGamma/NANOAOD/22Sep2023-v2/2560000/ebd1518a-972f-4561-b722-01d69ae2219e.root"]}' --output_format root