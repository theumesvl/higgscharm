#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022postEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022postEE/EGammaE --dataset EGammaE_18 --partition_fileset '{"EGammaE_18": ["root://maite.iihe.ac.be:1095//store/data/Run2022E/EGamma/NANOAOD/22Sep2023-v1/2530000/f9b3c6ae-aea0-481f-a592-26013b06950b.root", "root://redirector.t2.ucsd.edu:1095//store/data/Run2022E/EGamma/NANOAOD/22Sep2023-v1/2530000/f9e5b750-57a9-4a3b-9241-1172fbc990bd.root", "root://storage01.lcg.cscs.ch:1096//pnfs/lcg.cscs.ch/cms/trivcat/store/data/Run2022E/EGamma/NANOAOD/22Sep2023-v1/2530000/fa2ba6df-eac8-4ba0-8bca-377f68f494a8.root", "root://eos.cms.rcac.purdue.edu//store/data/Run2022E/EGamma/NANOAOD/22Sep2023-v1/2530000/fa746ab6-41b6-439e-b632-5d8f0eff101a.root", "root://xrootd-cms.infn.it:1194//store/data/Run2022E/EGamma/NANOAOD/22Sep2023-v1/2530000/fabf483b-fb8e-40a9-9db6-ef764c09e87d.root", "root://grid-cms-xrootd.physik.rwth-aachen.de:1094//store/data/Run2022E/EGamma/NANOAOD/22Sep2023-v1/2530000/fac4a564-84ef-43ff-85bf-f3f03cd3e493.root", "root://maite.iihe.ac.be:1095//store/data/Run2022E/EGamma/NANOAOD/22Sep2023-v1/2530000/fbd052c8-f0f0-4fb0-ac13-6aa1aa817265.root", "root://xrootd-cms.infn.it:1194//store/data/Run2022E/EGamma/NANOAOD/22Sep2023-v1/2530000/fdbb4bbb-c5d5-4ecf-857d-f5e267c0945f.root", "root://eos.cms.rcac.purdue.edu//store/data/Run2022E/EGamma/NANOAOD/22Sep2023-v1/2530000/fe44fa23-b0a5-4c5d-a515-d56e7d101ac5.root"]}' --output_format root