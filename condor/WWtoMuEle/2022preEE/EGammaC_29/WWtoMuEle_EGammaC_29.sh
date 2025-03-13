#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022preEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022preEE/EGammaC --dataset EGammaC_29 --partition_fileset '{"EGammaC_29": ["root://dcache-cms-xrootd.desy.de:1094//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/40000/fe1fee55-8d10-41ad-9e8e-50063e1640c0.root", "root://hactar01.crc.nd.edu//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/50000/1c6b3f1e-c8cd-43ef-83d0-c5854bb3c2cf.root", "root://redirector.t2.ucsd.edu:1095//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/50000/1f10c0a2-7770-471d-be46-6214738f34a2.root", "root://maite.iihe.ac.be:1095//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/50000/26ffd3fe-81f5-450d-9bd9-239d7234269e.root", "root://grid-cms-xrootd.physik.rwth-aachen.de:1094//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/50000/344a5e80-338c-472d-bdf4-933e66f9b658.root", "root://hactar01.crc.nd.edu//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/50000/38df5d5b-053a-48ae-aa7a-432e45a9fb8b.root", "root://eoscms.cern.ch//eos/cms/store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/50000/3b585eaa-b3ba-4bf6-abf4-fe424fa5c1e2.root", "root://maite.iihe.ac.be:1095//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/50000/3c15f464-10bd-4ab8-9064-a509a3470c8a.root", "root://dcache-cms-xrootd.desy.de:1094//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/50000/42521add-a1f8-45bb-9428-1bfaa9fb9552.root"]}' --output_format root