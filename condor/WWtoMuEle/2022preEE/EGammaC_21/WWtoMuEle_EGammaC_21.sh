#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022preEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022preEE/EGammaC --dataset EGammaC_21 --partition_fileset '{"EGammaC_21": ["root://hactar01.crc.nd.edu//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/40000/8a1f94bc-d623-4107-a54f-bc8f7532a757.root", "root://maite.iihe.ac.be:1095//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/40000/8c41e42e-7e98-4eb5-8f80-9a674b86bfdc.root", "root://eoscms.cern.ch//eos/cms/store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/40000/8db8f401-388e-47e3-96e5-1de89aee74da.root", "root://t3se01.psi.ch:1094//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/40000/90c94494-d854-41dc-8ef3-829950fca9c5.root", "root://hactar01.crc.nd.edu//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/40000/90f667e4-6eb4-4a4d-ac1f-026deb15d109.root", "root://eoscms.cern.ch//eos/cms/store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/40000/9327925d-b457-4af2-908c-0c2633b4c2ff.root", "root://t3se01.psi.ch:1094//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/40000/976c5678-a9e8-4645-a86f-b48631eb7612.root", "root://grid143.kfki.hu:1094//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/40000/995e6556-488b-4c78-aa5c-ebe61964d265.root", "root://gaexrdoor.ciemat.es:1094//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/40000/9a8535b5-0ec1-44b9-a5cd-8f210404e3af.root", "root://t3se01.psi.ch:1094//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/40000/9b4a8982-fae2-4824-ab20-b954aa89400c.root"]}' --output_format root