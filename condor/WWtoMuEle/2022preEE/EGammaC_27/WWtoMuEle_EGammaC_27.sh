#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022preEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022preEE/EGammaC --dataset EGammaC_27 --partition_fileset '{"EGammaC_27": ["root://gaexrdoor.ciemat.es:1094//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/40000/e574d414-b22a-435f-863c-ba4a04bd0d32.root", "root://eos.cms.rcac.purdue.edu//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/40000/e58b370f-495c-49cb-a1c7-149541a005cf.root", "root://redirector.t2.ucsd.edu:1095//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/40000/e5a34237-6488-4b62-8f2c-34348335b076.root", "root://grid143.kfki.hu:1094//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/40000/e7d9bb3b-a424-4319-a62b-207f23b188b9.root", "root://eoscms.cern.ch//eos/cms/store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/40000/eb3c86c8-053f-4e37-b296-36bc3edab4fc.root", "root://hactar01.crc.nd.edu//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/40000/ebde7b7f-8855-41fb-ac15-d223c8451d56.root", "root://cmsdcadisk.fnal.gov//dcache/uscmsdisk/store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/40000/edb5e9ed-cda9-4d8a-aef0-641e3080e4f8.root", "root://xrootd-vanderbilt.sites.opensciencegrid.org:1094//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/40000/edc968c5-c8ae-4646-abde-94797f8f37bc.root", "root://eoscms.cern.ch//eos/cms/store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/40000/f1321844-176c-49cc-83c8-215c7507c636.root"]}' --output_format root