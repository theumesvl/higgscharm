#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022preEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022preEE/EGammaC --dataset EGammaC_30 --partition_fileset '{"EGammaC_30": ["root://grid143.kfki.hu:1094//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/50000/551681da-d923-49a2-9bd7-0803fdbcec96.root", "root://xrootd-vanderbilt.sites.opensciencegrid.org:1094//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/50000/64624550-d985-4dce-9a65-84d5c9f177b7.root", "root://cmsdcadisk.fnal.gov//dcache/uscmsdisk/store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/50000/66644155-9330-4757-962c-a8608f1d0ce3.root", "root://cmsdcadisk.fnal.gov//dcache/uscmsdisk/store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/50000/6acb7c3c-fd21-41e0-9d2b-c2aa7163395c.root", "root://storage01.lcg.cscs.ch:1096//pnfs/lcg.cscs.ch/cms/trivcat/store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/50000/7769ec43-5822-4efe-8b4e-c04ca3348bd4.root", "root://grid143.kfki.hu:1094//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/50000/7c150ab1-2683-4163-8aac-94e5f8b7463e.root", "root://xrootd-vanderbilt.sites.opensciencegrid.org:1094//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/50000/7cae62d1-6346-404f-8b37-7232479b55f1.root", "root://dcache-cms-xrootd.desy.de:1094//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/50000/7dfe7893-2219-48e8-a1d9-67f361d1adc2.root", "root://cmsio2.rc.ufl.edu:1094//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/50000/7fcc3554-d3f7-4919-b971-6eb21bd862f7.root"]}' --output_format root