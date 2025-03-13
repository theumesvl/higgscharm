#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022preEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022preEE/EGammaD --dataset EGammaD_11 --partition_fileset '{"EGammaD_11": ["root://osg-se.sprace.org.br:1094//store/data/Run2022D/EGamma/NANOAOD/22Sep2023-v1/30000/56ee9e9c-31a7-4c7c-b22a-ff006cda7500.root", "root://cmsdcadisk.fnal.gov//dcache/uscmsdisk/store/data/Run2022D/EGamma/NANOAOD/22Sep2023-v1/30000/74a041b0-968d-4c82-bc6e-20232744d5f0.root", "root://hactar01.crc.nd.edu//store/data/Run2022D/EGamma/NANOAOD/22Sep2023-v1/30000/9ae84f95-1154-49cc-90ff-5d10dc6303c7.root", "root://storage01.lcg.cscs.ch:1096//pnfs/lcg.cscs.ch/cms/trivcat/store/data/Run2022D/EGamma/NANOAOD/22Sep2023-v1/30000/9df50dcd-f7ea-4c58-9073-b0f1a294e1ff.root", "root://eos.cms.rcac.purdue.edu//store/data/Run2022D/EGamma/NANOAOD/22Sep2023-v1/30000/c7eea1da-0974-44e3-bf33-b40c14859fb9.root", "root://eos.cms.rcac.purdue.edu//store/data/Run2022D/EGamma/NANOAOD/22Sep2023-v1/30000/cb8d11cd-a69b-451e-ab12-fd8a94199385.root", "root://osg-se.sprace.org.br:1094//store/data/Run2022D/EGamma/NANOAOD/22Sep2023-v1/30000/df29deb4-e128-4807-8d2f-06ee07fd7f19.root", "root://hactar01.crc.nd.edu//store/data/Run2022D/EGamma/NANOAOD/22Sep2023-v1/30000/eea71cb8-190c-4ad6-9bf0-3d7ae512a764.root", "root://maite.iihe.ac.be:1095//store/data/Run2022D/EGamma/NANOAOD/22Sep2023-v1/410000/8374e6b0-c97e-4bc8-b095-9eaee52632e6.root"]}' --output_format root