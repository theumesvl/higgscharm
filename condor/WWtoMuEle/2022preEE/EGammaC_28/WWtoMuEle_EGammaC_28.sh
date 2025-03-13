#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022preEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022preEE/EGammaC --dataset EGammaC_28 --partition_fileset '{"EGammaC_28": ["root://xrootd-vanderbilt.sites.opensciencegrid.org:1094//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/40000/f13df640-b2bb-455a-b12e-591bdcf478c9.root", "root://cmsio2.rc.ufl.edu:1094//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/40000/f23c8368-51bd-4cb6-83b1-71ac89d4b878.root", "root://maite.iihe.ac.be:1095//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/40000/f2c2d092-dab9-4fcb-b50d-30611bc2e9b0.root", "root://maite.iihe.ac.be:1095//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/40000/f3c4420c-dc5a-4d56-8f05-29174ff1b924.root", "root://dcache-cms-xrootd.desy.de:1094//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/40000/f5c6ef9a-a623-44be-ad5b-bd537972572d.root", "root://cmsio2.rc.ufl.edu:1094//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/40000/f6037a2e-3254-45c8-88cd-23d21e2505a8.root", "root://dcache-cms-xrootd.desy.de:1094//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/40000/f9e9872e-435c-480f-b0f7-06b3427d865e.root", "root://gaexrdoor.ciemat.es:1094//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/40000/fd0c6ff9-0e84-4b0c-8834-2886941b3740.root", "root://xrootd-vanderbilt.sites.opensciencegrid.org:1094//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/40000/fdfef697-4d08-47d9-adf4-96a5b57d53ba.root"]}' --output_format root