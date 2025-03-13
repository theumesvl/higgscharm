#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022postEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022postEE/EGammaG --dataset EGammaG_7 --partition_fileset '{"EGammaG_7": ["root://t3se01.psi.ch:1094//store/data/Run2022G/EGamma/NANOAOD/22Sep2023-v2/2560000/6dec84e1-445e-418d-9fff-aef17411271b.root", "root://cmsio2.rc.ufl.edu:1094//store/data/Run2022G/EGamma/NANOAOD/22Sep2023-v2/2560000/7d21927f-4ae4-4fc9-9b9c-be922dc3d507.root", "root://maite.iihe.ac.be:1095//store/data/Run2022G/EGamma/NANOAOD/22Sep2023-v2/2560000/9054ff39-9c15-4f5a-9a68-1acc2d4cad52.root", "root://maite.iihe.ac.be:1095//store/data/Run2022G/EGamma/NANOAOD/22Sep2023-v2/2560000/90fbe00a-0029-41c9-b62a-e0c8b9a4c682.root", "root://xrootd.hep.kbfi.ee:1094//store/data/Run2022G/EGamma/NANOAOD/22Sep2023-v2/2560000/941e4af5-473f-4ddd-a5b0-de5bdc5a5b9d.root", "root://cmsio2.rc.ufl.edu:1094//store/data/Run2022G/EGamma/NANOAOD/22Sep2023-v2/2560000/a2004054-f037-4894-ac32-26d17cfed9f4.root", "root://xrootd-cms.infn.it:1194//store/data/Run2022G/EGamma/NANOAOD/22Sep2023-v2/2560000/b16542ab-0a93-4e97-8aac-6e638e30099c.root", "root://hactar01.crc.nd.edu//store/data/Run2022G/EGamma/NANOAOD/22Sep2023-v2/2560000/b51dc8c0-ece7-4748-aee3-35b0a6abd3e8.root", "root://xrootd-cms.infn.it:1194//store/data/Run2022G/EGamma/NANOAOD/22Sep2023-v2/2560000/b5a85a95-5cce-4887-807d-40e301f3cd1b.root"]}' --output_format root