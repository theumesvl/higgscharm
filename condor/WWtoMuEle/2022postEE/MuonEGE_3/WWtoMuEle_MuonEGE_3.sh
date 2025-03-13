#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022postEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022postEE/MuonEGE --dataset MuonEGE_3 --partition_fileset '{"MuonEGE_3": ["root://cmsio2.rc.ufl.edu:1094//store/data/Run2022E/MuonEG/NANOAOD/22Sep2023-v1/2520000/bb2e2dd4-c350-422a-b6cc-34b6d4d79755.root", "root://cmsio2.rc.ufl.edu:1094//store/data/Run2022E/MuonEG/NANOAOD/22Sep2023-v1/2520000/d0da2e07-d1cc-461a-9ae1-8ac39e25a7e3.root", "root://maite.iihe.ac.be:1095//store/data/Run2022E/MuonEG/NANOAOD/22Sep2023-v1/2520000/e415a881-5467-435c-9c43-cfcc42e04dda.root", "root://t2dsk0011.cmsaf.mit.edu:1094//store/data/Run2022E/MuonEG/NANOAOD/22Sep2023-v1/2520000/e9e01a6f-9fc5-4b3a-8304-d4f4f90901e1.root", "root://hactar01.crc.nd.edu//store/data/Run2022E/MuonEG/NANOAOD/22Sep2023-v1/2520000/ef8f4d53-804d-469a-828c-e855d047a652.root", "root://cmsio2.rc.ufl.edu:1094//store/data/Run2022E/MuonEG/NANOAOD/22Sep2023-v1/2520000/f01c071c-3e7a-43db-ab0d-2db72dc3c10f.root", "root://maite.iihe.ac.be:1095//store/data/Run2022E/MuonEG/NANOAOD/22Sep2023-v1/2520000/f82989b5-553d-4684-9899-fb5f7604d3ee.root", "root://t2dsk0011.cmsaf.mit.edu:1094//store/data/Run2022E/MuonEG/NANOAOD/22Sep2023-v1/2520000/fe22e256-9805-4bed-9945-227660abd39b.root", "root://cmsio2.rc.ufl.edu:1094//store/data/Run2022E/MuonEG/NANOAOD/22Sep2023-v1/2520000/fe4cc0ea-dad6-4ef1-8d5f-cfd64eda892e.root"]}' --output_format root