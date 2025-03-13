#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022preEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022preEE/MuonC --dataset MuonC_11 --partition_fileset '{"MuonC_11": ["root://xrootd-local.unl.edu:1094//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/30000/e4e95f0d-dce1-4351-9f49-dec9f5220cb0.root", "root://eos.cms.rcac.purdue.edu//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/30000/e57f0355-6c09-413d-8e78-c7e5ee6cedcf.root", "root://redirector.t2.ucsd.edu:1095//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/30000/e6b11689-df1b-44ca-9d0f-358a749a633e.root", "root://xrootd-cms.infn.it:1194//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/30000/e8b13317-01c6-4e27-aea2-1b0e662f330a.root", "root://xrootd-cms.infn.it:1194//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/30000/ebadb3ba-2c55-4686-a45d-a3ee26a8c953.root", "root://k8s-redir.ultralight.org:1094//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/30000/ebe6442f-f09e-4f50-8901-1816ff4e8177.root", "root://hactar01.crc.nd.edu//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/30000/ef3d2d5d-85da-4d35-ac50-f765c7faed69.root", "root://xrootd-local.unl.edu:1094//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/30000/f87d9ef4-52f5-4f9f-b466-6c93cab299f4.root", "root://eos.cms.rcac.purdue.edu//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/30000/faf64614-6cbf-4573-b9d3-57b00d762cc8.root"]}' --output_format root