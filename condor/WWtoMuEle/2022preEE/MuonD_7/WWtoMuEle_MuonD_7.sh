#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022preEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022preEE/MuonD --dataset MuonD_7 --partition_fileset '{"MuonD_7": ["root://hactar01.crc.nd.edu//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/a60b91ba-9901-427c-8b55-3b14e0ce0da6.root", "root://hactar01.crc.nd.edu//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/a7f0cc2c-3014-49da-a91a-ce299d2f8eb2.root", "root://cmsio2.rc.ufl.edu:1094//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/a850a486-caac-4fa1-a4cb-146c35736965.root", "root://xrootd-vanderbilt.sites.opensciencegrid.org:1094//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/a98f8f69-5552-4d37-a3a8-ebd9caef3dd0.root", "root://grid-cms-xrootd.physik.rwth-aachen.de:1094//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/aee48c5d-76de-4173-979c-120b9c2f62ed.root", "root://hactar01.crc.nd.edu//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/af844fe7-965b-4cf3-8f2a-e788d3df2d0f.root", "root://cmsio2.rc.ufl.edu:1094//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/c6a0d3dd-6607-481d-94cb-71511a457446.root", "root://cceos.ihep.ac.cn:1094//eos/ihep/cms/store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/d167d633-2245-4e65-b3f7-04a5b8c4d00e.root", "root://cceos.ihep.ac.cn:1094//eos/ihep/cms/store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/d7167490-5e81-44c8-b6db-86112197aa8c.root"]}' --output_format root