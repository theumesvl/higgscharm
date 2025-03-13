#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022preEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022preEE/MuonD --dataset MuonD_2 --partition_fileset '{"MuonD_2": ["root://dcache-cms-xrootd.desy.de:1094//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/203cd43e-8323-438c-8274-dd7bdb4b13e5.root", "root://xrootd-local.unl.edu:1094//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/26c4bce2-79ed-4af0-81fc-6b79721c24a3.root", "root://eos.cms.rcac.purdue.edu//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/284a0433-31e7-46a9-8db7-abf9be9b1d35.root", "root://dcache-cms-xrootd.desy.de:1094//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/29f4b95b-802e-4daa-b95d-7cbee41c5315.root", "root://hactar01.crc.nd.edu//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/2b14ab6a-6b70-4ec1-bf48-d04f87f5173d.root", "root://xrootd-vanderbilt.sites.opensciencegrid.org:1094//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/2c21c20f-e7df-4f80-9c0c-f2c7dc170ca6.root", "root://cmsio2.rc.ufl.edu:1094//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/2fb19461-7c85-43ff-94b3-733369e4e6fa.root", "root://redirector.t2.ucsd.edu:1095//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/2fe0b89f-4fef-4cd0-b81a-eb5c8964194e.root", "root://cmsdcache-kit-disk.gridka.de:1094//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/3209afad-e53d-4540-869a-fe03f71914a5.root"]}' --output_format root