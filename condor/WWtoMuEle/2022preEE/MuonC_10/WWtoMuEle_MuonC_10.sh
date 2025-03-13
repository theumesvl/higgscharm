#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022preEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022preEE/MuonC --dataset MuonC_10 --partition_fileset '{"MuonC_10": ["root://gaexrdoor.ciemat.es:1094//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/30000/c9f4728c-3d2b-4b46-93e5-05ccc6b2162e.root", "root://redirector.t2.ucsd.edu:1095//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/30000/cee8344a-a711-4091-b239-7a69a527317f.root", "root://redirector.t2.ucsd.edu:1095//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/30000/d3fdafa0-6984-4aa9-b055-123db3b334a5.root", "root://cmsio2.rc.ufl.edu:1094//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/30000/d6c07ab7-1c95-43e9-b39a-9d6ec8384740.root", "root://eos.cms.rcac.purdue.edu//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/30000/d7841f16-35f0-4343-a4a6-ecc2d4ae25ff.root", "root://maite.iihe.ac.be:1095//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/30000/d7a3603d-2278-4baf-a256-e68246e62cff.root", "root://grid-cms-xrootd.physik.rwth-aachen.de:1094//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/30000/dea86b0c-949c-4bf7-95c3-5b6a12e4aba5.root", "root://xrootd-cms.infn.it:1194//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/30000/e1db2569-fef1-4c77-9d28-8cf69aa7c46b.root", "root://dcache-cms-xrootd.desy.de:1094//store/data/Run2022C/Muon/NANOAOD/22Sep2023-v1/30000/e32a6f0c-bacb-413f-92d4-8a12ce4e571d.root"]}' --output_format root