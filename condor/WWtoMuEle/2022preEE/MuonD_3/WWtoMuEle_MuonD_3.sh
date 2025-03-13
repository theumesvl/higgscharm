#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022preEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022preEE/MuonD --dataset MuonD_3 --partition_fileset '{"MuonD_3": ["root://maite.iihe.ac.be:1095//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/3458aef6-ca7f-4bd5-a254-e91ece23155f.root", "root://cceos.ihep.ac.cn:1094//eos/ihep/cms/store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/3e0140a0-ae18-4633-a34b-fa436cb19bba.root", "root://redirector.t2.ucsd.edu:1095//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/40645b8d-dded-4ac0-a5f9-fb9e0b6ad7cb.root", "root://hactar01.crc.nd.edu//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/438e3957-889e-42d2-80d3-19b9eed364df.root", "root://xrootd-vanderbilt.sites.opensciencegrid.org:1094//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/45ccc208-8883-49d3-81d9-78a13fb2244d.root", "root://dcache-cms-xrootd.desy.de:1094//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/48351160-e06d-4585-9980-9616220d9884.root", "root://cmsio2.rc.ufl.edu:1094//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/49306e91-5df0-4f05-855f-5f9b40d0277d.root", "root://gaexrdoor.ciemat.es:1094//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/4d377162-8af2-41cf-b963-e8b15d6acdbf.root", "root://cmsdcache-kit-disk.gridka.de:1094//store/data/Run2022D/Muon/NANOAOD/22Sep2023-v1/2520000/50eb4dec-d0c0-400a-a1eb-4561043b266c.root"]}' --output_format root