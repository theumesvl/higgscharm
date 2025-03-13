#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022preEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022preEE/MuonEGC --dataset MuonEGC_3 --partition_fileset '{"MuonEGC_3": ["root://xrootd-vanderbilt.sites.opensciencegrid.org:1094//store/data/Run2022C/MuonEG/NANOAOD/22Sep2023-v1/2540000/952e84ef-081c-4d54-bbdf-9ecd107f3a85.root", "root://cmsio2.rc.ufl.edu:1094//store/data/Run2022C/MuonEG/NANOAOD/22Sep2023-v1/2540000/a7be7f7d-8db9-45f7-bfe1-2ca587d059fe.root", "root://cmsio2.rc.ufl.edu:1094//store/data/Run2022C/MuonEG/NANOAOD/22Sep2023-v1/2540000/b739604e-dbba-4cdb-9292-9830df150e0c.root", "root://rdr.echo.stfc.ac.uk//store/data/Run2022C/MuonEG/NANOAOD/22Sep2023-v1/2540000/c6d101ce-dcbf-45fc-883f-8e6cde5b5105.root", "root://dcache-cms-xrootd.desy.de:1094//store/data/Run2022C/MuonEG/NANOAOD/22Sep2023-v1/2540000/d682eeb2-286d-46a1-828a-da6d2279c0ee.root", "root://hactar01.crc.nd.edu//store/data/Run2022C/MuonEG/NANOAOD/22Sep2023-v1/2540000/e2ad338b-eb1d-491b-9573-fccd2abfd7a0.root", "root://dcache-cms-xrootd.desy.de:1094//store/data/Run2022C/MuonEG/NANOAOD/22Sep2023-v1/2540000/e535669e-4df3-4741-92eb-68ad25b63d59.root", "root://cmsio2.rc.ufl.edu:1094//store/data/Run2022C/MuonEG/NANOAOD/22Sep2023-v1/2540000/ea2a50e1-ba82-4c84-91f0-1d8b8b6ec7fc.root", "root://xrootd-vanderbilt.sites.opensciencegrid.org:1094//store/data/Run2022C/MuonEG/NANOAOD/22Sep2023-v1/2540000/f534a64f-a7b5-45a5-9ae0-8c4ef2e5a2bb.root"]}' --output_format root