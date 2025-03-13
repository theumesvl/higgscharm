#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022preEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022preEE/MuonEGC --dataset MuonEGC_2 --partition_fileset '{"MuonEGC_2": ["root://gaexrdoor.ciemat.es:1094//store/data/Run2022C/MuonEG/NANOAOD/22Sep2023-v1/2540000/36a2f450-209b-4638-af77-5d2fe5ba639c.root", "root://dcache-cms-xrootd.desy.de:1094//store/data/Run2022C/MuonEG/NANOAOD/22Sep2023-v1/2540000/40f4d272-1551-49a2-8e72-0aa8b0d2ae28.root", "root://hactar01.crc.nd.edu//store/data/Run2022C/MuonEG/NANOAOD/22Sep2023-v1/2540000/51fff01e-8657-437c-8806-43b70d5f14b8.root", "root://xrootd-vanderbilt.sites.opensciencegrid.org:1094//store/data/Run2022C/MuonEG/NANOAOD/22Sep2023-v1/2540000/56be846c-875f-4d48-b214-54ba81f5d5d2.root", "root://dcache-cms-xrootd.desy.de:1094//store/data/Run2022C/MuonEG/NANOAOD/22Sep2023-v1/2540000/5e968b63-14c6-48e7-b19f-238a2bd1d1a4.root", "root://xrootd-vanderbilt.sites.opensciencegrid.org:1094//store/data/Run2022C/MuonEG/NANOAOD/22Sep2023-v1/2540000/61c79be4-45e0-4c5c-926b-3ad15a71b2b7.root", "root://gaexrdoor.ciemat.es:1094//store/data/Run2022C/MuonEG/NANOAOD/22Sep2023-v1/2540000/8913c9e9-9552-4d00-b40c-c348df1888a0.root", "root://dcache-cms-xrootd.desy.de:1094//store/data/Run2022C/MuonEG/NANOAOD/22Sep2023-v1/2540000/91ebcb9d-19a8-4328-9095-7d6d13f87f76.root", "root://rdr.echo.stfc.ac.uk//store/data/Run2022C/MuonEG/NANOAOD/22Sep2023-v1/2540000/92c52dd8-4272-413c-a6cc-ae6ea8b1218e.root"]}' --output_format root