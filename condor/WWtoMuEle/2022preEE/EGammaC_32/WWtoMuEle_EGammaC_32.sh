#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022preEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022preEE/EGammaC --dataset EGammaC_32 --partition_fileset '{"EGammaC_32": ["root://grid-cms-xrootd.physik.rwth-aachen.de:1094//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/50000/db5d3133-3687-44b7-8117-6b661b5e02c2.root", "root://redirector.t2.ucsd.edu:1095//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/50000/df908dec-7397-4687-b9ff-8cef1f04832c.root", "root://dcache-cms-xrootd.desy.de:1094//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/50000/df9d30f9-3fb5-4fb0-be0b-bcd091db0581.root", "root://redirector.t2.ucsd.edu:1095//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/50000/e2775a27-a1d5-4909-ba4d-e86526cae1f6.root", "root://xrootd-vanderbilt.sites.opensciencegrid.org:1094//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/50000/e9dcd001-1a52-4a1b-b251-1cead74883a5.root", "root://cmsdcadisk.fnal.gov//dcache/uscmsdisk/store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/50000/f48fecb3-769e-49c0-9379-4c07fa52223d.root", "root://cmsdcache-kit-disk.gridka.de:1094//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/50000/fa5c75dc-fde0-4198-87f0-88d7d189c29b.root", "root://grid-cms-xrootd.physik.rwth-aachen.de:1094//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/510000/61ef47ab-8a9e-464b-8702-1cffa353639b.root", "root://eoscms.cern.ch//eos/cms/store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/510000/dc006db3-45ea-4199-b760-b0e498fd8b52.root"]}' --output_format root