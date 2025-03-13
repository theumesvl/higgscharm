#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022postEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022postEE/EGammaG --dataset EGammaG_9 --partition_fileset '{"EGammaG_9": ["root://grid-cms-xrootd.physik.rwth-aachen.de:1094//store/data/Run2022G/EGamma/NANOAOD/22Sep2023-v2/2560000/ed644111-b3b0-4b17-93ff-85409d01c8ba.root", "root://rdr.echo.stfc.ac.uk//store/data/Run2022G/EGamma/NANOAOD/22Sep2023-v2/2560000/ee3f48f8-f0b1-4127-8e99-8e6cfc45a1a4.root", "root://gaexrdoor.ciemat.es:1094//store/data/Run2022G/EGamma/NANOAOD/22Sep2023-v2/2560000/f309ccdc-fbab-4533-91ad-156911a675d2.root", "root://grid-cms-xrootd.physik.rwth-aachen.de:1094//store/data/Run2022G/EGamma/NANOAOD/22Sep2023-v2/2560000/f8810663-01aa-4dbc-8b2c-e870d32267f9.root", "root://ruhex-osgce.rutgers.edu//store/data/Run2022G/EGamma/NANOAOD/22Sep2023-v2/2560000/fc70d472-858a-481e-b717-72d755f76d63.root", "root://xrootd.hep.kbfi.ee:1094//store/data/Run2022G/EGamma/NANOAOD/22Sep2023-v2/70000/1420884d-0a93-4d02-8167-e8576cef4aaa.root", "root://redirector.t2.ucsd.edu:1095//store/data/Run2022G/EGamma/NANOAOD/22Sep2023-v2/70000/73239860-5c98-4b4e-8d94-206587ff661f.root", "root://storage01.lcg.cscs.ch:1096//pnfs/lcg.cscs.ch/cms/trivcat/store/data/Run2022G/EGamma/NANOAOD/22Sep2023-v2/70000/90e703ae-a5b2-49c3-9470-73016480967e.root", "root://t3se01.psi.ch:1094//store/data/Run2022G/EGamma/NANOAOD/22Sep2023-v2/80000/4630fcad-899f-4ece-bc34-46affd20f268.root"]}' --output_format root