#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022preEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022preEE/EGammaC --dataset EGammaC_31 --partition_fileset '{"EGammaC_31": ["root://redirector.t2.ucsd.edu:1095//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/50000/8ec80483-9458-479a-9704-827fff82b6b2.root", "root://cmsdcache-kit-disk.gridka.de:1094//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/50000/9036fd69-6ea4-4936-a62c-9d207f0d442c.root", "root://maite.iihe.ac.be:1095//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/50000/9913c83c-32a9-4576-85d5-9744ea6590bf.root", "root://grid-cms-xrootd.physik.rwth-aachen.de:1094//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/50000/a1918aba-bfdb-44da-b6b7-9a27aa89d121.root", "root://storage01.lcg.cscs.ch:1096//pnfs/lcg.cscs.ch/cms/trivcat/store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/50000/bcde1751-7ed4-44b4-bd13-612bcd3829c0.root", "root://storage01.lcg.cscs.ch:1096//pnfs/lcg.cscs.ch/cms/trivcat/store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/50000/bfefa69d-4855-4417-83b7-42fcf89a7b2e.root", "root://eos.cms.rcac.purdue.edu//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/50000/ce0e7453-9827-40f7-9b4c-dd950a850b88.root", "root://t3se01.psi.ch:1094//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/50000/ce823b4f-ab9b-43ed-9668-5e6a7c3bccfe.root", "root://hactar01.crc.nd.edu//store/data/Run2022C/EGamma/NANOAOD/22Sep2023-v1/50000/d3b12439-e261-4233-b9b5-b329d1183e15.root"]}' --output_format root