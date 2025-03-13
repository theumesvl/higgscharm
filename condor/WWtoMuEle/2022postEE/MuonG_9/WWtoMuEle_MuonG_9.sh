#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022postEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022postEE/MuonG --dataset MuonG_9 --partition_fileset '{"MuonG_9": ["root://maite.iihe.ac.be:1095//store/data/Run2022G/Muon/NANOAOD/22Sep2023-v1/30000/f72e8c02-1e5d-4ce5-869a-36ebc7231052.root", "root://redirector.t2.ucsd.edu:1095//store/data/Run2022G/Muon/NANOAOD/22Sep2023-v1/30000/f7dd8213-c0bb-42a6-a455-4cf0bd64f34b.root", "root://gfe02.grid.hep.ph.ic.ac.uk:1094//pnfs/hep.ph.ic.ac.uk/data/cms/store/data/Run2022G/Muon/NANOAOD/22Sep2023-v1/30000/f841ee8c-b4a3-4399-9519-c78aabb929be.root", "root://redirector.t2.ucsd.edu:1095//store/data/Run2022G/Muon/NANOAOD/22Sep2023-v1/30000/fc541fb8-743e-4b11-8621-fe6a48ea505e.root", "root://hactar01.crc.nd.edu//store/data/Run2022G/Muon/NANOAOD/22Sep2023-v1/30000/fca0ed98-4d5f-46e7-b0ff-267b1f33f41f.root", "root://maite.iihe.ac.be:1095//store/data/Run2022G/Muon/NANOAOD/22Sep2023-v1/40000/0916c89c-b5e7-4877-a428-659ce97e226f.root", "root://xrootd-cms.infn.it:1194//store/data/Run2022G/Muon/NANOAOD/22Sep2023-v1/40000/2856ad1c-e520-4368-aac4-2fe38012b96f.root", "root://k8s-redir.ultralight.org:1094//store/data/Run2022G/Muon/NANOAOD/22Sep2023-v1/40000/e1345568-3861-4238-aae2-1957524bf53b.root", "root://gfe02.grid.hep.ph.ic.ac.uk:1094//pnfs/hep.ph.ic.ac.uk/data/cms/store/data/Run2022G/Muon/NANOAOD/22Sep2023-v1/410000/e4cc8ccb-d6fb-4145-bc86-035d9d01f87f.root"]}' --output_format root