#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022preEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022preEE/MuonEGD --dataset MuonEGD_2 --partition_fileset '{"MuonEGD_2": ["root://osg-se.sprace.org.br:1094//store/data/Run2022D/MuonEG/NANOAOD/22Sep2023-v1/40000/73fc2757-5d10-4835-b7a0-03eaf208aa4a.root", "root://grid143.kfki.hu:1094//store/data/Run2022D/MuonEG/NANOAOD/22Sep2023-v1/40000/7a7b7510-00a8-46e1-b456-8ac4f443cd05.root", "root://hactar01.crc.nd.edu//store/data/Run2022D/MuonEG/NANOAOD/22Sep2023-v1/40000/a74255fa-4b70-4ef0-9d47-1ee2651ac525.root", "root://hactar01.crc.nd.edu//store/data/Run2022D/MuonEG/NANOAOD/22Sep2023-v1/40000/b01af2f6-9d7a-4ea1-815a-6386d6842bca.root", "root://osg-se.sprace.org.br:1094//store/data/Run2022D/MuonEG/NANOAOD/22Sep2023-v1/40000/ba818ba7-3867-4041-9fc8-645a088637e4.root", "root://eos.cms.rcac.purdue.edu//store/data/Run2022D/MuonEG/NANOAOD/22Sep2023-v1/40000/dcd5fd93-5c65-4483-8f22-c9a2e284ff75.root", "root://cmsdcache-kit-disk.gridka.de:1094//store/data/Run2022D/MuonEG/NANOAOD/22Sep2023-v1/40000/f370bcdb-1bb6-4df3-9c6f-c86f733f79b4.root", "root://eos.cms.rcac.purdue.edu//store/data/Run2022D/MuonEG/NANOAOD/22Sep2023-v1/40000/f84ea3df-28a3-46d4-a999-d76309826592.root"]}' --output_format root