#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=/afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/t/tvanlaer/private/x509up_u158952
cd /afs/cern.ch/user/t/tvanlaer/Hc/higgscharm

python3 submit.py --processor WWtoMuEle --year 2022preEE --output_path /eos/user/t/tvanlaer/higgscharm/outputs/WWtoMuEle/2022preEE/WGtoLNuG-PTG-400to600 --dataset WGtoLNuG-PTG-400to600_2 --partition_fileset '{"WGtoLNuG-PTG-400to600_2": ["root://cmsdcache-kit-disk.gridka.de:1094//store/mc/Run3Summer22NanoAODv12/WGtoLNuG-1Jets_PTG-400to600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_v5-v1/70000/cf79542f-2826-4011-b76c-9c32657fdb70.root", "root://gaexrdoor.ciemat.es:1094//store/mc/Run3Summer22NanoAODv12/WGtoLNuG-1Jets_PTG-400to600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_v5-v1/70000/d0d121a9-4406-47fc-9231-e05d5ab931d7.root", "root://grid-cms-xrootd.physik.rwth-aachen.de:1094//store/mc/Run3Summer22NanoAODv12/WGtoLNuG-1Jets_PTG-400to600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_v5-v1/70000/d6573b13-7f76-4945-a462-ddd22ccc244f.root", "root://grid-cms-xrootd.physik.rwth-aachen.de:1094//store/mc/Run3Summer22NanoAODv12/WGtoLNuG-1Jets_PTG-400to600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_v5-v1/70000/d7847c40-b57f-4232-a1bf-bf7afb76a89d.root", "root://grid-cms-xrootd.physik.rwth-aachen.de:1094//store/mc/Run3Summer22NanoAODv12/WGtoLNuG-1Jets_PTG-400to600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_v5-v1/70000/e8f31c08-58c9-4509-8256-25c452aac6f4.root"]}' --output_format root