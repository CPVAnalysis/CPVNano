#!bin/bash

cd $CMSSW_BASE/src
scram b -j 8
  
cd $CMSSW_BASE/src/PhysicsTools/CPVNano/test

#cmsRun run_nano_hnl_cfg.py isMC=False maxEvents=50
#cmsRun run_nano_phi_cfg.py isMC=False maxEvents=100

#edmFileUtil bparknano.root

#python miniAnalyser.py

#display tmp.png

  




