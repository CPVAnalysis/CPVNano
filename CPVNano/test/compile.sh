#!bin/bash

cd $CMSSW_BASE/src
scram b -j 8
  
cd $CMSSW_BASE/src/PhysicsTools/CPVNano/test

