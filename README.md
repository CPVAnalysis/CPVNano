# NanoAOD tool for CPV analysis

## Installation

Setup the environment

```
cmsrel CMSSW_15_0_10
cd CMSSW_15_0_10/src
cmsenv
git cms-init
```

Add the CPVNano framework

```
git clone git@github.com:CPVAnalysis/CPVNano.git ./PhysicsTools
```

Add the NanoAOD CMSSW package

```
git cms-addpkg PhysicsTools/NanoAOD
```

Build everything

```
scram b -j 8
```

## After first installation

```
cd CMSSW_15_0_10/src/PhysicsTools/CPVNano
source setup.sh
```

Activate the proxy

```
voms-proxy-init --voms cms --valid 186:00
```

