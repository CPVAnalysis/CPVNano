# NanoAOD tool for CPV analysis

## Installation

Setup the environment

```
cmssw-el7
cmsrel CMSSW_10_2_15
cd CMSSW_10_2_15/src
cmsenv
git cms-init
```

Add the NanoAOD tool

```
git clone git@github.com:CPVAnalysis/CPVNano.git ./PhysicsTools
```

Import custom CMSSW modifications

```
git cms-merge-topic -u amlyon:BHNLNano
```

Build everything
```
scram b -j 8
```

