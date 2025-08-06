import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.genparticles_cff import *
from PhysicsTools.NanoAOD.simpleGenParticleFlatTableProducer_cfi import simpleGenParticleFlatTableProducer


# for BParkPark start with merged particles (pruned + packed),
# where pruned contain K* states, but not final states, 
# and packed contain final states (K pi).
# then you save also final states (granddaughters)
finalGenParticlesBPark = finalGenParticles.clone(
  src = cms.InputTag("mergedGenParticles"),
  select = cms.vstring(
        "drop *",
        "keep++ (abs(pdgId) == 511 || abs(pdgId) == 521 || abs(pdgId)==531)",  #keep all B0(=511) and B+/-(521) + their daughters and granddaughters
   )
)

genParticleBParkTable = simpleGenParticleFlatTableProducer.clone(
  src = cms.InputTag("finalGenParticlesBPark"),
  name = cms.string("GenPart"),
  doc = cms.string("interesting gen particles for BPark"),  
  variables = cms.PSet(
      genParticleTable.variables,
      vx = Var("vx", float, doc="x coordinate of the production vertex position, in cm", precision=12),
      vy = Var("vy", float, doc="y coordinate of the production vertex position, in cm", precision=12),
      vz = Var("vz", float, doc="z coordinate of the production vertex position, in cm", precision=12),
  )
)


genParticleBParkSequence = cms.Sequence(finalGenParticlesBPark)
genParticleBParkTables = cms.Sequence(genParticleBParkTable)

