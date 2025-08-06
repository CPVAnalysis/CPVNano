from __future__ import print_function
import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.globals_cff import *
from PhysicsTools.NanoAOD.nano_cff import *
from PhysicsTools.NanoAOD.vertices_cff import *
from PhysicsTools.NanoAOD.NanoAODEDMEventContent_cff import *



##for gen and trigger muon
from PhysicsTools.CPVNano.genparticlesBPark_cff import finalGenParticlesBPark, genParticleBParkTable, genParticleBParkSequence, genParticleBParkTables
from PhysicsTools.CPVNano.particlelevelBPark_cff import *
from PhysicsTools.CPVNano.primaryverticesBPark_cff import *
from PhysicsTools.CPVNano.muonsBPark_cff import * 

## filtered input collections
from PhysicsTools.CPVNano.tracksBPark_cff import *

## B collections
from PhysicsTools.CPVNano.BsToPhiPhiTo4K_cff import *


vertexTable.svSrc = cms.InputTag("slimmedSecondaryVertices")

nanoSequence = cms.Sequence(nanoMetadata + 
                            cms.Sequence(vertexTask) +
                            cms.Sequence(globalTablesTask)+ 
                            cms.Sequence(vertexTablesTask) +
                            cms.Sequence(pVertexTable) 
                           )

nanoSequenceMC = cms.Sequence(particleLevelBParkSequence + 
                              genParticleBParkSequence + 
                              #globalTablesMC + 
                              genWeightsTable + 
                              genParticleBParkTables + 
                              lheInfoTable
                             ) 

# from BPHNano
#def nanoAOD_customizeMC(process):
#    process.nanoSequence = cms.Sequence(process.nanoSequence +particleLevelBPHSequence + genParticleBPHSequence+ genParticleBPHTables )
#    return process

def nanoAOD_customizeMuonTriggerBPark(process, addTriggerMuonCollection=False):
    if addTriggerMuonCollection:
      process.nanoSequence = cms.Sequence( process.nanoSequence + muonBParkSequence + muonTriggerMatchedTables + muonBParkTables)
    else:
      process.nanoSequence = cms.Sequence( process.nanoSequence + muonBParkSequence + muonBParkTables)
    return process

def nanoAOD_customizeTrackFilteredBPark(process, addProbeTracksCollection=False):
    if addProbeTracksCollection:
      process.nanoSequence = cms.Sequence( process.nanoSequence + tracksBParkSequence + tracksBParkTables)
    else:
      process.nanoSequence = cms.Sequence( process.nanoSequence + tracksBParkSequence)
    return process

def nanoAOD_customizeBsToPhiPhiTo4K(process, isMC=False):
    if isMC == False:
      process.nanoBsToPhiPhiTo4KSequence = cms.Sequence( PhiToKKSequence + BsToPhiPhiTo4KSequence + PhiToKKTable + BsToPhiPhiTo4KTable + CountPhiToKK)
    else:
      process.nanoBsToPhiPhiTo4KSequence = cms.Sequence( PhiToKKSequenceMC + BsToPhiPhiTo4KSequenceMC + PhiToKKTable + BsToPhiPhiTo4KTable )
    return process

from FWCore.ParameterSet.MassReplace import massSearchReplaceAnyInputTag
def nanoAOD_customizeMC(process, ancestor_particles=[531, 333], addTriggerMuonCollection=False, addProbeTracksCollection=False):  
    for name, path in process.paths.iteritems():
        # replace all the non-match embedded inputs with the matched ones
        massSearchReplaceAnyInputTag(path, 'muonTrgSelector:SelectedMuons', 'selectedMuonsMCMatchEmbedded')
        massSearchReplaceAnyInputTag(path, 'muonTrgSelector:trgMuons', 'triggerMuonsMCMatchEmbedded')
        massSearchReplaceAnyInputTag(path, 'tracksBPark:SelectedTracks', 'tracksBParkMCMatchEmbedded')

        # make the PhiToKKTable/count talk to the correct producer
        massSearchReplaceAnyInputTag(path, 'PhiToKK', 'PhiToKKMC')

        # make the BsToPhiPhiTo4KTable/count talk to the correct producer
        massSearchReplaceAnyInputTag(path, 'BsToPhiPhiTo4K', 'BsToPhiPhiTo4KMC')

        # save the all descendants of ancestor_particles
        to_save = ' || '.join(['abs(pdgId) == %d'%ipdg for ipdg in ancestor_particles])
        to_save = 'keep++ (%s)'%to_save
        finalGenParticlesBPark.select.append(to_save)

        # modify the path to include mc-specific info
        path.insert(0, nanoSequenceMC)
        if addTriggerMuonCollection:
          path.replace(process.muonBParkSequence, process.muonBParkMCWithTriggerMuon)
        else:
          path.replace(process.muonBParkSequence, process.muonBParkMC)
        if addProbeTracksCollection:
          path.replace(process.tracksBParkSequence, process.tracksBParkMCWithTable)
        else:
          path.replace(process.tracksBParkSequence, process.tracksBParkMC)
