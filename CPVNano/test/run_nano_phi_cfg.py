# E.g.: cmsRun run_nano_hnl_cfg.py isMC=true

# info on GT: https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVAnalysisSummaryTable

from FWCore.ParameterSet.VarParsing import VarParsing
import FWCore.ParameterSet.Config as cms
from glob import glob

options = VarParsing('python')

options.register('isMC'                    , False            , VarParsing.multiplicity.singleton, VarParsing.varType.bool  , "Run this on real data"                  )
options.register('doSignal'                , True            , VarParsing.multiplicity.singleton, VarParsing.varType.bool  , "Run the BToMuMuPiBuilder"               )
#options.register('doControl'               , False           , VarParsing.multiplicity.singleton, VarParsing.varType.bool  , "Run the BToKMuMuBuilder"                )
#options.register('doHNL'                   , False           , VarParsing.multiplicity.singleton, VarParsing.varType.bool  , "Run the HNLToMuPiBuilder"               )
#options.register('doTagAndProbe'           , False           , VarParsing.multiplicity.singleton, VarParsing.varType.bool  , "Run the TagAndProbeJpsiToMuMu"          )
options.register('doGeneral'               , False           , VarParsing.multiplicity.singleton, VarParsing.varType.bool  , "Run without builder"                    )
options.register('addTriggerMuonCollection', False           , VarParsing.multiplicity.singleton, VarParsing.varType.bool  , "Add the TriggerMuon_* branches"         )
options.register('addProbeTracksCollection', False           , VarParsing.multiplicity.singleton, VarParsing.varType.bool  , "Add the ProbeTracks_* branches"         )
options.register('skipDuplicated'          , True            , VarParsing.multiplicity.singleton, VarParsing.varType.bool  , "Skip duplicated events. True by default")
options.register('globalTag'               ,'NOTSET'         , VarParsing.multiplicity.singleton, VarParsing.varType.string, "Set global tag"                         )
options.register('wantSummary'             , True            , VarParsing.multiplicity.singleton, VarParsing.varType.bool  , "Run this on real data"                  )
options.register('reportEvery'             , 1000            , VarParsing.multiplicity.singleton, VarParsing.varType.int   , "report every N events"                  )
#options.register('reportEvery'             , 1            , VarParsing.multiplicity.singleton, VarParsing.varType.int   , "report every N events"                  )
options.register('skip'                    ,  0              , VarParsing.multiplicity.singleton, VarParsing.varType.int   , "skip first N events"                    )
options.register('inputFile'               , None            , VarParsing.multiplicity.singleton, VarParsing.varType.string, "inputFile name"                         )
options.register('outFile'                 , 'bparknano.root', VarParsing.multiplicity.singleton, VarParsing.varType.string, "outputFile name"                        )
#options.register('outFile'                 , '/scratch/anlyon/tmp/nanophi_D4_nj97.root', VarParsing.multiplicity.singleton, VarParsing.varType.string, "outputFile name"                        )

options.setDefault('maxEvents', -1)

options.parseArguments()

#TODO check GT
#FIXME! update GT to ultralegacy one? Use '106X_dataRun2_v37' according to https://twiki.cern.ch/twiki/bin/view/CMS/PdmVRun2LegacyAnalysis ?
#FIXME! Update golden JSON!
globaltag = '102X_dataRun2_v11' if not options.isMC else '102X_upgrade2018_realistic_v15'
if options._beenSet['globalTag']:
    globaltag = options.globalTag

extension = {False : 'data', True : 'mc'}
outputFileNANO = cms.untracked.string('_'.join(['BParkNANO', extension[options.isMC], options.tag])+'.root')
outputFileFEVT = cms.untracked.string('_'.join(['BParkFullEvt', extension[options.isMC], options.tag])+'.root')


if not options.inputFiles:
    options.inputFiles = ['/store/data/Run2018D/ParkingBPH1/MINIAOD/05May2019promptD-v1/270000/F85EA23D-7ACA-CC47-ABA5-3F0D8DFFE32E.root'] if not options.isMC else \
                         ['file:/eos/user/a/anlyon/CPVGen/102X_crab_v0/BsToPhiPhiTo4K/crab_102X_crab_v0_BsToPhiPhiTo4K_20250130_141241/250130_131709/0000/step4_115.root']

annotation = '%s nevts:%d' % (outputFileNANO, options.maxEvents)

from Configuration.StandardSequences.Eras import eras
process = cms.Process('BParkNANO',eras.Run2_2018)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('PhysicsTools.CPVNano.nanoBPark_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# Input source
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles) if not options.inputFile else cms.untracked.vstring('file:{}'.format(options.inputFile)),
    secondaryFileNames = cms.untracked.vstring(),
    skipEvents=cms.untracked.uint32(options.skip),
    duplicateCheckMode = cms.untracked.string('checkEachFile' if options.skipDuplicated else 'noDuplicateCheck'),
)


process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(options.wantSummary),
)

process.nanoMetadata.strings.tag = annotation
# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string(annotation),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition
process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    ),
    fileName = outputFileFEVT,
    outputCommands = (cms.untracked.vstring('keep *',
                                            'drop *_*_SelectedTransient*_*',
                     )),
    splitLevel = cms.untracked.int32(0)
)

process.NANOAODoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('NANOAOD'),
        filterName = cms.untracked.string('')
    ),
    #fileName = outputFileNANO,
    fileName = outputFileNANO if not options.outFile else cms.untracked.string('file:{}'.format(options.outFile)),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep nanoaodFlatTable_*Table_*_*',     # event data
        'keep nanoaodUniqueString_nanoMetadata_*_*',   # basic metadata
        'keep nanoaodMergeableCounterTable_*Table_*_*',  # includes gentables
    )

)


# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globaltag, '')

from PhysicsTools.CPVNano.nanoBPark_cff import *
process = nanoAOD_customizeMuonTriggerBPark      (process, addTriggerMuonCollection=options.addTriggerMuonCollection)
process = nanoAOD_customizeTrackFilteredBPark    (process, addProbeTracksCollection=options.addProbeTracksCollection)
process = nanoAOD_customizeBsToPhiPhiTo4K             (process, isMC=options.isMC)
#process = nanoAOD_customizeBToMuMuPi             (process, isMC=options.isMC)
#process = nanoAOD_customizeBToKMuMu              (process, isMC=options.isMC) 
#process = nanoAOD_customizeHNLToMuPi             (process, isMC=options.isMC)
#process = nanoAOD_customizeTagAndProbeJPsiToMuMu (process, isMC=options.isMC) 
process = nanoAOD_customizeTriggerBitsBPark      (process)

# Path and EndPath definitions
process.nanoAOD_general_step = cms.Path(process.nanoSequence)
process.nanoAOD_BsToPhiPhiTo4K_step = cms.Path(process.nanoSequence + process.nanoBsToPhiPhiTo4KSequence + CountBsToPhiPhiTo4K )
#process.nanoAOD_MuMuPi_step = cms.Path(process.nanoSequence + process.nanoBMuMuPiSequence + CountBToMuMuPi )
#process.nanoAOD_KMuMu_step  = cms.Path(process.nanoSequence + process.nanoBKMuMuSequence + CountBToKmumu ) 
#process.nanoAOD_HNLToMuPi_step = cms.Path(process.nanoSequence + process.nanoHNLToMuPiSequence + CountHNLToMuPi )
#process.nanoAOD_JPsiToMuMu_step  = cms.Path(process.nanoSequence + process.nanoJPsiToMuMuSequence + CountJPsiToMuMu ) 

# customisation of the process.
if options.isMC:
    from PhysicsTools.CPVNano.nanoBPark_cff import nanoAOD_customizeMC
    nanoAOD_customizeMC(process, ancestor_particles=[531, 333], addTriggerMuonCollection=options.addTriggerMuonCollection, addProbeTracksCollection=options.addProbeTracksCollection) 

if not options.isMC:
  # apply lumi mask, see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideGoodLumiSectionsJSONFile#cmsRun
  import FWCore.PythonUtilities.LumiList as LumiList
  process.source.lumisToProcess = LumiList.LumiList(url='https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt').getVLuminosityBlockRange()

process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)
process.NANOAODoutput_step = cms.EndPath(process.NANOAODoutput)

# Schedule definition
process.schedule = cms.Schedule(
)
if options.doGeneral:
  process.schedule += cms.Schedule(
      process.nanoAOD_general_step,
  )
if options.doSignal:
  process.schedule += cms.Schedule(
      process.nanoAOD_BsToPhiPhiTo4K_step,
  )
#if options.doSignal:
#  process.schedule += cms.Schedule(
#      process.nanoAOD_MuMuPi_step,
#  )
#if options.doControl:
#  process.schedule += cms.Schedule(
#      process.nanoAOD_KMuMu_step, 
#  )
#if options.doHNL:
#  process.schedule += cms.Schedule(
#      process.nanoAOD_HNLToMuPi_step, 
#  )
#if options.doTagAndProbe:
#  process.schedule += cms.Schedule(
#      process.nanoAOD_JPsiToMuMu_step, 
#  )
process.schedule += cms.Schedule(
    process.endjob_step, 
    process.NANOAODoutput_step
)
    
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

process_string = cms.vstring()
if options.doGeneral: process_string.append('nanoAOD_general_step')
if options.doSignal: process_string.append('nanoAOD_BsToPhiPhiTo4K_step')
#if options.doSignal: process_string.append('nanoAOD_MuMuPi_step')
#if options.doControl: process_string.append('nanoAOD_KMuMu_step')
#if options.doHNL: process_string.append('nanoAOD_HNLToMuPi_step')
#if options.doTagAndProbe: process_string.append('nanoAOD_JPsiToMuMu_step')

process.NANOAODoutput.SelectEvents = cms.untracked.PSet(
    SelectEvents = process_string
)

### from https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/3287/1/1/1/1/1.html
process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))
process.NANOAODoutput.fakeNameForCrab=cms.untracked.bool(True)    

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
