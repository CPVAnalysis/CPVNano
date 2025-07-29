import os
import sys
sys.path.append('../data/samples')
from bparkingdata_samples_2018 import bpark_samples_2018
from bparkingdata_samples_2022 import bpark_samples_2022
from bparkingdata_samples_2024 import bpark_samples_2024


def getOptions():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='Script to launch the nanoAOD tool on top of miniAOD files', add_help=True)
  parser.add_argument('--pl'      , type=str, dest='pl'          , help='label of the sample file'                                                            , default=None)
  parser.add_argument('--ds'      , type=str, dest='ds'          , help='[optional | data/mccentral] run on specify data set. e.g "--ds D1"'                  , default=None)
  parser.add_argument('--year'    , type=str, dest='year'        , help='year to process in data'                                                             , default=None)
  parser.add_argument('--data'              , dest='data'        , help='run the nano tool on a data sample'                             , action='store_true', default=False)
  return parser.parse_args()
  

def checkParser(opt):
  if opt.pl==None:
    raise RuntimeError('Please indicate the production label: for --mcprivate, it has to correspond to the label of the miniAOD')


class NanoLauncher(object):
  def __init__(self, opt):
    self.prodlabel = vars(opt)['pl']
    self.year      = vars(opt)['year']
    self.ds        = vars(opt)['ds']
    self.data      = vars(opt)['data']

    self.do_all = False
    if self.data and self.ds == None:
      self.do_all = True

    if self.data:
      if self.year not in ['2018', '2022', '2024']:
        raise RuntimeError("Wrong year '{}'. Please choose amongst ['2018', '2022', '2024']'".format(self.year))

      if self.year == '2018':
        self.bpark_samples = bpark_samples_2018
      elif self.year == '2022':
        self.bpark_samples = bpark_samples_2022
      elif self.year == '2024':
        self.bpark_samples = bpark_samples_2024

      if self.ds == None:
        self.keys = self.bpark_samples.keys()

      else: # process a specific data set
        if self.ds not in self.bpark_samples.keys():
          raise RuntimeError('Please indicate on which period of the BParking dataset you want to run. Label "{}" not recognised. Choose among {}'.format(self.ds, self.bpark_samples.keys()))

        self.keys = [self.ds]


  def prepare_nano_config(self):
    doSignal = 'True' 
    doGeneral = 'False'
    addTriggerMuonCollection = 'True'
    addProbeTracksCollection = 'False'

    if self.data:
      isMC = 'False'
      if self.year == '2018':
        json_file = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'
        gt = '102X_dataRun2_v11'
        era = 'Run2_2018'
      else:
        raise RuntimeError('Please insert GT, json and era for year {}'.format(self.year))

    config = [
      "import FWCore.ParameterSet.Config as cms",
      "",
      "isMC = {}".format(isMC),
      "globaltag = '{}'".format(gt),
      "json_file = '{}'".format(json_file),
      "doSignal = {}".format(doSignal),
      "doGeneral = {}".format(doGeneral),
      "addTriggerMuonCollection = cms.untracked.bool({})".format(addTriggerMuonCollection),
      "addProbeTracksCollection = cms.untracked.bool({})".format(addProbeTracksCollection),
      "reportEvery = 1000",
      "",
      "from Configuration.StandardSequences.Eras import eras",
      "process = cms.Process('BParkNANO',eras.{})".format(era),
      "",
      "process.load('Configuration.StandardSequences.Services_cff')",
      "process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')",
      "process.load('FWCore.MessageService.MessageLogger_cfi')",
      "process.load('Configuration.EventContent.EventContent_cff')",
      "process.load('Configuration.StandardSequences.GeometryRecoDB_cff')",
      "process.load('Configuration.StandardSequences.MagneticField_cff')",
      "process.load('PhysicsTools.CPVNano.nanoBPark_cff')",
      "process.load('Configuration.StandardSequences.EndOfProcess_cff')",
      "process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')",
      "",
      "process.MessageLogger.cerr.FwkReport.reportEvery = reportEvery",
      "process.maxEvents = cms.untracked.PSet(",
      "    input = cms.untracked.int32(-1)",
      ")",
      "",
      "process.source = cms.Source(",
      "    'PoolSource',",
      "    fileNames = cms.untracked.vstring(),",
      "    secondaryFileNames = cms.untracked.vstring(),",
      "    skipEvents=cms.untracked.uint32(0),",
      "    duplicateCheckMode = cms.untracked.string('checkEachFile'),",
      ")",
      "",
      "process.options = cms.untracked.PSet(",
      "    wantSummary = cms.untracked.bool(True),",
      ")",
      "",
      "outputFileNANO = cms.untracked.string('BParkNANO.root')",
      "annotation = '%s nevts:%d' % (outputFileNANO, -1)",
      "process.nanoMetadata.strings.tag = annotation",
      "process.configurationMetadata = cms.untracked.PSet(",
      "    annotation = cms.untracked.string(annotation),",
      "    name = cms.untracked.string('Applications'),",
      "    version = cms.untracked.string('$Revision: 1.19 $')",
      ")",
      "",
      "outputFileFEVT = cms.untracked.string('BParkFullEvt.root')",
      "process.FEVTDEBUGHLToutput = cms.OutputModule('PoolOutputModule',",
      "    dataset = cms.untracked.PSet(",
      "        dataTier = cms.untracked.string('GEN-SIM-RECO'),",
      "        filterName = cms.untracked.string('')",
      "    ),",
      "    fileName = outputFileFEVT,",
      "    outputCommands = (cms.untracked.vstring('keep *',",
      "                                            'drop *_*_SelectedTransient*_*',",
      "                     )),",
      "    splitLevel = cms.untracked.int32(0)",
      ")",
      "",
      "process.NANOAODoutput = cms.OutputModule('NanoAODOutputModule',",
      "    compressionAlgorithm = cms.untracked.string('LZMA'),",
      "    compressionLevel = cms.untracked.int32(9),",
      "    dataset = cms.untracked.PSet(",
      "        dataTier = cms.untracked.string('NANOAOD'),",
      "        filterName = cms.untracked.string('')",
      "    ),",
      "    fileName = cms.untracked.string('file:bparknano.root'),",
      "    outputCommands = cms.untracked.vstring(",
      "        'drop *',",
      "        'keep nanoaodFlatTable_*Table_*_*',",
      "        'keep nanoaodUniqueString_nanoMetadata_*_*',",
      "        'keep nanoaodMergeableCounterTable_*Table_*_*',",
      "    )",
      "",
      ")",
      "",
      "from Configuration.AlCa.GlobalTag import GlobalTag",
      "process.GlobalTag = GlobalTag(process.GlobalTag, globaltag, '')",
      "",
      "from PhysicsTools.CPVNano.nanoBPark_cff import *",
      "process = nanoAOD_customizeMuonTriggerBPark      (process, addTriggerMuonCollection=addTriggerMuonCollection)",
      "process = nanoAOD_customizeTrackFilteredBPark    (process, addProbeTracksCollection=addProbeTracksCollection)",
      "process = nanoAOD_customizeBsToPhiPhiTo4K        (process, isMC=isMC)",
      "process = nanoAOD_customizeTriggerBitsBPark      (process)",
      "",
      "process.nanoAOD_general_step = cms.Path(process.nanoSequence)",
      "process.nanoAOD_BsToPhiPhiTo4K_step = cms.Path(process.nanoSequence + process.nanoBsToPhiPhiTo4KSequence + CountBsToPhiPhiTo4K)",
      "",
      "if isMC:",
      "    from PhysicsTools.CPVNano.nanoBPark_cff import nanoAOD_customizeMC",
      "    nanoAOD_customizeMC(process, ancestor_particles=[531, 333], addTriggerMuonCollection=addTriggerMuonCollection, addProbeTracksCollection=addProbeTracksCollection)",
      "",
      "if not isMC:",
      "  import FWCore.PythonUtilities.LumiList as LumiList",
      "  process.source.lumisToProcess = LumiList.LumiList(url=json_file).getVLuminosityBlockRange()",
      "",
      "process.endjob_step = cms.EndPath(process.endOfProcess)",
      "process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)",
      "process.NANOAODoutput_step = cms.EndPath(process.NANOAODoutput)",
      "",
      "process.schedule = cms.Schedule(",
      ")",
      "if doGeneral:",
      "  process.schedule += cms.Schedule(",
      "      process.nanoAOD_general_step,",
      "  )",
      "if doSignal:",
      "  process.schedule += cms.Schedule(",
      "      process.nanoAOD_BsToPhiPhiTo4K_step,",
      "  )",
      "",
      "process.schedule += cms.Schedule(",
      "    process.endjob_step,", 
      "    process.NANOAODoutput_step",
      ")",
      "",    
      "from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask",
      "associatePatAlgosToolsTask(process)",
      "",
      "process_string = cms.vstring()",
      "if doGeneral: process_string.append('nanoAOD_general_step')",
      "if doSignal: process_string.append('nanoAOD_BsToPhiPhiTo4K_step')",
      "",
      "process.NANOAODoutput.SelectEvents = cms.untracked.PSet(",
      "    SelectEvents = process_string",
      ")",
      "",
      "process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))",
      "process.NANOAODoutput.fakeNameForCrab=cms.untracked.bool(True)",    
      "",
      "process.load('TrackingTools/TransientTrack/TransientTrackBuilder_cfi')",
      "from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete",
      "process = customiseEarlyDelete(process)",
      ]

    config = '\n'.join(config)

    f_out = open('the_nano_config.py', 'w+')
    f_out.write(config)
    f_out.close()

    print('\t--> the_nano_config.py created')


  def prepare_CRAB_config(self, key):
    dataset = self.bpark_samples[key] 

    label = '{}_{}'.format(self.prodlabel, key)

    if self.data:
      outdir = '/store/group/phys_bphys/anlyon/CPVGen/data/{a}/{b}/{c}'.format(
          a = self.prodlabel,
          b = self.year, 
          c = key,
          )

    units_per_job = 5 # 10

    config = [
      "from CRABClient.UserUtilities import config, ClientException",
      "import yaml",
      "import datetime",
      "from fnmatch import fnmatch",
      "from argparse import ArgumentParser",
      "import time",

      "label = '{}'".format(label),

      "ts = time.time()",
      "date_time = datetime.datetime.fromtimestamp(ts).strftime('%Y%m%d_%H%M%S')",

      "config = config()",
      "config.section_('General')",
      "config.General.requestName = '{}_{}'.format(label, date_time)",
      "config.General.transferOutputs = True",        
      "config.General.transferLogs    = True",        
      "config.General.workArea        = 'crab_workdir/{}/{}'.format(label, date_time)",   

      "config.section_('Data')",
      "config.Data.inputDataset       = '{}'".format(dataset),
      "config.Data.publication        = False",
      "config.Data.outLFNDirBase      = '{}'".format(outdir),
      "config.Data.inputDBS           = 'global'",
      "config.Data.allowNonValidInputDataset = True",
      "config.Data.splitting          = 'FileBased'",
      "config.Data.unitsPerJob        = {}".format(units_per_job),         
      #"config.Data.totalUnits = 1 ", #FIXME

      "config.section_('JobType')",
      "config.JobType.pluginName      = 'Analysis'",       
      "config.JobType.psetName        =  'the_nano_config.py'", 
      "config.JobType.inputFiles      = ['the_nano_config.py']", 
      "config.JobType.allowUndistributedCMSSW = True", 
      "config.JobType.outputFiles     = ['bparknano.root']", 
      "config.JobType.maxMemoryMB = 2500",
      "#config.JobType.maxJobRuntimeMin = 3000 #TODO add",

      "config.section_('Site')",
      "config.Site.storageSite = 'T2_CH_CERN'",
      ]

    config = '\n'.join(config)

    f_out = open('the_crab_config.py', 'w+')
    f_out.write(config)
    f_out.close()

    print('\t--> the_crab_config.py created')


  def submit(self):
      command = 'crab submit -c the_crab_config.py'
      os.system(command)

      command_rm_nano = 'rm the_nano_config.py'
      os.system(command_rm_nano)

      command_rm_crab = 'rm the_crab_config.py'
      os.system(command_rm_crab)


  def process(self):
    print('\n------------')
    print(' Processing NanoLauncher on production {} '.format(self.prodlabel))
    print('------------')

    if self.data:
      for key in self.keys:
        dataset = self.bpark_samples[key]
        print('\n-> Processing data set {}'.format(dataset)) 

        print('\n   -> Preparing nano config') 
        self.prepare_nano_config() 

        print('\n   -> Preparing CRAB config') 
        self.prepare_CRAB_config(key=key) 

        print('\n   -> Submitting...') 
        self.submit()

        print('\n   -> Submission completed')

    print('\nDone')



if __name__ == "__main__":
  opt = getOptions()
  checkParser(opt)
  NanoLauncher(opt).process()


  
