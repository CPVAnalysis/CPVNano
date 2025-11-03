import os
import sys
sys.path.append('../data/samples')
from os import path
from glob import glob
from bparkingdata_samples_2018 import bpark_samples_2018
from bparkingdata_samples_2022 import bpark_samples_2022
from bparkingdata_samples_2024 import bpark_samples_2024
from signal_samples_2018 import signal_samples_2018
from signal_samples_2024 import signal_samples_2024


def getOptions():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='Script to launch the nanoAOD tool on top of miniAOD files', add_help=True)
  parser.add_argument('--pl'      , type=str, dest='pl'          , help='label of the sample file'                                                            , default=None)
  parser.add_argument('--ds'      , type=str, dest='ds'          , help='[optional | data/sigcentral] run on specify data set. e.g "--ds D1"'                  , default=None)
  parser.add_argument('--year'    , type=str, dest='year'        , help='year to process in data'                                                             , default=None)
  parser.add_argument('--tagnano' , type=str, dest='tagnano'     , help='[optional] tag to be added on the outputfile name of the nano sample'                , default=None)
  parser.add_argument('--tagflat' , type=str, dest='tagflat'     , help='[optional] tag to be added on the outputfile name of the flat sample'                , default=None)
  parser.add_argument('--data'              , dest='data'        , help='run the nano tool on a data sample'                             , action='store_true', default=False)
  parser.add_argument('--mcprivate'         , dest='mcprivate'   , help='run the BParking nano tool on a private MC sample'              , action='store_true', default=False)
  parser.add_argument('--sigcentral'        , dest='sigcentral'  , help='run the BParking nano tool on a central MC sample'              , action='store_true', default=False)
  parser.add_argument('--donano'            , dest='donano'      , help='launch the nano tool on top of the minifile'                    , action='store_true', default=False)
  parser.add_argument('--doflat'            , dest='doflat'      , help='launch the ntupliser on top of the nanofile'                    , action='store_true', default=False)
  parser.add_argument('--dosignal'          , dest='dosignal'    , help='run the BToMuMuPi process'                                      , action='store_true', default=False)
  parser.add_argument('--dotageprobe'       , dest='dotageprobe' , help='run the JpsiToMuMu process (tag and probe study)'               , action='store_true', default=False)
  parser.add_argument('--dogeneral'         , dest='dogeneral'   , help='run without process'                                            , action='store_true', default=False)
  parser.add_argument('--docondor'          , dest='docondor'    , help='submit job on HTCondor'                                         , action='store_true', default=False)
  parser.add_argument('--dosubmit'          , dest='dosubmit'    , help='submit the jobs'                                                , action='store_true', default=False)
  return parser.parse_args()
  

def checkParser(opt):
  if opt.pl==None:
    raise RuntimeError('Please indicate the production label: for --mcprivate, it has to correspond to the label of the miniAOD')

  if opt.mcprivate==False and opt.sigcentral==False and opt.data==False:
    raise RuntimeError('Please indicate if you want to run on data or MC by adding either --data or --mcprivate or --sigcentral to the command line')

  if opt.donano==False and opt.doflat==False:
    raise RuntimeError('Please indicate if you want to run the nano tool (--donano) and/or the ntupliser (--doflat)')

  if opt.dosignal==False and opt.dotageprobe==False and opt.dogeneral==False:
    raise RuntimeError('Please indicate the process you want to run (--dosignal and/or --dotageprobe and/or --dogeneral)')


class NanoLauncher(object):
  def __init__(self, opt):
    self.prodlabel   = vars(opt)['pl']
    self.year        = vars(opt)['year']
    self.tagnano     = vars(opt)['tagnano']
    self.tagflat     = vars(opt)['tagflat']
    self.ds          = vars(opt)['ds']
    self.data        = vars(opt)['data']
    self.mcprivate   = vars(opt)['mcprivate']
    self.sigcentral  = vars(opt)['sigcentral']
    self.donano      = vars(opt)["donano"]
    self.doflat      = vars(opt)["doflat"]
    self.dosignal    = vars(opt)["dosignal"]
    self.dotageprobe = vars(opt)["dotageprobe"]
    self.dogeneral   = vars(opt)["dogeneral"]
    self.docondor    = vars(opt)['docondor']
    self.dosubmit    = vars(opt)["dosubmit"]

    self.user = os.environ["USER"]
    self.outfilename_nano = 'bparknano.root' if self.tagnano == None else 'bparknano_{}.root'.format(self.tagnano)

    if self.year not in ['2018', '2022', '2024']:
      raise RuntimeError("Wrong year '{}'. Please choose amongst ['2018', '2022', '2024']'".format(self.year))

    if self.data:
      if self.year == '2018':
        self.samples = bpark_samples_2018
      elif self.year == '2022':
        self.samples = bpark_samples_2022
      elif self.year == '2024':
        self.samples = bpark_samples_2024

      if self.ds == None:
        self.keys = self.samples.keys()
      else: # process a specific data set
        if self.ds not in self.samples.keys():
          raise RuntimeError('Please indicate on which period of the BParking dataset you want to run. Label "{}" not recognised. Choose among {}'.format(self.ds, self.samples.keys()))
        self.keys = [self.ds]

    if self.sigcentral:
      if self.year == '2018':
        self.samples = signal_samples_2018
      elif self.year == '2022':
        self.samples = signal_samples_2022
      elif self.year == '2024':
        self.samples = signal_samples_2024

      if self.ds == None:
        self.keys = self.samples.keys()
      else: # process a specific data set
        if self.ds not in self.samples.keys():
          raise RuntimeError('Please indicate on which signal sample you want to run. Label "{}" not recognised. Choose among {}'.format(self.ds, self.samples.keys()))
        self.keys = [self.ds]


  def get_chunk_id(self, chunk):
      idx = chunk.rfind('/')
      chunk_id = int(chunk[idx+1:len(chunk)])

      return chunk_id


  def create_flat_directory(self, chunk):
    outdir = chunk + '/flat'
    if not path.exists(outdir):
        os.makedirs(outdir)


  def prepare_nano_config(self, chunk=None, chunk_id=0, key=None):
    doSignal = 'True' if self.dosignal else 'False' 
    doTagAndProbe = 'True' if self.dotageprobe else 'False'
    doGeneral = 'True' if self.dogeneral else 'False'
    addTriggerMuonCollection = 'True'
    addProbeTracksCollection = 'True'

    if self.data:
      isMC = 'False'
      isMCCentral = 'False'
      isMCPrivate = 'False'
      if self.year == '2018':
        json_file = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'
        gt = '106X_dataRun2_v37'
        era = 'Run2_2018'
      elif self.year == '2022':
        json_file = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions22/Cert_Collisions2022_355100_362760_Golden.json'
        # see https://docs.google.com/presentation/d/1F4ndU7DBcyvrEEyLfYqb29NGkBPs20EAnBxe_l7AEII/edit?slide=id.g289f499aa6b_2_52#slide=id.g289f499aa6b_2_52
        if 'C' in key or 'D' in key: gt = '124X_dataRun3_PromptAnalysis_v1'
        elif 'E' in key: gt = '124X_dataRun3_Prompt_v10'
        if 'F' in key or 'G' in key: gt = '124X_dataRun3_PromptAnalysis_v2'
        era = 'Run3'
      elif self.year == '2024':
        json_file = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions24/Cert_Collisions2024_378981_386951_Golden.json'
        gt = '150X_dataRun3_v2'
        era = 'Run3_2024'
      else:
        raise RuntimeError('Please insert GT, json and era for year {}'.format(self.year))

    elif self.mcprivate or self.sigcentral:
      isMC = 'True'
      isMCCentral = 'True' if self.sigcentral else 'False'
      isMCPrivate = 'True' if self.mcprivate else 'False'
      if self.year == '2018':
        json_file = ''
        gt = '106X_upgrade2018_realistic_v16_L1v1'
        era = 'Run2_2018'
      # my understanding from https://docs.google.com/presentation/d/1F4ndU7DBcyvrEEyLfYqb29NGkBPs20EAnBxe_l7AEII/edit?slide=id.g289f499aa6b_2_52#slide=id.g289f499aa6b_2_52 is that preEE corresponds to data C/D and postEE E/F/G
      #elif self.year == '2022':
      #  json_file = ''
      #  gt = '130X_mcRun3_2022_realistic_v5' # preEE
      #  gt = '130X_mcRun3_2022_realistic_postEE_v6' # postEE
      #  era = 'Run3,run3_miniAOD_12X' #postEE
      elif self.year == '2024':
        json_file = ''
        gt = '150X_mcRun3_2024_realistic_v2'
        era = 'Run3_2024'
      else:
        raise RuntimeError('Please insert GT, json and era for year {}'.format(self.year))

    config = [
      "import FWCore.ParameterSet.Config as cms",
      "from glob import glob",
      "",
      "isMC = {}".format(isMC),
      "isMCCentral = {}".format(isMCCentral),
      "isMCPrivate = {}".format(isMCPrivate),
      "globaltag = '{}'".format(gt),
      "json_file = '{}'".format(json_file),
      "doSignal = {}".format(doSignal),
      "doTagAndProbe = {}".format(doTagAndProbe),
      "doGeneral = {}".format(doGeneral),
      "addTriggerMuonCollection = cms.untracked.bool({})".format(addTriggerMuonCollection),
      "addProbeTracksCollection = cms.untracked.bool({})".format(addProbeTracksCollection),
      "reportEvery = 1000",
      "",
      "if isMCPrivate:",
      "   inputFiles = ['file:%s' %i for i in glob('{}/step4_*.root')]".format(chunk),
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
      "    fileNames = cms.untracked.vstring() if not isMCPrivate else cms.untracked.vstring(inputFiles),",
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
      "    fileName = cms.untracked.string('file:{}'),".format(self.outfilename_nano),
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
      "process = nanoAOD_customizeTagAndProbeJPsiToMuMu (process, isMC=isMC)",
      "",
      "process.nanoAOD_general_step = cms.Path(process.nanoSequence)",
      "process.nanoAOD_BsToPhiPhiTo4K_step = cms.Path(process.nanoSequence + process.nanoBsToPhiPhiTo4KSequence + CountBsToPhiPhiTo4K)",
      "process.nanoAOD_JPsiToMuMu_step  = cms.Path(process.nanoSequence + process.nanoJPsiToMuMuSequence + CountJPsiToMuMu)",
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
      "if doTagAndProbe:",
      "  process.schedule += cms.Schedule(",
      "      process.nanoAOD_JPsiToMuMu_step,", 
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
      "if doTagAndProbe: process_string.append('nanoAOD_JPsiToMuMu_step')",
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

    if self.data or self.sigcentral:
        config_name = 'the_nano_config.py' 
    else:
        config_name = 'the_nano_config_{}.py'.format(chunk_id)

    f_out = open(config_name, 'w+')
    f_out.write(config)
    f_out.close()

    print('\t--> {} created'.format(config_name))


  def prepare_CRAB_config(self, key):
    dataset = self.samples[key] 

    label = '{}_{}'.format(self.prodlabel, key)

    if self.data:
      outdir = '/store/group/phys_bphys/anlyon/CPVGen/data/{a}/{b}/{c}'.format(
          a = self.prodlabel,
          b = self.year, 
          c = key,
          )
    elif self.sigcentral:
      outdir = '/store/group/phys_bphys/anlyon/CPVGen/signal_central/{a}/{b}/{c}'.format(
          a = self.prodlabel,
          b = self.year, 
          c = key,
          )
    else:
      raise RuntimeError('Please insert info for that data type')


    units_per_job = 4 # 3 # 10
    max_time_min = 130
    max_mem_MB = 2100

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
      "config.JobType.outputFiles     = ['{}']".format(self.outfilename_nano), 
      "config.JobType.maxMemoryMB = {}".format(max_mem_MB),
      "config.JobType.maxJobRuntimeMin = {}".format(max_time_min),

      "config.section_('Site')",
      "config.Site.storageSite = 'T2_CH_CERN'",
      ]

    config = '\n'.join(config)

    f_out = open('the_crab_config.py', 'w+')
    f_out.write(config)
    f_out.close()

    print('\t--> the_crab_config.py created')


  def prepare_dumper_starter(self, chunk, starter_name):
     infilename = 'bparknano_*.root' if self.tagnano == None else 'bparknano_{}_*.root'.format(self.tagnano) 
     list_files = [f for f in glob('{}/{}'.format(chunk, infilename))]

     outname = 'flat_bparknano'
     if self.tagnano != None: outname += '_' + self.tagnano
     if self.tagflat != None: outname += '_' + self.tagflat
     outfilename = '{}/flat/{}.root'.format(chunk, outname)

     #chunk_id = self.get_chunk_id(chunk=chunk)
     #if not path.exists('./condor'): os.makedirs('./condor')
     #starter_name = './condor/starter_{}_{}'.format(self.prodlabel, chunk_id)

     event_chain = []
     event_chain.append('TChain * c = new TChain("Events");')
     #for f in list_files:
     for ifile, f in enumerate(list_files):
       #if ifile > 1: continue #FIXME
       event_chain.append('  c->Add("{}");'.format(f))
     if self.dosignal: event_chain.append('  c->Process("BsToPhiPhiTo4KDumper.C+", outFileName);')
     event_chain = '\n'.join(event_chain)

     if self.data:
       addMC = ''
     else:
       addMC = 'outFileName += "_isSignalMC";'
  
     content = [
        '#include "TChain.h"',
        '#include <iostream>',
        'void starter(){',
        #'void {strn}'.format(strn=starter_name + '(){'),
           '  TString outFileName = "{outf}";'.format(outf=outfilename),
           '  {addMC}'.format(addMC = addMC),
           '  {chain}'.format(chain=event_chain),
           '}',
           ]
  
     content = '\n'.join(content)
  
     starter = './condor/{}.C'.format(starter_name)
     f_starter = open(starter, 'w+')
     f_starter.write(content)
     f_starter.close()
     print('\t--> {} created'.format(starter))
  
     #command_cp = 'cp {} {}/flat'.format(starter, path)
     #os.system(command_cp)
  
     #return starter


  def prepare_condor_config(self, chunk_id, submitter_name):

    logdir = './log/{}'.format(self.prodlabel)
    if self.tagnano != None: logdir += '_{}'.format(self.tagnano)
    if not path.exists(logdir):
        os.makedirs(logdir)

    #TODO implement option
    if self.donano:
        maxtime = 12 * 60 * 60 # in seconds
    elif self.doflat:
        maxtime = 12 * 60 * 60 # in seconds
    #maxtime = 0.3 * 60 * 60 # in seconds #FIXME

    config = [
       "universe              = vanilla",
       "executable            = ./condor/{}.sh".format(submitter_name),
       "mylogfile             = {}/job_{}_$(ClusterId)_$(ProcId).log".format(logdir, chunk_id),
       "log                   = $(mylogfile)",
       "output                = $(mylogfile)",
       "error                 = $(mylogfile)",
       "should_transfer_files = Yes",
       "use_x509userproxy     = true",
       "getenv                = True",
       "environment           = 'LS_SUBCWD={}'".format(os.getcwd()),
       "request_memory        = 2500",
       "+MaxRuntime           = {}".format(maxtime),
       "+AccountingGroup      = 'group_u_CMST3.all'",
       "queue",
      ]

    config = '\n'.join(config)

    f_out = open('the_condor_config.sub', 'w+')
    f_out.write(config)
    f_out.close()

    print('\t--> the_condor_config.sub created')


  def prepare_submitter(self, chunk_id):
                
    workdir = '/tmp/{}/{}_{}/'.format(self.user, self.prodlabel, chunk_id)

    submitter = [
      '#!/bin/bash',
      'cmsenv',
      'workdir="{}"'.format(workdir),
      'echo "creating workdir "$workdir',
      'mkdir -p $workdir',
      'echo "copying driver to workdir"',
      'cp {}/the_nano_config_{}.py $workdir'.format(os.getcwd(), chunk_id),
      'cd $workdir',
      'echo "going to run nano step"',
      'DATE_START=`date +%s`',
      'cmsRun {}/the_nano_config_{}.py'.format(workdir, chunk_id),
      #'cmsRun the_nano_config.py maxEvents=1000', #FIXME!
      'DATE_END=`date +%s`',
      'echo "finished running nano step"',
      'echo "content of the workdir"',
      'ls -l',
      'outdir="/eos/cms/store/group/phys_bphys/{}/CPVGen/{}/BsToPhiPhiTo4K/nanoFiles/Chunk{}"'.format(self.user, self.prodlabel, chunk_id),
      'echo "creating outdir "$outdir',
      'mkdir -p $outdir',
      'echo "copying the file"',
      'mv {} $outdir'.format(self.outfilename_nano),
      'cd $CMSSW_BASE/src/PhysicsTools/CPVNano/test',
      'echo "clearing the workdir"',
      'rm -r $workdir',
      'runtime=$((DATE_END-DATE_START))',
      'echo "Wallclock running time: $runtime s"',
      ]

    submitter = '\n'.join(submitter)

    f_out = open('the_submitter_{}.sh'.format(chunk_id), 'w+')
    f_out.write(submitter)
    f_out.close()

    print('\t--> the_submitter_{}.sh created'.format(chunk_id))


  def prepare_submitter_dumper(self, chunk_id, starter_name, submitter_name):
                
    #workdir = '/tmp/{}/dumper_{}_{}/'.format(self.user, self.prodlabel, chunk_id)
    workdir = '/tmp/{}/dumper_{}/'.format(self.user, starter_name)
    starter = './condor/' + starter_name + '.C'

    submitter = [
      '#!/bin/bash',
      'cmsenv',
      'startdir={}'.format(os.getcwd()),
      'workdir="{}"'.format(workdir),
      'echo "creating workdir "$workdir',
      'mkdir -p $workdir',
      'echo "copying scripts to workdir"',
      'mv $startdir/{} $workdir/starter.C'.format(starter),
      'cp $startdir/../plugins/dumper/utils.C $workdir',
      'cp $startdir/../plugins/dumper/BsToPhiPhiTo4KDumper.C $workdir',
      'cp $startdir/../plugins/dumper/BsToPhiPhiTo4KDumper.h $workdir',
      'cd $workdir',
      'echo "going to run the ntupliser"',
      'DATE_START_DUMP=`date +%s`',
      'root -l -q -b starter.C+', 
      'DATE_END_DUMP=`date +%s`',
      'echo "finished running the ntupliser"',
      'echo "content of the workdir"',
      'ls -l',
      'cd $startdir/',
      'echo "clearing the workdir"',
      'rm -r $workdir',
      'rm ./condor/{}.sh'.format(submitter_name),
      'runtime=$((DATE_END-DATE_START))',
      'echo "Wallclock running time: $runtime s"',
      ]

    submitter = '\n'.join(submitter)

    outname = './condor/' + submitter_name + '.sh'
    f_out = open(outname, 'w+')
    f_out.write(submitter)
    f_out.close()

    print('\t--> {} created'.format(outname))


  def submit_crab(self):
      command = 'crab submit -c the_crab_config.py'
      os.system(command)

      command_rm_nano = 'rm the_nano_config.py'
      os.system(command_rm_nano)

      command_rm_crab = 'rm the_crab_config.py'
      os.system(command_rm_crab)


  def submit_condor(self, job_name):
      command = 'condor_submit -batch-name {} the_condor_config.sub'.format(job_name)
      os.system(command)

      command_rm = 'rm the_condor_config.sub'
      os.system(command_rm)


  def submit_local(self):
      command = 'sh the_submitter.sh'
      os.system(command)

      command_rm_nano = 'rm the_nano_config.py'
      os.system(command_rm_nano)

      command_rm_submitter = 'rm the_submitter.sh'
      os.system(command_rm_submitter)


  def process(self):
    print('\n------------')
    print(' Processing NanoLauncher on production {} '.format(self.prodlabel))
    print('------------')

    if self.data:
      for key in self.keys:
        dataset = self.samples[key]
        print('\n-> Processing data set {}'.format(dataset)) 

        if self.donano:
            print('\n   -> Preparing nano config') 
            self.prepare_nano_config(key=key) 

            print('\n   -> Preparing CRAB config') 
            self.prepare_CRAB_config(key=key) 

            if self.dosubmit:
                print('\n   -> Submitting...') 
                self.submit_crab()

            print('\n   -> Submission completed')

        elif self.doflat:
            if not path.exists('./condor'): os.makedirs('./condor')

            path_dir = '/eos/cms/store/group/phys_bphys/{}/CPVGen/data/{}/{}/{}/*/*/*/*'.format(self.user, self.prodlabel, self.year, key)
            chunks = [f for f in glob(path_dir)]

            for chunk in chunks:
                chunk_id = self.get_chunk_id(chunk=chunk)
                print('\n   #-#-#- Chunk{} -#-#-#'.format(chunk_id)) 

                print('\n   -> Creating directory') 
                self.create_flat_directory(chunk=chunk)

                print('\n   -> Preparing dumper starter') 
                starter_name = 'starter_{}_{}_{}_{}'.format(self.year, self.prodlabel, key, chunk_id)
                if self.tagnano != None: starter_name += '_' + self.tagnano
                if self.tagflat != None: starter_name += '_' + self.tagflat

                self.prepare_dumper_starter(chunk=chunk, starter_name=starter_name)

                print('\n   -> Preparing dumper submitter') 
                submitter_name = 'the_submitter_dumper_{}_{}_{}_{}'.format(self.year, self.prodlabel, key, chunk_id)
                if self.tagnano != None: submitter_name += '_' + self.tagnano
                if self.tagflat != None: submitter_name += '_' + self.tagflat

                self.prepare_submitter_dumper(chunk_id=chunk_id, starter_name=starter_name, submitter_name=submitter_name)

                print('\n   -> Preparing condor config') 
                self.prepare_condor_config(chunk_id=chunk_id, submitter_name=submitter_name)

                if self.dosubmit:
                    print('\n   -> Submitting...') 
                    job_name = 'dumperstep_{}_{}_{}'.format(self.prodlabel, key, chunk_id) 
                    self.submit_condor(job_name=job_name)

    elif self.sigcentral:
      for key in self.keys:
        dataset = self.samples[key]
        print('\n-> Processing signal sample {}'.format(dataset)) 

        if self.donano:
            print('\n   -> Preparing nano config') 
            self.prepare_nano_config() 

            print('\n   -> Preparing CRAB config') 
            self.prepare_CRAB_config(key=key) 

            if self.dosubmit:
                print('\n   -> Submitting...') 
                self.submit_crab()

            print('\n   -> Submission completed')

        elif self.doflat:
            if not path.exists('./condor'): os.makedirs('./condor')

            path_dir = '/eos/cms/store/group/phys_bphys/{}/CPVGen/signal_central/{}/{}/{}/*/*/*/*'.format(self.user, self.prodlabel, self.year, key)
            chunks = [f for f in glob(path_dir)]

            for chunk in chunks:
                chunk_id = self.get_chunk_id(chunk=chunk)
                print('\n   #-#-#- Chunk{} -#-#-#'.format(chunk_id)) 

                print('\n   -> Creating directory') 
                self.create_flat_directory(chunk=chunk)

                print('\n   -> Preparing dumper starter') 
                starter_name = 'starter_{}_{}_{}_{}'.format(self.year, self.prodlabel, key, chunk_id)
                if self.tagnano != None: starter_name += '_' + self.tagnano
                if self.tagflat != None: starter_name += '_' + self.tagflat

                self.prepare_dumper_starter(chunk=chunk, starter_name=starter_name)

                print('\n   -> Preparing dumper submitter') 
                submitter_name = 'the_submitter_dumper_{}_{}_{}_{}'.format(self.year, self.prodlabel, key, chunk_id)
                if self.tagnano != None: submitter_name += '_' + self.tagnano
                if self.tagflat != None: submitter_name += '_' + self.tagflat

                self.prepare_submitter_dumper(chunk_id=chunk_id, starter_name=starter_name, submitter_name=submitter_name)

                print('\n   -> Preparing condor config') 
                self.prepare_condor_config(chunk_id=chunk_id, submitter_name=submitter_name)

                if self.dosubmit:
                    print('\n   -> Submitting...') 
                    job_name = 'dumperstep_{}_{}_{}_{}'.format(self.year, self.prodlabel, key, chunk_id) 
                    self.submit_condor(job_name=job_name)


    elif self.mcprivate:
        print('\n-> Processing private mc sample {}'.format(self.prodlabel)) 

        path_dir = '/eos/cms/store/group/phys_bphys/{}/CPVGen/{}/BsToPhiPhiTo4K/'.format(self.user, self.prodlabel)
        directories = [f for f in glob('{}/*'.format(path_dir)) if os.path.basename(f).startswith('crab')]

        for directory in directories:
            # get chunks
            chunks = [f for f in glob('{}/*/*'.format(directory))]

            for chunk in chunks:
                chunk_id = self.get_chunk_id(chunk=chunk)
                print('\n   #-#-#- Chunk{} -#-#-#'.format(chunk_id)) 

                if self.donano:
                    print('\n   -> Preparing nano config') 
                    self.prepare_nano_config(chunk=chunk, chunk_id=chunk_id) 

                    print('\n   -> Preparing submitter') 
                    self.prepare_submitter(chunk_id=chunk_id) 

                    if self.docondor:
                        print('\n   -> Preparing condor config') 
                        self.prepare_condor_config(chunk_id=chunk_id) #FIXME add submitter_name

                        if self.dosubmit:
                            print('\n   -> Submitting...') 
                            self.submit_condor() #FIXME add jobname

                    else:
                        if self.dosubmit:
                            print('\n   -> Processing...') 
                            self.submit_local()

    print('\nDone')



if __name__ == "__main__":
  opt = getOptions()
  checkParser(opt)
  NanoLauncher(opt).process()


  
