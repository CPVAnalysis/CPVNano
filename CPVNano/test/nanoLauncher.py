import os
import sys
sys.path.append('../data/samples')
from os import path
from glob import glob
from bparkingdata_samples_2018 import bpark_samples_2018
from bparkingdata_samples_2022 import bpark_samples_2022
from bparkingdata_samples_2024 import bpark_samples_2024


def getOptions():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='Script to launch the nanoAOD tool on top of miniAOD files', add_help=True)
  parser.add_argument('--pl'      , type=str, dest='pl'          , help='label of the sample file'                                                            , default=None)
  parser.add_argument('--ds'      , type=str, dest='ds'          , help='[optional | data/mccentral] run on specify data set. e.g "--ds D1"'                  , default=None)
  parser.add_argument('--year'    , type=str, dest='year'        , help='year to process in data'                                                             , default=None)
  parser.add_argument('--tagnano' , type=str, dest='tagnano'     , help='[optional] tag to be added on the outputfile name of the nano sample'                , default=None)
  parser.add_argument('--data'              , dest='data'        , help='run the nano tool on a data sample'                             , action='store_true', default=False)
  parser.add_argument('--mcprivate'         , dest='mcprivate'   , help='run the BParking nano tool on a private MC sample'              , action='store_true', default=False)
  parser.add_argument('--dosignal'          , dest='dosignal'    , help='run the BToMuMuPi process'                                      , action='store_true', default=False)
  parser.add_argument('--dotageprobe'       , dest='dotageprobe' , help='run the JpsiToMuMu process (tag and probe study)'               , action='store_true', default=False)
  parser.add_argument('--dogeneral'         , dest='dogeneral'   , help='run without process'                                            , action='store_true', default=False)
  parser.add_argument('--docondor'          , dest='docondor'    , help='submit job on HTCondor'                                         , action='store_true', default=False)
  return parser.parse_args()
  

def checkParser(opt):
  if opt.pl==None:
    raise RuntimeError('Please indicate the production label: for --mcprivate, it has to correspond to the label of the miniAOD')

  if opt.dosignal==False and opt.dotageprobe==False and opt.dogeneral==False:
    raise RuntimeError('Please indicate the process you want to run (--dosignal and/or --dotageprobe and/or --dogeneral)')


class NanoLauncher(object):
  def __init__(self, opt):
    self.prodlabel = vars(opt)['pl']
    self.year      = vars(opt)['year']
    self.tagnano   = vars(opt)['tagnano']
    self.ds        = vars(opt)['ds']
    self.data      = vars(opt)['data']
    self.mcprivate = vars(opt)['mcprivate']
    self.dosignal    = vars(opt)["dosignal"]
    self.dotageprobe = vars(opt)["dotageprobe"]
    self.dogeneral   = vars(opt)["dogeneral"]
    self.docondor  = vars(opt)['docondor']

    self.user = os.environ["USER"]
    self.outfilename_nano = 'bparknano.root' if self.tagnano == None else 'bparknano_{}.root'.format(self.tagnano)

    if self.year not in ['2018', '2022', '2024']:
      raise RuntimeError("Wrong year '{}'. Please choose amongst ['2018', '2022', '2024']'".format(self.year))

    if self.data:
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


  def get_chunk_id(self, chunk):
      idx = chunk.rfind('/')
      chunk_id = int(chunk[idx+1:len(chunk)])

      return chunk_id


  def prepare_nano_config(self, chunk=None, chunk_id=0):
    doSignal = 'True' if self.dosignal else 'False' 
    doTagAndProbe = 'True' if self.dotageprobe else 'False'
    doGeneral = 'True' if self.dogeneral else 'False'
    addTriggerMuonCollection = 'True'
    addProbeTracksCollection = 'True'

    if self.data:
      isMC = 'False'
      if self.year == '2018':
        json_file = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'
        gt = '102X_dataRun2_v11'
        era = 'Run2_2018'
      elif self.year == '2022':
        json_file = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions22/Cert_Collisions2022_355100_362760_Golden.json'
        gt = '' # see https://docs.google.com/presentation/d/1F4ndU7DBcyvrEEyLfYqb29NGkBPs20EAnBxe_l7AEII/edit?slide=id.g289f499aa6b_2_52#slide=id.g289f499aa6b_2_52
        era = ''
      elif self.year == '2024':
        json_file = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions24/Cert_Collisions2024_378981_386951_Golden.json'
        gt = '150X_dataRun3_v2'
        era = 'Run3_2024'
      else:
        raise RuntimeError('Please insert GT, json and era for year {}'.format(self.year))

    elif self.mcprivate:
      isMC = 'True'
      if self.year == '2018':
        json_file = ''
        gt = '102X_upgrade2018_realistic_v15'
        era = 'Run2_2018'
      else:
        raise RuntimeError('Please insert GT, json and era for year {}'.format(self.year))

    config = [
      "import FWCore.ParameterSet.Config as cms",
      "from glob import glob",
      "",
      "isMC = {}".format(isMC),
      "globaltag = '{}'".format(gt),
      "json_file = '{}'".format(json_file),
      "doSignal = {}".format(doSignal),
      "doTagAndProbe = {}".format(doTagAndProbe),
      "doGeneral = {}".format(doGeneral),
      "addTriggerMuonCollection = cms.untracked.bool({})".format(addTriggerMuonCollection),
      "addProbeTracksCollection = cms.untracked.bool({})".format(addProbeTracksCollection),
      "reportEvery = 1000",
      "",
      "if isMC:",
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
      "    fileNames = cms.untracked.vstring() if not isMC else cms.untracked.vstring(inputFiles),",
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

    if self.data:
        config_name = 'the_nano_config.py' 
    else:
        config_name = 'the_nano_config_{}.py'.format(chunk_id)

    f_out = open(config_name, 'w+')
    f_out.write(config)
    f_out.close()

    print('\t--> {} created'.format(config_name))


  def prepare_CRAB_config(self, key):
    dataset = self.bpark_samples[key] 

    label = '{}_{}'.format(self.prodlabel, key)

    if self.data:
      outdir = '/store/group/phys_bphys/anlyon/CPVGen/data/{a}/{b}/{c}'.format(
          a = self.prodlabel,
          b = self.year, 
          c = key,
          )

    units_per_job = 3 # 10

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


  def prepare_condor_config(self, chunk_id):

    logdir = './log/{}'.format(self.prodlabel)
    if self.tagnano != None: logdir += '_{}'.format(self.tagnano)
    if not path.exists(logdir):
        os.makedirs(logdir)

    #TODO implement option
    maxtime = 12 * 60 * 60 # in seconds
    #maxtime = 0.3 * 60 * 60 # in seconds #FIXME

    config = [
       "universe              = vanilla",
       "executable            = the_submitter_{}.sh".format(chunk_id),
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


  def submit_crab(self):
      command = 'crab submit -c the_crab_config.py'
      os.system(command)

      command_rm_nano = 'rm the_nano_config.py'
      os.system(command_rm_nano)

      command_rm_crab = 'rm the_crab_config.py'
      os.system(command_rm_crab)


  def submit_condor(self):
      command = 'condor_submit -batch-name nano_{} the_condor_config.sub'.format(self.prodlabel)
      os.system(command)


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
        dataset = self.bpark_samples[key]
        print('\n-> Processing data set {}'.format(dataset)) 

        print('\n   -> Preparing nano config') 
        self.prepare_nano_config() 

        print('\n   -> Preparing CRAB config') 
        self.prepare_CRAB_config(key=key) 

        print('\n   -> Submitting...') 
        self.submit_crab()

        print('\n   -> Submission completed')


    if self.mcprivate:
        print('\n-> Processing private mc sample {}'.format(self.prodlabel)) 

        path_dir = '/eos/cms/store/group/phys_bphys/{}/CPVGen/{}/BsToPhiPhiTo4K/'.format(self.user, self.prodlabel)
        directories = [f for f in glob('{}/*'.format(path_dir)) if os.path.basename(f).startswith('crab')]

        for directory in directories:
            # get chunks
            chunks = [f for f in glob('{}/*/*'.format(directory))]

            for chunk in chunks:
                chunk_id = self.get_chunk_id(chunk=chunk)
                print('\n   #-#-#- Chunk{} -#-#-#'.format(chunk_id)) 

                print('\n   -> Preparing nano config') 
                self.prepare_nano_config(chunk=chunk, chunk_id=chunk_id) 

                print('\n   -> Preparing submitter') 
                self.prepare_submitter(chunk_id=chunk_id) 

                if self.docondor:
                    print('\n   -> Preparing condor config') 
                    self.prepare_condor_config(chunk_id=chunk_id)

                    print('\n   -> Submitting...') 
                    self.submit_condor()

                else:
                    print('\n   -> Processing...') 
                    self.submit_local()

    print('\nDone')



if __name__ == "__main__":
  opt = getOptions()
  checkParser(opt)
  NanoLauncher(opt).process()


  
