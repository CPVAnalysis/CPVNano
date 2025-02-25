import ROOT
import argparse
import numpy as np
from time import time
from datetime import datetime, timedelta
from array import array
from glob import glob
from collections import OrderedDict
from DataFormats.FWLite import Events, Handle
from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi, bestMatch
from itertools import product, combinations

handles = OrderedDict()
handles['muons'  ] = ('slimmedMuons'                 , Handle('std::vector<pat::Muon>')                   )
handles['trg_res'] = (('TriggerResults', '', 'HLT' ) , Handle('edm::TriggerResults'        )              )
handles['trg_ps' ] = (('patTrigger'    , '')         , Handle('pat::PackedTriggerPrescales')              )
handles['tobjs'  ] = ('slimmedPatTrigger'            , Handle('std::vector<pat::TriggerObjectStandAlone>'))

files = [
          '/scratch/anlyon/samples_tmp/data/4682963C-2EFF-FF4D-B234-8ED5973F70E4.root',
        ]

paths = ['HLT_Mu7_IP4', 'HLT_Mu8_IP6', 'HLT_Mu8_IP5', 'HLT_Mu8_IP3', 'HLT_Mu8p5_IP3p5', 'HLT_Mu9_IP6', 'HLT_Mu9_IP5', 'HLT_Mu9_IP4', 'HLT_Mu10p5_IP3p5', 'HLT_Mu12_IP6']
branches = [
  'HLT_Mu7_IP4',
  'HLT_Mu8_IP6',
  'HLT_Mu8_IP5',
  'HLT_Mu8_IP3',
  'HLT_Mu8p5_IP3p5',
  'HLT_Mu9_IP6',
  'HLT_Mu9_IP5',
  'HLT_Mu9_IP4',
  'HLT_Mu10p5_IP3p5',
  'HLT_Mu12_IP6',
  'HLT_Mu7_IP4_ps',
  'HLT_Mu8_IP6_ps',
  'HLT_Mu8_IP5_ps',
  'HLT_Mu8_IP3_ps',
  'HLT_Mu8p5_IP3p5_ps',
  'HLT_Mu9_IP6_ps',
  'HLT_Mu9_IP5_ps',
  'HLT_Mu9_IP4_ps',
  'HLT_Mu10p5_IP3p5_ps',
  'HLT_Mu12_IP6_ps',
  ]

tofill = OrderedDict(zip(branches, [-99.]*len(branches)))

events = Events(files)
maxevents = 5000 #-1

# start the stopwatch
start = time()

for i, event in enumerate(events):

    if (i+1) > maxevents:
        break
            
    if i%100 == 0:
        percentage = float(i) / maxevents * 100.
        speed = float(i) / (time() - start)
        eta = datetime.now() + timedelta(seconds=(maxevents-i) / max(0.1, speed))
        print('\t===> processing %d / %d event \t completed %.1f%s \t %.1f ev/s \t ETA %s s' %(i, maxevents, percentage, '%', speed, eta.strftime('%Y-%m-%d %H:%M:%S')))

    # access the handles
    for k, v in handles.iteritems():
        event.getByLabel(v[0], v[1])
        setattr(event, k, v[1].product())

    run  = event.eventAuxiliary().run()
    lumi = event.eventAuxiliary().luminosityBlock()
    iev  = event.eventAuxiliary().event()

    ######################################################################################
    #####      RECO PART HERE (GEN PART REMOVED FOR NOW)
    ######################################################################################

    # yeah, fire some trigger at least! For now, I've hard coded HLT_Mu7_IP4_part0
    trg_names = event.object().triggerNames(event.trg_res)

    hlt_passed = False

    for iname in trg_names.triggerNames():
        #if 'part0' not in iname: continue
        for ipath in paths:
            idx = len(trg_names)
            if iname.startswith(ipath):
                idx = trg_names.triggerIndex(iname)
                tofill[ipath        ] = ( idx < len(trg_names)) * (event.trg_res.accept(idx))
                tofill[ipath + '_ps'] = event.trg_ps.getPrescaleForIndex(idx)
                #if ipath=='HLT_Mu7_IP4' and event.trg_ps.getPrescaleForIndex(idx)>0 and ( idx < len(trg_names)) * (event.trg_res.accept(idx)):
                if ipath=='HLT_Mu12_IP6' and event.trg_ps.getPrescaleForIndex(idx)>0 and ( idx < len(trg_names)) * (event.trg_res.accept(idx)):
                    hlt_passed = True

    #print '{} {} {} {}'.format(run, lumi, iev, hlt_passed)

    if not hlt_passed:
        continue            

    # trigger matching
    # these are the filters, MAYBE!! too lazy to check confDB. Or, more appropriately: confDB sucks
    # https://github.com/cms-sw/cmssw/blob/6d2f66057131baacc2fcbdd203588c41c885b42c/Configuration/Skimming/python/pwdgSkimBPark_cfi.py#L11-L18 
    good_tobjs = []
    for to in [to for to in event.tobjs if to.pt()>6.5 and abs(to.eta())<2.]:
        to.unpackFilterLabels(event.object(), event.trg_res)
        #if to.hasFilterLabel('hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered7IP4Q'):
        if to.hasFilterLabel('hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered12Q'):
            good_tobjs.append(to)

    if len(good_tobjs) == 0: print len(good_tobjs)

    #muons = [mu for mu in event.muons if mu.pt()>1. and abs(mu.eta())<2.5 and mu.isPFMuon() and mu.isGlobalMuon()]
    muons = [mu for mu in event.muons if mu.pt()>1.5 and abs(mu.eta())<2.5]
    muons.sort(key = lambda x : x.pt(), reverse = True)

    #import pdb ; pdb.set_trace()


