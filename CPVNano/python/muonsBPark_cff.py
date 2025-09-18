import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

#FIXME add lines for 2022 and 2024
Path=["HLT_Mu7_IP4", "HLT_Mu8_IP6", "HLT_Mu8_IP5", "HLT_Mu8_IP3", "HLT_Mu8p5_IP3p5", "HLT_Mu9_IP6", "HLT_Mu9_IP5", "HLT_Mu9_IP4", "HLT_Mu10p5_IP3p5", "HLT_Mu12_IP6"]

muonTrgSelector = cms.EDProducer("MuonTriggerSelector",
                                 muonCollection = cms.InputTag("slimmedMuons"), #same collection as in NanoAOD                                                           
                                 displacedStandaloneMuonCollection = cms.InputTag("displacedStandAloneMuons"), #same collection as in NanoAOD                                                           
                                 vertexCollection = cms.InputTag('offlineSlimmedPrimaryVertices'), 

                                 # added for restoring first trigger matching method
                                 bits = cms.InputTag("TriggerResults", "", "HLT"), # format: module, label, process
                                 triggerObjects = cms.InputTag("slimmedPatTrigger"),
                                 prescales = cms.InputTag("patTrigger"),
                                 beamSpot = cms.InputTag("offlineBeamSpot"),

                                 # trigger muon matching conditions
                                 max_deltaR_trigger_matching = cms.double(0.03),
                                 max_deltaPtRel_trigger_matching = cms.double(0.1),
                                 
                                 # selection for the selected and trigger muon
                                 selmu_ptMin = cms.double(5.8),
                                 selmu_absEtaMax = cms.double(2.5), # increased for tag and probe study
                                 HLTPaths=cms.vstring(Path),
                             )
#cuts minimun number in B both mu and e, min number of trg, dz muon, dz and dr track, 
countTrgMuons = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("muonTrgSelector", "trgMuons")
)


muonBParkTable = cms.EDProducer("SimplePATMuonFlatTableProducer",
    src = cms.InputTag("muonTrgSelector:SelectedMuons"),
    cut = cms.string(""), #we should not filter on cross linked collections
    name = cms.string("Muon"),
    doc  = cms.string("slimmedMuons for BPark after basic selection"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(
        CandVars,
        ptErr = Var("bestTrack().ptError()", float, doc="ptError of the muon track", precision=12),
        vx = Var("vx()", float, doc="x coordinate of vertex position, in cm"),
        vy = Var("vy()", float, doc="y coordinate of vertex position, in cm"),
        vz = Var("vz()", float, doc="z coordinate of vertex position, in cm"),

        covQopQop   = Var("bestTrack().covariance(0, 0)", float, doc="Cov of q/p with q/p", precision=12),
        covLamLam   = Var("bestTrack().covariance(1, 1)", float, doc="Cov of lambda with lambda", precision=12),
        covPhiPhi   = Var("bestTrack().covariance(2, 2)", float, doc="Cov of phi with phi", precision=12),
        covQopLam   = Var("bestTrack().covariance(0, 1)", float, doc="Cov of q/p with lambda", precision=12),
        covQopPhi   = Var("bestTrack().covariance(0, 2)", float, doc="Cov of q/p with phi", precision=12),
        covLamPhi   = Var("bestTrack().covariance(1, 2)", float, doc="Cov of lambda with phi", precision=12),

        dz = Var("userFloat('dz')", float, doc="dz (with sign) wrt first PV, in cm", precision=10),
        dzS = Var("userFloat('dzS')", float, doc="dz (with sign) significance wrt first PV", precision=10),
        dz_alaRdst = Var("userFloat('dz_alaRdst')", float, doc="dz (innerTrack) (with sign) wrt first PV (a la R(D*)), in cm", precision=10),
        dxy = Var("userFloat('dxy')", float, doc="dxy (with sign) wrt first PV, in cm", precision=10),
        dxyS = Var("userFloat('dxyS')", float, doc="dxy (with sign) significance wrt first PV", precision=10),
        dxy_BS_alaRdst = Var("userFloat('dxy_BS_alaRdst')", float, doc="dxy (with sign) wrt first BS (a la R(D*)), in cm", precision=10),
        dxyS_BS_alaRdst = Var("userFloat('dxyS_BS_alaRdst')", float, doc="dxy (with sign) significance wrt BS (a la R(D*))", precision=10),
        dxy_BS = Var("userFloat('dxy_BS')", float, doc="dxy (with sign) wrt BS, in cm", precision=10),
        dxyS_BS = Var("userFloat('dxyS_BS')", float, doc="dxy (with sign) significance wrt BS", precision=10),
        ip3d = Var("userFloat('ip3d')", float, doc="3D impact parameter wrt first PV, in cm", precision=10),
        sip3d = Var("userFloat('sip3d')", float, doc="3D impact parameter significance wrt first PV", precision=10),

        looseID = Var("isLooseMuon()", int, doc="reco muon is Loose"),
        mediumID = Var("passed('CutBasedIdMedium')", bool, doc="cut-based ID, medium WP"),
        tightID = Var("passed('CutBasedIdTight')", bool, doc="cut-based ID, tight WP"),
        softID = Var("userInt('softID')", bool, doc="soft cut-based ID"),

        isPF = Var("isPFMuon()", int, doc="muon is PF candidate"),
        isGlobalMuon = Var("isGlobalMuon()", int, doc="muon is global muon"),
        isTrackerMuon = Var("isTrackerMuon()", int, doc="muon is tracker muon"),

        inTimeMuon = Var("passed('InTimeMuon')", bool, doc="inTimeMuon ID"),
        segmentCompatibility = Var("userFloat('segmentCompatibility')", float, doc = "muon segment compatibility: propagating the tracker tracks to the muon system and evaluate the number of matched segments and the closeness of the matching", precision=14), # keep higher precision since people have cuts with 3 digits on this
        caloCompatibility = Var("userFloat('caloCompatibility')", float, doc = "calorimetric compatibility"),
        validHitFraction = Var("userFloat('validHitFraction')", float, doc = "fraction of hits a tracker track uses (among inner tracker layers it traverses)"),
        kinkFinderChi2 = Var("userFloat('kinkFinderChi2')", float, doc = "chi2 of kink-finding algorithm: how likely it is that a track is made of more than one single track"),
        globalNormalisedChi2 = Var("userFloat('globalNormalisedChi2')", float, doc = "chi2/ndof of global fit"),
        localPositionChi2 = Var("userFloat('localPositionChi2')", float, doc = "chi2 of the position match between the tracker muon and the standalone muon"),
        trackerHighPurityFlag = Var("userInt('trackerHighPurityFlag')", int, doc = "tracker high-purity flag"), # int? 
        numberOfValidMuonHits = Var("userInt('numberOfValidMuonHits')", int, doc = "number of hits in the muon stations"),
        numberOfValidPixelHits = Var("userInt('numberOfValidPixelHits')", int, doc = "number of pixel hits"),
        numberOfTrackerLayers = Var("userInt('numberOfTrackerLayers')", int, doc = "number of tracker layers with hits"),
        numberOfPixelLayers = Var("userInt('numberOfPixelLayers')", int, doc = "number of pixel layers with hits"),
        numberOfStations = Var("userInt('nStations')", int, doc = "number of matched stations with default arbitration (segment & track)"),

        isTriggering = Var("userInt('isTriggering')", int, doc="flag the reco muon is also triggering"),
        isTriggeringBPark = Var("userInt('isTriggeringBPark')", int, doc="flag the reco muon is also triggering (only for BPark lines)"),
        matched_dr = Var("userFloat('DR')", float, doc="dr with the matched triggering muon" ),
        matched_dpt = Var("userFloat('DPT')", float, doc="dpt/pt with the matched triggering muon" ),   
        #TODO add Run 3 trigger lines
        fired_HLT_Mu7_IP4 = Var("userInt('HLT_Mu7_IP4')", int, doc="reco muon fired this trigger"),
        fired_HLT_Mu8_IP6 = Var("userInt('HLT_Mu8_IP6')", int, doc="reco muon fired this trigger"),
        fired_HLT_Mu8_IP5 = Var("userInt('HLT_Mu8_IP5')", int, doc="reco muon fired this trigger"),
        fired_HLT_Mu8_IP3 = Var("userInt('HLT_Mu8_IP3')", int, doc="reco muon fired this trigger"),
        fired_HLT_Mu8p5_IP3p5 = Var("userInt('HLT_Mu8p5_IP3p5')", int, doc="reco muon fired this trigger"),
        fired_HLT_Mu9_IP6 = Var("userInt('HLT_Mu9_IP6')", int, doc="reco muon fired this trigger"),
        fired_HLT_Mu9_IP5 = Var("userInt('HLT_Mu9_IP5')", int, doc="reco muon fired this trigger"),
        fired_HLT_Mu9_IP4 = Var("userInt('HLT_Mu9_IP4')", int, doc="reco muon fired this trigger"),
        fired_HLT_Mu10p5_IP3p5 = Var("userInt('HLT_Mu10p5_IP3p5')", int, doc="reco muon fired this trigger"),
        fired_HLT_Mu12_IP6 = Var("userInt('HLT_Mu12_IP6')", int, doc="reco muon fired this trigger"),
        prescale_HLT_Mu7_IP4 = Var("userInt('HLT_Mu7_IP4_prescale')", int, doc="reco muon prescale this trigger"),
        prescale_HLT_Mu8_IP6 = Var("userInt('HLT_Mu8_IP6_prescale')", int, doc="reco muon prescale this trigger"),
        prescale_HLT_Mu8_IP5 = Var("userInt('HLT_Mu8_IP5_prescale')", int, doc="reco muon prescale this trigger"),
        prescale_HLT_Mu8_IP3 = Var("userInt('HLT_Mu8_IP3_prescale')", int, doc="reco muon prescale this trigger"),
        prescale_HLT_Mu8p5_IP3p5 = Var("userInt('HLT_Mu8p5_IP3p5_prescale')", int, doc="reco muon prescale this trigger"),
        prescale_HLT_Mu9_IP6 = Var("userInt('HLT_Mu9_IP6_prescale')", int, doc="reco muon prescale this trigger"),
        prescale_HLT_Mu9_IP5 = Var("userInt('HLT_Mu9_IP5_prescale')", int, doc="reco muon prescale this trigger"),
        prescale_HLT_Mu9_IP4 = Var("userInt('HLT_Mu9_IP4_prescale')", int, doc="reco muon prescale this trigger"),
        prescale_HLT_Mu10p5_IP3p5 = Var("userInt('HLT_Mu10p5_IP3p5_prescale')", int, doc="reco muon prescale this trigger"),
        prescale_HLT_Mu12_IP6 = Var("userInt('HLT_Mu12_IP6_prescale')", int, doc="reco muon prescale this trigger"),
    ),
)

muonsBParkMCMatchForTable = cms.EDProducer("MCMatcher",       # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src         = muonBParkTable.src,                         # final reco collection
    matched     = cms.InputTag("finalGenParticlesBPark"),     # final mc-truth particle collection
    mcPdgId     = cms.vint32(13),                             # one or more PDG ID (13 = mu); absolute values (see below)
    checkCharge = cms.bool(False),                            # True = require RECO and MC objects to have the same charge
    mcStatus    = cms.vint32(1),                              # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR   = cms.double(0.25),                           # Minimum deltaR for the match
    maxDPtRel   = cms.double(0.25),                           # Minimum deltaPt/Pt for the match
    resolveAmbiguities    = cms.bool(False),                   # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True),                   # False = just match input in order; True = pick lowest deltaR pair first
    motherPdgId = cms.vint32(9900015, 443, 511, 521, 531, 541),
)

muonBParkMCTable = cms.EDProducer("CandMCMatchTableProducerBPark",
    src     = muonBParkTable.src,
    mcMap   = cms.InputTag("muonsBParkMCMatchForTable"),
    objName = muonBParkTable.name,
    objType = muonBParkTable.name, 
    branchName = cms.string("genPart"),
    docString = cms.string("MC matching to status==1 muons"),
)

selectedMuonsMCMatchEmbedded = cms.EDProducer(
    'MuonMatchEmbedder',
    src = cms.InputTag('muonTrgSelector:SelectedMuons'),
    matching = cms.InputTag('muonsBParkMCMatchForTable')
)


muonTriggerBParkTable = muonBParkTable.clone(
    src = cms.InputTag("muonTrgSelector:trgMuons"),
    name = cms.string("TriggerMuon"),
    doc  = cms.string("HLT Muons matched with reco muons"), #reco muon matched to triggering muon"),
    variables = cms.PSet(CandVars,
        vx = Var("vx()",float,doc="x coordinate of vertex position, in cm"),
        vy = Var("vy()",float,doc="y coordinate of vertex position, in cm"),
        vz = Var("vz()",float,doc="z coordinate of vertex position, in cm"),
        #trgMuonIndex = Var("userInt('trgMuonIndex')", int,doc="index in trigger muon collection"),
   )
)

# not used in the end
muonsTriggerBParkMCMatchForTable = cms.EDProducer("MCMatcher",# cut on deltaR, deltaPt/Pt; pick best by deltaR
    src         = muonTriggerBParkTable.src,                  # final reco collection
    matched     = cms.InputTag("finalGenParticlesBPark"),     # final mc-truth particle collection
    mcPdgId     = cms.vint32(13),                             # one or more PDG ID (13 = mu); absolute values (see below)
    checkCharge = cms.bool(False),                            # True = require RECO and MC objects to have the same charge
    mcStatus    = cms.vint32(1),                              # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR   = cms.double(0.1),                            # Minimum deltaR for the match
    maxDPtRel   = cms.double(0.25),                           # Minimum deltaPt/Pt for the match
    resolveAmbiguities    = cms.bool(True),                   # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True),                   # False = just match input in order; True = pick lowest deltaR pair first
    motherPdgId = cms.vint32(511, 521, 531, 541),
)

muonTriggerBParkMCTable = cms.EDProducer("CandMCMatchTableProducerBPark",
    src     = muonTriggerBParkTable.src,
    mcMap   = cms.InputTag("muonsTriggerBParkMCMatchForTable"),
    objName = muonTriggerBParkTable.name,
    objType = muonTriggerBParkTable.name, 
    branchName = cms.string("genPart"),
    docString = cms.string("MC matching to status==1 muons"),
)

triggerMuonsMCMatchEmbedded = cms.EDProducer(
    'TriggerMuonMatchEmbedder',
    src = cms.InputTag('muonTrgSelector', 'trgMuons'),
    matching = cms.InputTag('muonsTriggerBParkMCMatchForTable')
)

# saving events with trigger muon only
#muonBParkSequence = cms.Sequence(muonTrgSelector * countTrgMuons)
#muonBParkMC = cms.Sequence(muonTrgSelector + muonsBParkMCMatchForTable + selectedMuonsMCMatchEmbedded + muonBParkMCTable * countTrgMuons)
#muonBParkMCWithTriggerMuon = cms.Sequence(muonBParkMC + muonsTriggerBParkMCMatchForTable + triggerMuonsMCMatchEmbedded + muonTriggerBParkMCTable * countTrgMuons)

# saving events with or without trigger muon
muonBParkSequence = cms.Sequence(muonTrgSelector)
muonBParkMC = cms.Sequence(muonTrgSelector + muonsBParkMCMatchForTable + selectedMuonsMCMatchEmbedded + muonBParkMCTable)
muonBParkMCWithTriggerMuon = cms.Sequence(muonBParkMC + muonsTriggerBParkMCMatchForTable + triggerMuonsMCMatchEmbedded + muonTriggerBParkMCTable)

# saving table
muonBParkTables = cms.Sequence(muonBParkTable)
muonTriggerMatchedTables = cms.Sequence(muonTriggerBParkTable)
