import FWCore.ParameterSet.Config as cms


# need to redo the tracking particles
from Configuration.StandardSequences.MixingNoPileUp_cff import *
import Validation.RecoTrack.cutsTPEffic_cfi

#tracking particles
tpMuon = Validation.RecoTrack.cutsTPEffic_cfi.cutsTPEffic.clone()

# associator
import SimTracker.TrackAssociation.TrackAssociatorByPosition_cfi
from SimTracker.TrackAssociation.TrackAssociatorByHits_cfi import *
from SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi import *

TrackAssociatorByChi2ESProducer.beamSpot = cms.InputTag("hltOfflineBeamSpot")

assByHit = cms.string('TrackAssociatorByHits')
#assByHit = cms.string('TrackAssociatorByChi2')

AssociatorByDeltaR0pt1 = SimTracker.TrackAssociation.TrackAssociatorByPosition_cfi.TrackAssociatorByPosition.clone(propagator = cms.string('SteppingHelixPropagatorAlong'))
AssociatorByDeltaR0pt2 = SimTracker.TrackAssociation.TrackAssociatorByPosition_cfi.TrackAssociatorByPosition.clone(propagator = cms.string('SteppingHelixPropagatorAlong'))


DQMStore = cms.Service("DQMStore")

from RecoMuon.MuonXRay.IDconverttoBinNum_cfi import IDconverttoBinNum
from RecoMuon.MuonXRay.DQMHelper_muonXRay_cfi import DQMHelper_muonXRay

muonXray = cms.EDFilter("MuonXRay",
    DQMHelper = DQMHelper_muonXRay,
    IDconverttoBinNum = IDconverttoBinNum
)

muonXrayAtL1 = cms.EDFilter("MuonXRay",
    DQMHelper = DQMHelper_muonXRay,
    IDconverttoBinNum = IDconverttoBinNum
)

muonXrayAtL2 = cms.EDFilter("MuonXRay",
    DQMHelper = DQMHelper_muonXRay,
    IDconverttoBinNum = IDconverttoBinNum
)

muonXrayAtL2iso = cms.EDFilter("MuonXRay",
    DQMHelper = DQMHelper_muonXRay,
    IDconverttoBinNum = IDconverttoBinNum
)

muonXrayAtL3 = cms.EDFilter("MuonXRay",
    DQMHelper = DQMHelper_muonXRay,
    IDconverttoBinNum = IDconverttoBinNum
)

muonXrayAtL3iso = cms.EDFilter("MuonXRay",
    DQMHelper = DQMHelper_muonXRay,
    IDconverttoBinNum = IDconverttoBinNum
)


#mytrackingtruthprod = SimGeneral.TrackingAnalysis.trackingParticles_cfi.mergedtruth.clone()
mytrackingtruthprod = cms.EDProducer("TrackingTruthProducer",
    discardOutVolume = cms.bool(False),
    DiscardHitsFromDeltas = cms.bool(True),
    simHitLabel = cms.string('g4SimHits'),
    volumeRadius = cms.double(1200.0),
    vertexDistanceCut = cms.double(0.003),
                                     mergedBremsstrahlung = cms.bool(True),
    HepMCDataLabels = cms.vstring('VtxSmeared', 
        'PythiaSource', 
        'source'),
    TrackerHitLabels = cms.vstring(
                                   'g4SimHitsTrackerHitsPixelBarrelLowTof',
                                   'g4SimHitsTrackerHitsPixelBarrelHighTof',
                                   'g4SimHitsTrackerHitsPixelEndcapLowTof',
                                   'g4SimHitsTrackerHitsPixelEndcapHighTof',
                                   'g4SimHitsTrackerHitsTIBLowTof',
                                   'g4SimHitsTrackerHitsTIBHighTof',
                                   'g4SimHitsTrackerHitsTIDLowTof',
                                   'g4SimHitsTrackerHitsTIDHighTof',
                                   'g4SimHitsTrackerHitsTOBLowTof',
                                   'g4SimHitsTrackerHitsTOBHighTof',
                                   'g4SimHitsTrackerHitsTECLowTof',
                                   'g4SimHitsTrackerHitsTECHighTof'
                                   ),
    volumeZ = cms.double(3000.0)
)

# add the muon hits to the trackingParticle
mytrackingtruthprodWithmuHits = cms.EDProducer("TrackingTruthProducer",
    discardOutVolume = cms.bool(False),
    DiscardHitsFromDeltas = cms.bool(True),
    simHitLabel = cms.string('g4SimHits'),
    volumeRadius = cms.double(1200.0),
    vertexDistanceCut = cms.double(0.003),
                                               mergedBremsstrahlung = cms.bool(True),
    HepMCDataLabels = cms.vstring('VtxSmeared', 
        'PythiaSource', 
        'source'),
    TrackerHitLabels = cms.vstring(
                                   'g4SimHitsTrackerHitsPixelBarrelLowTof',
                                   'g4SimHitsTrackerHitsPixelBarrelHighTof',
                                   'g4SimHitsTrackerHitsPixelEndcapLowTof',
                                   'g4SimHitsTrackerHitsPixelEndcapHighTof',
                                   'g4SimHitsTrackerHitsTIBLowTof',
                                   'g4SimHitsTrackerHitsTIBHighTof',
                                   'g4SimHitsTrackerHitsTIDLowTof',
                                   'g4SimHitsTrackerHitsTIDHighTof',
                                   'g4SimHitsTrackerHitsTOBLowTof',
                                   'g4SimHitsTrackerHitsTOBHighTof',
                                   'g4SimHitsTrackerHitsTECLowTof',
                                   'g4SimHitsTrackerHitsTECHighTof',
                                   'g4SimHitsMuonDTHits',
                                   'g4SimHitsMuonCSCHits',
                                   'g4SimHitsMuonRPCHits'
                                   ),
    volumeZ = cms.double(3000.0)
)

# add the muon hits to the trackingParticle
#replace mytrackingtruthprodWithmuHits.MuonHitLabels +={ g4SimHits:MuonDTHits, g4SimHits:MuonCSCHits, g4SimHits:MuonRPCHits }
mytrackingtruthprodTEST = cms.EDFilter("TrackingTruthTest")

from RecoMuon.MuonXRay.DQMHelper_associatedMuonXRay_cfi import DQMHelper_associatedMuonXRay

associatedMuonXrayAtL2 = cms.EDFilter("AssociatedMuonXRay",
    IDconverttoBinNum = IDconverttoBinNum,
    DQMHelper = DQMHelper_associatedMuonXRay,                                     
    trackingParticleLabel = cms.InputTag("tpMuon"),
    associatorName = cms.string('AssociatorByDeltaR0.2'),
    trackLabel = cms.InputTag("hltL2Muons","UpdatedAtVtx")
)

associatedMuonXrayAtL2noUpdate = cms.EDFilter("AssociatedMuonXRay",
    IDconverttoBinNum = IDconverttoBinNum,
    DQMHelper = DQMHelper_associatedMuonXRay,
    trackingParticleLabel = cms.InputTag("tpMuon"),
    associatorName = cms.string('AssociatorByDeltaR0.2'),
    trackLabel = cms.InputTag("hltL2Muons")
)

associatedAnyXrayAtL2 = cms.EDFilter("AssociatedMuonXRay",
    IDconverttoBinNum = IDconverttoBinNum,
    DQMHelper = DQMHelper_associatedMuonXRay,
    trackingParticleLabel = cms.InputTag("mytrackingtruthprodWithmuHits"),
    associatorName = cms.string('AssociatorByDeltaR0.2'),
    trackLabel = cms.InputTag("hltL2Muons","UpdatedAtVtx"),
    isolationLabel = cms.InputTag("hltL2MuonIsolations")
)

associatedAnyXrayAtL2Cand = cms.EDFilter("AssociatedMuonXRay",
    associatorName = cms.string('AssociatorByDeltaR0.2'),
    confirmLabel = cms.InputTag("hltSingleMuNoIsoL2PreFiltered"),
    IDconverttoBinNum = IDconverttoBinNum,
    DQMHelper = DQMHelper_associatedMuonXRay,
    trackLabel = cms.InputTag("hltL2Muons","UpdatedAtVtx"),
    trackingParticleLabel = cms.InputTag("mytrackingtruthprodWithmuHits")
)

associatedMuonXrayAtL3 = cms.EDFilter("AssociatedMuonXRay",
    IDconverttoBinNum = IDconverttoBinNum,
    DQMHelper = DQMHelper_associatedMuonXRay,
    trackingParticleLabel = cms.InputTag("tpMuon"),
    associatorName = cms.string('AssociatorByDeltaR0.1'),
    trackLabel = cms.InputTag("hltL3Muons"),
    isolationLabel = cms.InputTag("hltL3MuonIsolations")                                      
)

associatedAnyXrayAtL3 = cms.EDFilter("AssociatedMuonXRay",
    IDconverttoBinNum = IDconverttoBinNum,
    DQMHelper = DQMHelper_associatedMuonXRay,
    trackingParticleLabel = cms.InputTag("mytrackingtruthprodWithmuHits"),
    associatorName = cms.string('AssociatorByDeltaR0.1'),
    trackLabel = cms.InputTag("hltL3Muons")
)

associatedAnyXrayAtL3Cand = cms.EDFilter("AssociatedMuonXRay",
    associatorName = cms.string('AssociatorByDeltaR0.1'),
    confirmLabel = cms.InputTag("hltSingleMuNoIsoL3PreFiltered"),
    IDconverttoBinNum = IDconverttoBinNum,
    DQMHelper = DQMHelper_associatedMuonXRay,
    trackLabel = cms.InputTag("hltL3Muons"),
    trackingParticleLabel = cms.InputTag("mytrackingtruthprodWithmuHits")
)

associatedAnyXrayAtL3byHits = cms.EDFilter("AssociatedMuonXRay",
    IDconverttoBinNum = IDconverttoBinNum,
    DQMHelper = DQMHelper_associatedMuonXRay,
    trackingParticleLabel = cms.InputTag("mytrackingtruthprod"),
    associatorName = assByHit,
    trackLabel = cms.InputTag("hltL3Muons","L2Seeded")
)

associatedAnyXrayAtL3CandbyHits = cms.EDFilter("AssociatedMuonXRay",
    associatorName = assByHit,
    confirmLabel = cms.InputTag("hltSingleMuNoIsoL3PreFiltered"),
    IDconverttoBinNum = IDconverttoBinNum,
    linkLabel = cms.InputTag("hltL3Muons"),
    DQMHelper = DQMHelper_associatedMuonXRay,
    trackLabel = cms.InputTag("hltL3Muons","L2Seeded"),
    trackingParticleLabel = cms.InputTag("mytrackingtruthprod")
)

dqmSaver = cms.EDAnalyzer("DQMRootFile",
                          DQMRootFileName = cms.string('DQM.root')
                          )

#dqmSaver = cms.EDFilter("DQMFileSaver",
#    prescaleEvt = cms.untracked.int32(-1),
#    workflow = cms.untracked.string('/A/B/C'),
#    forceRunNumber = cms.untracked.int32(1),
#    prescaleLS = cms.untracked.int32(-1),
#    saveAtJobEnd = cms.untracked.bool(True),
#    fileName = cms.untracked.string('dqmRootFile'),
#    #        untracked string environment = "Online"
#    environment = cms.untracked.string('Offline'),
#    saveAtRunEnd = cms.untracked.bool(False),
#    prescaleTime = cms.untracked.int32(-1)
#)
#


#redo the DigiLinks too
from SimTracker.Configuration.SimTracker_cff import *
#from IOMC.RandomEngine.IOMC_cff import *
RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                           restoreStateLabel = cms.untracked.string('randomEngineStateProducer'),
                                           simSiPixelDigis = cms.PSet(
    initialSeed = cms.untracked.uint32(1234567),
    engineName = cms.untracked.string('HepJamesRandom')
    ),
                                           simSiStripDigis = cms.PSet(
    initialSeed = cms.untracked.uint32(1234567),
    engineName = cms.untracked.string('HepJamesRandom')
    )
                                           )

simParticle = cms.Sequence(mix*mytrackingtruthprod*mytrackingtruthprodWithmuHits*tpMuon*trDigi)


from Workspace.MuonHLTTreeUtility.muonHLTTreeUtility_cfi import hltMuonTreeMaker
treePath = cms.Path(simParticle)
treeEndPath = cms.EndPath( hltMuonTreeMaker)

xRay = cms.Path(muonXray)
#include "HLTrigger/Muon/data/PathSingleMu_1032_NoIso.cff"

from HLTrigger.Configuration.HLT_2E30_cff import hltL1sSingleMuNoIso10,hltSingleMuNoIsoL1Filtered10,hltSingleMuNoIsoL2PreFiltered11,hltSingleMuNoIsoL3PreFiltered13
from HLTrigger.Configuration.HLT_2E30_cff import hltL1sSingleMuIso7,hltSingleMuIsoL1Filtered,hltSingleMuIsoL2PreFiltered7,hltSingleMuIsoL2IsoFiltered7,hltSingleMuIsoL3PreFiltered9,hltSingleMuIsoL3IsoFiltered9

xRayL1 = cms.Path(hltL1sSingleMuNoIso10+hltSingleMuNoIsoL1Filtered10+muonXrayAtL1)
xRayL2 = cms.Path(hltL1sSingleMuNoIso10+hltSingleMuNoIsoL1Filtered10+hltSingleMuNoIsoL2PreFiltered11+muonXrayAtL2)
xRayL2iso = cms.Path(hltL1sSingleMuIso7+hltSingleMuIsoL1Filtered+hltSingleMuIsoL2PreFiltered7+hltSingleMuIsoL2IsoFiltered7+muonXrayAtL2iso)
xRayL3 = cms.Path(hltL1sSingleMuNoIso10+hltSingleMuNoIsoL1Filtered10+hltSingleMuNoIsoL2PreFiltered11+hltSingleMuNoIsoL3PreFiltered13+muonXrayAtL3)
xRayL3iso = cms.Path(hltL1sSingleMuIso7+hltSingleMuIsoL1Filtered+hltSingleMuIsoL2PreFiltered7+hltSingleMuIsoL2IsoFiltered7+hltSingleMuIsoL3PreFiltered9+hltSingleMuIsoL3IsoFiltered9+muonXrayAtL3)

l2Xray = cms.Path(simParticle+hltL1sSingleMuNoIso10+hltSingleMuNoIsoL1Filtered10+associatedMuonXrayAtL2+associatedAnyXrayAtL2+associatedMuonXrayAtL2noUpdate+hltSingleMuNoIsoL2PreFiltered11+associatedAnyXrayAtL2Cand)
#l2XrayIso = cms.Path(simParticle+hltL1sSingleMuIso7+hltSingleMuIsoL1Filtered+  +hltSingleMuIsoL2PreFiltered7+associatedAnyXrayAtL2Cand)

l3Xray = cms.Path(simParticle+hltL1sSingleMuNoIso10+hltSingleMuNoIsoL1Filtered10+hltSingleMuNoIsoL2PreFiltered11+associatedMuonXrayAtL3+associatedAnyXrayAtL3+associatedAnyXrayAtL3byHits+hltSingleMuNoIsoL3PreFiltered13+associatedAnyXrayAtL3Cand+associatedAnyXrayAtL3CandbyHits)
#l3XrayIso = cms.Path(simParticle+hltL1sSingleMuIso7+hltSingleMuIsoL1Filtered+hltSingleMuIsoL2PreFiltered7+associatedMuonXrayAtL3+associatedAnyXrayAtL3+associatedAnyXrayAtL3byHits+hltSingleMuIsoL3PreFiltered9+associatedAnyXrayAtL3Cand+associatedAnyXrayAtL3CandbyHits)

# 
# include "HLTrigger/Muon/data/PathSingleMu_1032_Iso.cff"
# path xRayL1 = { hltSingleMuIsoLevel1 & hltSingleMuIsoL1Filtered & muonXrayAtL1 }
# path xRayL2 = { hltSingleMuIsoL2PreFiltered7 & muonXrayAtL2 }
# path xRayL2iso = { hltSingleMuIsoL2PreFiltered7 & muonXrayAtL2iso }
# path xRayL3 = { hltSingleMuIsoL3PreFiltered9 & muonXrayAtL3 }
# path xRayL3iso = { hltSingleMuIsoL3PreFiltered9 & muonXrayAtL3iso }
#
saver = cms.EndPath(dqmSaver)

from HLTrigger.Configuration.HLT_2E30_cff import HLTBeginSequence,HLTL2muonrecoSequence,HLTL2muonisorecoSequence,HLTL3muonrecoSequence,HLTL3muonisorecoSequence

pReco = cms.Path(HLTBeginSequence+HLTL2muonrecoSequence+HLTL2muonisorecoSequence+HLTL3muonrecoSequence+HLTL3muonisorecoSequence)

#allow decays in flight to be there
mytrackingtruthprod.vertexDistanceCut = 1000
#allow decays in flight to be there
mytrackingtruthprodWithmuHits.vertexDistanceCut = 1000
#only muons
tpMuon.pdgId = [13, -13]
#allows decays in flight
tpMuon.tip = 10000
tpMuon.lip = 10000
tpMuon.src = 'mytrackingtruthprodWithmuHits'

AssociatorByDeltaR0pt1.method = 'posdr'
AssociatorByDeltaR0pt1.QCut = 0.1
AssociatorByDeltaR0pt1.ComponentName = 'AssociatorByDeltaR0.1'

AssociatorByDeltaR0pt2.method = 'posdr'
AssociatorByDeltaR0pt2.QCut = 0.2
AssociatorByDeltaR0pt2.ComponentName = 'AssociatorByDeltaR0.2'

associatedMuonXrayAtL2.DQMHelper.H1s.h_mu_d0.xAxis.Min = -10
associatedMuonXrayAtL2.DQMHelper.H1s.h_mu_d0.xAxis.Max = 10
associatedMuonXrayAtL2noUpdate.DQMHelper.H1s.h_mu_d0.xAxis.Min = -10
associatedMuonXrayAtL2noUpdate.DQMHelper.H1s.h_mu_d0.xAxis.Max = 10
associatedAnyXrayAtL2.DQMHelper.H1s.h_mu_d0.xAxis.Min = -10
associatedAnyXrayAtL2.DQMHelper.H1s.h_mu_d0.xAxis.Max = 10
associatedAnyXrayAtL2Cand.DQMHelper.H1s.h_mu_d0.xAxis.Min = -10
associatedAnyXrayAtL2Cand.DQMHelper.H1s.h_mu_d0.xAxis.Max = 10


