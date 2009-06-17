import FWCore.ParameterSet.Config as cms

#migration to tfile service at some point
TFileService = cms.Service("TFileService", fileName = cms.string('HLTbit.root') )

#associators
import SimTracker.TrackAssociation.TrackAssociatorByPosition_cfi
# associator
AssociatorByDeltaR0pt1 = SimTracker.TrackAssociation.TrackAssociatorByPosition_cfi.TrackAssociatorByPosition.clone(
    method = cms.string('posdr'),
    QCut = cms.double(0.1),
    ComponentName = cms.string('AssociatorByDeltaR0.1')
    )

AssociatorByDeltaR0pt2 = SimTracker.TrackAssociation.TrackAssociatorByPosition_cfi.TrackAssociatorByPosition.clone(
    method = cms.string('posdr'),
    QCut = cms.double(0.2),
    ComponentName = cms.string('AssociatorByDeltaR0.2')
    )

AssociatorByDeltaR1 = SimTracker.TrackAssociation.TrackAssociatorByPosition_cfi.TrackAssociatorByPosition.clone(
    method = cms.string('posdr'),
    QCut = cms.double(1.0),
    ComponentName = cms.string('AssociatorByDeltaR1.0')
    )



from Workspace.MuonHLTTreeUtility.muonHLTTreeUtility_cfi import *

#redo tracking particles
import SimGeneral.TrackingAnalysis.trackingParticles_cfi
mytrackingParticles = SimGeneral.TrackingAnalysis.trackingParticles_cfi.mergedtruth.clone()
mytrackingParticles.vertexDistanceCut = 1000

#muon tracking particles
import Validation.RecoTrack.cutsTPEffic_cfi
tpMuon = Validation.RecoTrack.cutsTPEffic_cfi.cutsTPEffic.clone()
#only muons
tpMuon.pdgId = cms.vint32(13,-13)
#allows decays in flight
tpMuon.tip = cms.double(10000.0)
tpMuon.lip = cms.double(10000.0)
tpMuon.src = cms.InputTag('mytrackingParticles')


#define the basic sequence of modules for muon HLT reconstruction
def muonHLTrecoSequence(process):
    #from HLTrigger.Configuration.HLT_2E30_cff import HLTBeginSequence,HLTL2muonrecoSequence,HLTL2muonisorecoSequence,HLTL3muonrecoSequence,HLTL3muonisorecoSequence
    muonHLTreco = cms.Sequence(process.HLTBeginSequence+process.HLTL2muonrecoSequence+process.HLTL2muonisorecoSequence+process.HLTL3muonrecoSequence+process.HLTL3muonisorecoSequence)
    return (muonHLTreco)

#redo the tracking particle
from Configuration.StandardSequences.MixingNoPileUp_cff import *
tpProduction = cms.Sequence ( mix * mytrackingParticles * tpMuon )
#redigi needed for association by hits only
tkSimDigiLinkAreThere = cms.EDFilter("IsProductAvailable",
                                     className = cms.string(''),
                                     src = cms.InputTag('')
                                     )
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

#from Validation.RecoMuon.associators_cff import tpToL3MuonAssociation,tpToL2MuonAssociation

#reDIGI_Path = cms.Paht( !tkSimDigiLinkAreThere + trdigi )


MHTU_Path = cms.Path( tpProduction ) 

#has to be in the endapth so that it can get the TriggerResults object
TimerService = cms.Service("TimerService",useCPUtime = cms.untracked.bool(True))
import HLTrigger.Timer.timer_cfi
hltTimer = HLTrigger.Timer.timer_cfi.myTimer.clone()
MHTU_EndPath = cms.EndPath( hltTimer * hltMuonTreeMaker )

#MHTUSchedule = cms.Schedule( reDIGI_Path + MHTU_Path )
MHTUSchedule = cms.Schedule( MHTU_Path )
#remember to extend with the EndPath

def insertMHTU(process):
    process.load('Workspace.MuonHLTTreeUtility.muonHLTTreeUtility_cff')
    process.muonHLTreco = muonHLTrecoSequence(process)
    process.MHTU_Path+=process.muonHLTreco
    import FWCore.ParameterSet.SequenceTypes
    for p in process.schedule:
        if (p.__class__==FWCore.ParameterSet.SequenceTypes.EndPath):
            process.schedule.insert(process.schedule.index(p), process.MHTU_Path )
            break
    process.schedule.append( process.MHTU_EndPath )

    ##actually do the --no_output option
    if (hasattr(process,"out_step")):
        process.schedule.remove(process.out_step)
