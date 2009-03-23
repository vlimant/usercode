import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'ERROR'
process.MessageLogger.categories.append('PATLayer0Summary')
process.MessageLogger.categories.append('Selections')
process.MessageLogger.infos = cms.untracked.PSet(
    threshold        = cms.untracked.string( 'INFO' ),
    INFO= cms.untracked.PSet( limit = cms.untracked.int32(0)  ),
    default          = cms.untracked.PSet( limit = cms.untracked.int32(0)  ),
    PATLayer0Summary = cms.untracked.PSet( limit = cms.untracked.int32(-1) ),
    Selections = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# source
process.source = cms.Source("PoolSource", 
     fileNames = cms.untracked.vstring(
   # 'file:/afs/cern.ch/cms/PRS/top/cmssw-data/relval200-for-pat-testing/FullSimTTBar-2_1_X_2008-07-08_STARTUP_V4-AODSIM.100.root'
#  'file:/afs/cern.ch/user/h/hegner/scratch0/PAT/testPatTuple_recHits_221.root'  
'/store/relval/CMSSW_2_2_3/RelValLM1_sfts/GEN-SIM-RECO/IDEAL_V11_v4/0003/4E82B5E1-BFCB-DD11-AEA2-001D09F24D8A.root'
#POOLSOURCE
    )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(MAXEVENTS) )

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = cms.string('STARTUP_V4::All')
process.GlobalTag.globaltag = 'IDEAL_V11::All'
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
#process.load("Configuration.StandardSequences.MagneticField_cff")

#--------------------------------------------------------------------
#needed for tcMET
#--------------------------------------------------------------------
# tracking geometry
process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
#det id
process.load("TrackingTools.TrackAssociator.DetIdAssociatorESProducer_cff")

# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patLayer0_cff")
process.load("PhysicsTools.PatAlgos.patLayer1_cff")
#process.content = cms.EDAnalyzer("EventContentAnalyzer")
#process.load("PhysicsTools.PatAlgos.patLayer1_cff")
process.load("Workspace.ConfigurableAnalysis.configurableAnalysis_cff")
process.load("SusyAnalysis.PatCrossCleaner.patCrossCleaner_cfi")

#-------------------------------------------------
#load met muon corrections
#-------------------------------------------------

process.load("JetMETCorrections.Type1MET.MetMuonCorrections_cff")

#-------------------------------------------------
# load tcMET producer
#-------------------------------------------------

process.load("RecoMET.METProducers.TCMET_cfi")

## Necessary fixes to run 2.2.X on 2.1.X data
from PhysicsTools.PatAlgos.tools.cmsswVersionTools import run22XonSummer08AODSIM
run22XonSummer08AODSIM(process)

# CaloTowerConstituentsMap needed for Electron/Photon-Jet cleaning
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.CaloTowerConstituentsMapBuilder = cms.ESProducer("CaloTowerConstituentsMapBuilder",
                                            MapFile = cms.untracked.string('Geometry/CaloTopology/data/CaloTowerEEGeometric.map.gz')
                                                         )

process.p = cms.Path(
                process.MetMuonCorrections*process.tcMet+
                process.patLayer0  
                #+ process.content # uncomment to get a dump of the output after layer 0
                + process.patLayer1
                + process.patcrosscleaner
                +
                process.configurableAnalysis
            )


# Output module configuration
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('PATLayer1_Output.fromAOD_full.root'),
    # save only events passing the full path
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    outputCommands = cms.untracked.vstring('keep *')
)
#process.outpath = cms.EndPath(process.out)

