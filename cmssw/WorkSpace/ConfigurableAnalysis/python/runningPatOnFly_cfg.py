#
#  SUSY-PAT configuration file
#
#  PAT configuration for the SUSY group - 42X series
#  More information here:
#  https://twiki.cern.ch/twiki/bin/view/CMS/SusyPatLayer1DefV10
#

# Starting with a skeleton process which gets imported with the following line
#from PhysicsTools.PatAlgos.patTemplate_cfg import *

import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")

#-- Message Logger ------------------------------------------------------------
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.PATSummaryTables = cms.untracked.PSet(
    limit = cms.untracked.int32(-1),
    reportEvery = cms.untracked.int32(100)
)

#-- Source information ------------------------------------------------------
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
			'file:/LVM/SATA/wto/RECO/Run166512/MuHad/FA3D4FB8-BA91-E011-97AB-003048673374.root',
    )
		#maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

## Standard PAT Configuration File
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.BFieldColl = cms.EDProducer('BFieldProducer')
process.JetCorrColl = cms.EDProducer('JetCorrProducer')

#Need this for L1 triggers with CMSSW >= 381
process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff")
process.patTrigger.addL1Algos = cms.bool( True )

## Output Module Configuration (expects a path 'p')
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out = cms.OutputModule("PoolOutputModule",
     verbose = cms.untracked.bool(True),
     fileName = cms.untracked.string('patTuple.root'),
     SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
     outputCommands = cms.untracked.vstring('drop *', "keep *_BFieldColl_*_*_","keep *_JetCorrColl_*_*_",*patEventContent )
)
#process.outpath = cms.EndPath(process.out)

#-- SUSYPAT and GlobalTag Settings -----------------------------------------------------------
#process.GlobalTag.globaltag = 'START42_V13::All' # MC Setting
from PhysicsTools.Configuration.SUSY_pattuple_cff import addDefaultSUSYPAT, getSUSY_pattuple_outputCommands

process.GlobalTag.globaltag = 'GR_R_42_V19::All' 	# Data Setting
addDefaultSUSYPAT(process,False,'HLT',['L1FastJet','L2Relative','L3Absolute','L2L3Residual'],'',['AK5PF','AK5JPT'])
process.metJESCorAK5PFTypeI.corrector = cms.string('ak5PFL2L3Residual') # Type1PFMET Residual for data only.
#addDefaultSUSYPAT(process,True,'HLT',['L1FastJet','L2Relative','L3Absolute'],'',['AK5PF','AK5JPT'])
SUSY_pattuple_outputCommands = getSUSY_pattuple_outputCommands( process )
############################## END SUSYPAT specifics ####################################

#Turn on trigger info
from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger
switchOnTrigger(process, triggerProducer='patTrigger', triggerEventProducer='patTriggerEvent', sequence='patDefaultSequence', hltProcess="HLT")

process.load("Workspace.ConfigurableAnalysis.configurableAnalysis_ForPattuple_cff")
process.load('CommonTools/RecoAlgos/HBHENoiseFilterResultProducer_cfi')

#-- Output module configuration -----------------------------------------------
process.out.fileName = "SUSYPAT.root" 
process.out.splitLevel = cms.untracked.int32(99)  # Turn on split level (smaller files???)
process.out.overrideInputFileSplitLevels = cms.untracked.bool(True)
process.out.dropMetaData = cms.untracked.string('DROPPED')   # Get rid of metadata related to dropped collections
process.out.outputCommands = cms.untracked.vstring('drop *',"keep *_HBHENoiseFilterResultProducer_*_*","keep *_BFieldColl_*_*","keep *_JetCorrectionColl_*_*", *SUSY_pattuple_outputCommands )

#-- Execution path ------------------------------------------------------------
# Full path
#This is to run on full sim or data
process.p = cms.Path(process.HBHENoiseFilterResultProducer + process.BFieldColl + process.susyPatDefaultSequence + process.JetCorrColl + process.configurableAnalysis)
#This is to run on FastSim
#process.p = cms.Path( process.BFieldColl + process.susyPatDefaultSequence + process.JetCorrColl +process.configurableAnalysis)

#-- Dump config ------------------------------------------------------------
file = open('SusyPAT_cfg.py','w')
file.write(str(process.dumpPython()))
file.close()
