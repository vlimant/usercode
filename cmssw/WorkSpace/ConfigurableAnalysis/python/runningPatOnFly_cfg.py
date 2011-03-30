# Creates cfA nTuple running on Reco->PAT->cfA
import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.MessageLogger.suppressWarning = cms.untracked.vstring(
	'patTriggerPF','patTrigger','patTriggerPF'
)

## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

## Source
process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(
		'file:/home/wto/cmssw/output_1_1_3IE.root'
	)
)
process.GlobalTag.globaltag = "GR_R_39X_V5::All" #for Dec22ReReco
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

## Standard PAT Configuration File
process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.BFieldColl = cms.EDProducer('BFieldProducer')

##Need this for L1 triggers with CMSSW >= 381
process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff")
process.patTrigger.addL1Algos = cms.bool( True )

#For electron ID
process.patElectrons.addElectronID = cms.bool(True)
process.patElectrons.electronIDSources = cms.PSet(
    simpleEleId95relIso= cms.InputTag("simpleEleId95relIso"),
    simpleEleId90relIso= cms.InputTag("simpleEleId90relIso"),
    simpleEleId85relIso= cms.InputTag("simpleEleId85relIso"),
    simpleEleId80relIso= cms.InputTag("simpleEleId80relIso"),
    simpleEleId70relIso= cms.InputTag("simpleEleId70relIso"),
    simpleEleId60relIso= cms.InputTag("simpleEleId60relIso"),
    simpleEleId95cIso= cms.InputTag("simpleEleId95cIso"),
    simpleEleId90cIso= cms.InputTag("simpleEleId90cIso"),
    simpleEleId85cIso= cms.InputTag("simpleEleId85cIso"),
    simpleEleId80cIso= cms.InputTag("simpleEleId80cIso"),
    simpleEleId70cIso= cms.InputTag("simpleEleId70cIso"),
    simpleEleId60cIso= cms.InputTag("simpleEleId60cIso"),
    eidLoose = cms.InputTag("eidLoose"),
    eidTight = cms.InputTag("eidTight"),
    eidRobustLoose = cms.InputTag("eidRobustLoose"),
    eidRobustTight = cms.InputTag("eidRobustTight"),
    eidRobustHighEnergy = cms.InputTag("eidRobustHighEnergy"),
)
process.load("ElectroWeakAnalysis.WENu.simpleEleIdSequence_cff")
process.patElectronIDs = cms.Sequence(process.simpleEleIdSequence)
process.makePatElectrons = cms.Sequence(process.patElectronIDs*process.electronMatch*process.patElectrons)


## Output Module Configuration (expects a path 'p')
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out = cms.OutputModule("PoolOutputModule",
	verbose = cms.untracked.bool(True),
	fileName = cms.untracked.string('patTuple.root'),
	SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
	outputCommands = cms.untracked.vstring('drop *', "keep *_BFieldColl_*_*_",*patEventContent ) 
)
#unComment to save PAT file.
#process.outpath = cms.EndPath(process.out)

#-- Message Logger ------------------------------------------------------------
process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.PATSummaryTables = cms.untracked.PSet(
    limit = cms.untracked.int32(-1),
    reportEvery = cms.untracked.int32(1)
    )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

############################# START SUSYPAT specifics ####################################
from PhysicsTools.Configuration.SUSY_pattuple_cff import addDefaultSUSYPAT, getSUSY_pattuple_outputCommands
#Apply SUSYPAT, parameters are: mcInfo, HLT menu, 
				#Jet energy corrections, 
				#mcVersion,JetCollections
addDefaultSUSYPAT(process,False,'HLT',
									['L2Relative','L3Absolute','L2L3Residual'],
									'',['AK5PF','AK5JPT'])
#addDefaultSUSYPAT(process,True,'HLT',['L2Relative','L3Absolute'],'',['AK5PF','AK5JPT'])
#addDefaultSUSYPAT(process,True,'REDIGI38X',['L2Relative','L3Absolute'],'',['AK5PF','AK5JPT']) 
SUSY_pattuple_outputCommands = getSUSY_pattuple_outputCommands( process )
############################## END SUSYPAT specifics ####################################


#-- configurableAnalysis stuff -------------------------------------------------------
process.load("Workspace.ConfigurableAnalysis.configurableAnalysis_ForPattuple_cff")

process.load('CommonTools/RecoAlgos/HBHENoiseFilterResultProducer_cfi')

#Only run this for data
process.metJESCorAK5PFTypeI.corrector = cms.string('ak5PFL2L3Residual')

#-- Output module configuration -----------------------------------------------
process.out.fileName = 'SUSYPAT.root'       

# Custom settings
process.out.splitLevel = cms.untracked.int32(99)  # Turn on split level (smaller files???)
process.out.overrideInputFileSplitLevels = cms.untracked.bool(True)
process.out.outputCommands = cms.untracked.vstring('drop *',"keep *_HBHENoiseFilterResultProducer_*_*","keep *_BFieldColl_*_*", *SUSY_pattuple_outputCommands )

#-- Execution path ------------------------------------------------------------
process.p = cms.Path(	process.HBHENoiseFilterResultProducer 
										+ process.BFieldColl 
										+ process.susyPatDefaultSequence 
										+ process.configurableAnalysis)


#-- Dump config ------------------------------------------------------------
file = open('SusyPAT_cfg.py','w')
file.write(str(process.dumpPython()))
file.close()
