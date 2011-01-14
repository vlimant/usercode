#
#  SUSY-PAT configuration file
#
#  PAT configuration for the SUSY group - 35X/36X series
#  More information here:
#  https://twiki.cern.ch/twiki/bin/view/CMS/SusyPatLayer1DefV8
#

# Starting with a skeleton process which gets imported with the following line
#from PhysicsTools.PatAlgos.patTemplate_cfg import *

import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")

## Options and Output Report
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

## Source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      #'/store/relval/CMSSW_3_6_0/RelValTTbar/GEN-SIM-RECO/START36_V4-v1/0013/306F945C-9A49-DF11-85F8-0018F3D0965A.root'
    )
)
## Maximal Number of Events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('START36_V4::All')
process.load("Configuration.StandardSequences.MagneticField_cff")

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
                               #fileName = cms.untracked.string('PATLayer1_Output.fromAOD_full.root'),
                               # save only events passing the full path
                               SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
                               # save PAT Layer 1 output; you need a '*' to
                               # unpack the list of commands 'patEventContent'
                               outputCommands = cms.untracked.vstring('drop *', "keep *_BFieldColl_*_*_",*patEventContent ) 
                               )
#process.outpath = cms.EndPath(process.out)


#-- Meta data to be logged in DBS ---------------------------------------------
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.12 $'),
    name = cms.untracked.string('$Source: /cvs_server/repositories/CMSSW/UserCode/JRVlimant/cmssw/WorkSpace/ConfigurableAnalysis/python/runningPatOnFly_cfg.py,v $'),
    annotation = cms.untracked.string('SUSY pattuple definition')
)

#-- Message Logger ------------------------------------------------------------
process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.PATSummaryTables = cms.untracked.PSet(
    limit = cms.untracked.int32(-1),
    reportEvery = cms.untracked.int32(1)
    )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

#-- Input Source --------------------------------------------------------------
process.source.fileNames = [
     #'/store/relval/CMSSW_3_8_2/RelValZmumuJets_Pt_20_300_GEN/GEN-SIM-RECO/MC_38Y_V9_PU_E7TeV_AVE_2_BX2808-v1/0019/EE6557E2-07B1-DF11-B5A9-0026189438D3.root'
     #'/store/relval/CMSSW_3_8_3/RelValTTbar/GEN-SIM-RECO/START38_V9-v1/0022/CA9763E0-EFBF-DF11-81C5-002618943845.root'
     #'/store/data/Run2010B/Jet/RECO/PromptReco-v2/000/146/331/16DFEFAD-DEC5-DF11-9E29-0030487CD6DA.root'
     #'file:/DataE/wto/Reco/Mu_Run2010B-PromptReco-v2_RECO/D47EAE08-ADC6-DF11-B1D8-0030487CD6D8.root'
     #'/store/data/Run2010B/Electron/RECO/PromptReco-v2/000/146/511/52C66503-9EC7-DF11-B04C-001D09F2A465.root'
     #'file:/DataF/pbgeff/temp_385_ntuple/LM0_SUSY_sftsht_7TeV-pythia6_Fall10-START38_V12-v1.root'
     'file:/LVM/SATA/pbgeff/temp_387_ntuple/Jet_Nov4ReReco_06C7ADF6-60EC-DF11-A937-003048D436C6.root'
    ]
process.maxEvents.input = 20
# Due to problem in production of LM samples: same event number appears multiple times
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

#-- Calibration tag -----------------------------------------------------------
# Should match input file's tag
#process.GlobalTag.globaltag = 'START38_V14::All' # for MC
process.GlobalTag.globaltag = 'GR_R_38X_V15::All' #for Nov4 rereco

############################# START SUSYPAT specifics ####################################
from PhysicsTools.Configuration.SUSY_pattuple_cff import addDefaultSUSYPAT, getSUSY_pattuple_outputCommands
#Apply SUSYPAT, parameters are: mcInfo, HLT menu, Jet energy corrections, mcVersion ('35x' for 35x samples, empty string for 36X samples),JetCollections
addDefaultSUSYPAT(process,False,'HLT',['L2Relative','L3Absolute','L2L3Residual'],'',['AK5PF','AK5JPT'])
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
process.out.fileName = 'SUSYPAT.root'       # <-- CHANGE THIS TO SUIT YOUR NEEDS

# Custom settings
process.out.splitLevel = cms.untracked.int32(99)  # Turn on split level (smaller files???)
process.out.overrideInputFileSplitLevels = cms.untracked.bool(True)
process.out.dropMetaData = cms.untracked.string('DROPPED')   # Get rid of metadata related to dropped collections
#process.out.outputCommands = cms.untracked.vstring('drop *', *SUSY_pattuple_outputCommands )
process.out.outputCommands = cms.untracked.vstring('drop *',"keep *_HBHENoiseFilterResultProducer_*_*","keep *_BFieldColl_*_*", *SUSY_pattuple_outputCommands )

#-- Execution path ------------------------------------------------------------
# Full path
#This is to run on full sim or data
process.p = cms.Path(process.HBHENoiseFilterResultProducer + process.BFieldColl + process.susyPatDefaultSequence + process.configurableAnalysis)
#This is to run on FastSim
#process.p = cms.Path( process.BFieldColl + process.susyPatDefaultSequence +process.configurableAnalysis)


#-- Dump config ------------------------------------------------------------
file = open('SusyPAT_cfg.py','w')
file.write(str(process.dumpPython()))
file.close()
