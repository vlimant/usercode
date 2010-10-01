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

#Need this for L1 triggers with CMSSW >= 381
process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff")
process.patTrigger.addL1Algos = cms.bool( True )

## Output Module Configuration (expects a path 'p')
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out = cms.OutputModule("PoolOutputModule",
     			       fileName = cms.untracked.string('patTuple.root'),
                               #fileName = cms.untracked.string('PATLayer1_Output.fromAOD_full.root'),
                               # save only events passing the full path
                               SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
                               # save PAT Layer 1 output; you need a '*' to
                               # unpack the list of commands 'patEventContent'
                               outputCommands = cms.untracked.vstring('drop *', *patEventContent ) 
                               )
#process.outpath = cms.EndPath(process.out)



#-- Meta data to be logged in DBS ---------------------------------------------
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.6 $'),
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
     '/store/relval/CMSSW_3_8_3/RelValTTbar/GEN-SIM-RECO/START38_V9-v1/0022/CA9763E0-EFBF-DF11-81C5-002618943845.root'
     #'/store/data/Run2010B/Jet/RECO/PromptReco-v2/000/146/331/16DFEFAD-DEC5-DF11-9E29-0030487CD6DA.root'
     #'file:/DataE/wto/Reco/Mu_Run2010B-PromptReco-v2_RECO/D47EAE08-ADC6-DF11-B1D8-0030487CD6D8.root'
    ]
process.maxEvents.input = 200
# Due to problem in production of LM samples: same event number appears multiple times
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

#-- Calibration tag -----------------------------------------------------------
# Should match input file's tag
#process.GlobalTag.globaltag = 'START38_V9::All'
#process.GlobalTag.globaltag = 'START38_V10::All'
#process.GlobalTag.globaltag = 'GR_R_38X_V9::All'
#process.GlobalTag.globaltag = 'GR_R_38X_V11::All'
#process.GlobalTag.globaltag = 'GR10_P_V9::All'
process.GlobalTag.globaltag = 'GR10_P_V10::All'
#process.GlobalTag.globaltag = 'MC_38Y_V9::All'

############################# START SUSYPAT specifics ####################################
from PhysicsTools.Configuration.SUSY_pattuple_cff import addDefaultSUSYPAT, getSUSY_pattuple_outputCommands
#Apply SUSYPAT, parameters are: mcInfo, HLT menu, Jet energy corrections, mcVersion ('35x' for 35x samples, empty string for 36X samples),JetCollections
addDefaultSUSYPAT(process,False,'HLT','Spring10','',['IC5Calo','AK5JPT']) 
#addDefaultSUSYPAT(process,True,'HLT','Spring10','',['IC5Calo','AK5JPT'])
SUSY_pattuple_outputCommands = getSUSY_pattuple_outputCommands( process )
############################## END SUSYPAT specifics ####################################


#-- configurableAnalysis stuff -------------------------------------------------------
process.load("Workspace.ConfigurableAnalysis.configurableAnalysis_ForPattuple_cff")


#-- Output module configuration -----------------------------------------------
process.out.fileName = 'SUSYPAT.root'       # <-- CHANGE THIS TO SUIT YOUR NEEDS

# Custom settings
process.out.splitLevel = cms.untracked.int32(99)  # Turn on split level (smaller files???)
process.out.overrideInputFileSplitLevels = cms.untracked.bool(True)
process.out.dropMetaData = cms.untracked.string('DROPPED')   # Get rid of metadata related to dropped collections
process.out.outputCommands = cms.untracked.vstring('drop *', *SUSY_pattuple_outputCommands )

#-- Execution path ------------------------------------------------------------
# Full path
process.p = cms.Path( process.susyPatDefaultSequence +process.configurableAnalysis)
#-- Dump config ------------------------------------------------------------
file = open('SusyPAT_cfg.py','w')
file.write(str(process.dumpPython()))
file.close()
