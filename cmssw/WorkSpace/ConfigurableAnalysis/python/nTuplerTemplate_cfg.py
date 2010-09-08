# NTuple Configuration file to run on RECO files in 36X. Updated 27May10
import FWCore.ParameterSet.Config as cms
process = cms.Process("NTUP")

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10

## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

## PAT and cfA ntuple configurations.
process.load("PhysicsTools.PatAlgos.patSequences_cff")
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
from PhysicsTools.Configuration.SUSY_pattuple_cff import addDefaultSUSYPAT, getSUSY_pattuple_outputCommands
process.load("Workspace.ConfigurableAnalysis.configurableAnalysis_ForPattuple_cff")

## Output Module
process.out = cms.OutputModule("PoolOutputModule",
	fileName = cms.untracked.string('SUSYPAT.root'),
	# save only events passing the full path
	SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
	# save PAT Layer 1 output; you need a '*' to unpack the list of commands 'patEventContent'
	outputCommands = cms.untracked.vstring('drop *', *patEventContent ) 
)
process.out.splitLevel = cms.untracked.int32(99)  # Turn on split level (smaller files???)
process.out.overrideInputFileSplitLevels = cms.untracked.bool(True)
process.out.dropMetaData = cms.untracked.string('DROPPED')   # Get rid of metadata related to dropped collections
SUSY_pattuple_outputCommands = getSUSY_pattuple_outputCommands( process )
process.out.outputCommands = cms.untracked.vstring('drop *', *SUSY_pattuple_outputCommands)
#process.outpath = cms.EndPath(process.out) # if you want to save the PAT files.

################################## Most changes occur here ##############################
# Change input filenames and setup options for data/MC.
process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(
	'file:/DataE/wto/Reco/Commissioning10-SD_JetMETTau-Jun9thSkim_v1RECO.root')
)

#process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck') #check for dupes.
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#nTuple output file
process.TFileService.fileName= 'nTuple2.root' #nTuple filename.

#Apply SUSYPAT, parameters are: mcInfo, HLT menu, Jet energy corrections, JetCollections

#### Settings for Running on  MC 36X
#addDefaultSUSYPAT(process,True,'HLT','Spring10','',['IC5Calo','AK5JPT'])
#process.GlobalTag.globaltag = 'START36_V10A::All'
#process.GlobalTag.globaltag = 'MC_36Y_V9::All'
#process.GlobalTag.globaltag = 'START36_V9::All'

### Setting for Running on Commissioing10 361p4 Data, also change secVertex Btag.
#addDefaultSUSYPAT(process,False,'HLT','Spring10','',['IC5Calo','AK5JPT'])
#process.GlobalTag.globaltag = 'GR10_P_V6::All'

### Setting for Running on 362 MC with REDIGI36X HLT tag.
#addDefaultSUSYPAT(process,True,'REDIGI36X','Spring10','',['IC5Calo','AK5JPT'])
#process.GlobalTag.globaltag = 'START36_V9::All'

### Setting for Running on Commissioing10 358p3 Data, also change secVertex Btag.
#addDefaultSUSYPAT(process,False,'HLT','Spring10','35x',['IC5Calo','AK5JPT'])
#process.GlobalTag.globaltag = 'GR10_P_V5::All'

### Setting for Running on Jun9thReReco data 
addDefaultSUSYPAT(process,False,'HLT','Spring10','',['IC5Calo','AK5JPT']) 
#process.GlobalTag.globaltag = 'GR_R_37X_V6A::All'
#needed to run in 38X... seems to work on 37X RECO also...
process.GlobalTag.globaltag = cms.string('START38_V8::All') 

### Setting for Running on 37X MC 
#addDefaultSUSYPAT(process,True,'HLT','Spring10','',['IC5Calo','AK5JPT']) 
#process.GlobalTag.globaltag = 'START37_V5::All'

### Setting for Running on Jul16thReReco data 
#addDefaultSUSYPAT(process,False,'HLT','Spring10','',['IC5Calo','AK5JPT']) 
#process.GlobalTag.globaltag = 'GR_R_37X_V6D::All'

##########################################################################################

#-- Execution path ------------------------------------------------------------
process.p = cms.Path( process.susyPatDefaultSequence +process.configurableAnalysis)

#-- Dump config ------------------------------------------------------------
file = open('SusyPAT_cfg.py','w')
file.write(str(process.dumpPython()))
file.close()
