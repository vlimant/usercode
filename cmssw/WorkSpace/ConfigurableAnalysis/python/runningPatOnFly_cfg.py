#
#  SUSY-PAT configuration file
#
#  PAT configuration for the SUSY group - 33X series
#  More information here:
#  https://twiki.cern.ch/twiki/bin/view/CMS/SusyPatLayer1DefV7
#

# Starting with a skeleton process which gets imported with the following line
#from PhysicsTools.PatAlgos.patTemplate_cfg import *


import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")

## Options and Output Report
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

## Source
process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
#            '/store/relval/CMSSW_3_3_0/RelValTTbar/GEN-SIM-RECO/STARTUP31X_V8-v1/0001/3291E09D-67B7-DE11-9ED6-003048678C9A.root'
'file:/tmp/rebassoo/F6E6E5CD-60D8-DE11-A93B-0019B9F705A3.root'
                )
                            )
## Maximal Number of Events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
## global tag for 33X
process.GlobalTag.globaltag = cms.string('STARTUP31X_V1::All')
## global tag for 34X
## process.GlobalTag.globaltag = cms.string('STARTUP3XY_V9::All')
process.load("Configuration.StandardSequences.MagneticField_cff")

## Standard PAT Configuration File
process.load("PhysicsTools.PatAlgos.patSequences_cff")

from PhysicsTools.PatAlgos.tools.jetTools import *
print "*********************************************************************"
print "Switching all processes to use the anti-kT algorithm by default."
print "Switch the jet collection to your desired algorithm if this is not"
print "what you want to use. Note that L7Parton correction are taken from"
print "SC5 instead of AK5. This is an intermediate solution for the time "
print "being."
print "*********************************************************************"
switchJetCollection(process,
                    cms.InputTag('ak5CaloJets'),
                    doJTA            = True,
                    doBTagging       = True,
                    jetCorrLabel     = ('AK5','Calo'),
                    doType1MET       = True,
                    genJetCollection = cms.InputTag("ak5GenJets")
                    )

# Output Module Configuration (expects a path 'p')
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('PATLayer1_Output.fromAOD_full.root'),
                               # save only events passing the full path
                               SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
                               # save PAT Layer 1 output; you need a '*' to
                               # unpack the list of commands 'patEventContent'
                               outputCommands = cms.untracked.vstring('drop *', *patEventContent )
                               )
#process.outpath = cms.EndPath(process.out)


#-- Meta data to be logged in DBS ---------------------------------------------
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.20 $'),
    name = cms.untracked.string('$Source: /cvs_server/repositories/CMSSW/CMSSW/PhysicsTools/Configuration/test/SUSY_pattuple_cfg.py,v $'),
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
     'rfio://?svcclass=cmscafuser&path=/castor/cern.ch/user/n/nmohr/QCDDiJet_Pt380to470_MC_31X_V9_ReReco332.root'
# 'rfio://?svcclass=cmscafuser&path=/castor/cern.ch/user/c/chiochia/09_beam_commissioning/BSCskim_123592_Express_bit40-41.root'
# 'rfio://?svcclass=cmscafuser&path=/castor/cern.ch/user/g/gpetrucc/900GeV/DATA/bit40-run123596-lumi130_143.root'
# 'rfio://?svcclass=cmscafuser&path=/castor/cern.ch/user/g/gpetrucc/900GeV/DATA/bit40-run123596-lumi_1_68.root'
#     'rfio://?svcclass=cmscafuser&path=/castor/cern.ch/user/g/gpetrucc/900GeV/DATA/bit40-run123596-lumi_68_129.root'
#         'rfio://?svcclass=cmscafuser&path=/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/123/596/FEE389F2-33E2-DE11-A62E-001617C3B76E.root'
#     'file:/tmp/rebassoo/F6E6E5CD-60D8-DE11-A93B-0019B9F705A3.root'
#    'rfio://?svcclass=cmscafuser&path=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/bit40or41skim.root'
#      'rfio://?svcclass=cmscafuser&path=/castor/cern.ch/user/c/chiochia/09_beam_commissioning/BSCskim_123151_Express.root'
#'rfio://?svcclass=cmscafuser&path=/castor/cern.ch/user/r/rebassoo/F6E6E5CD-60D8-DE11-A93B-0019B9F705A3.root',
#'rfio://?svcclass=cmscafuser&path=/castor/cern.ch/user/r/rebassoo/CECE0A1C-60D8-DE11-8B50-001D09F26C5C.root',
#'rfio://?svcclass=cmscafuser&path=/castor/cern.ch/user/r/rebassoo/C4F34BCE-60D8-DE11-BA15-0019B9F581C9.root'
#'rfio://?svcclass=cmscafuser&path=/castor/cern.ch/user/r/rebassoo/4259728C-61D8-DE11-BE08-000423D6A6F4.root',
#'rfio://?svcclass=cmscafuser&path=/castor/cern.ch/user/r/rebassoo/2E6A82EE-60D8-DE11-B655-001D09F2423B.root',
#'rfio://?svcclass=cmscafuser&path=/castor/cern.ch/user/r/rebassoo/245448AF-5ED8-DE11-9E67-001D09F2527B.root',
#'rfio://?svcclass=cmscafuser&path=/castor/cern.ch/user/r/rebassoo/2288321E-60D8-DE11-B57D-001D09F2AF1E.root',
#'rfio://?svcclass=cmscafuser&path=/castor/cern.ch/user/r/rebassoo/1653DF8D-61D8-DE11-A283-00304879FBB2.root',
#'rfio://?svcclass=cmscafuser&path=/castor/cern.ch/user/r/rebassoo/04E779D1-60D8-DE11-9811-001D09F24498.root'
#    'rfio://?svcclass=cmscafuser&path=/castor/cern.ch/user/r/rebassoo/C4F34BCE-60D8-DE11-BA15-0019B9F581C9.root'
#'/store/relval/CMSSW_3_3_3/RelValLM1_sfts/GEN-SIM-RECO/MC_31X_V9-v1/0003/465BA515-32D2-DE11-808A-000423D944FC.root'
    ]
process.maxEvents.input = -1
# Due to problem in production of LM samples: same event number appears multiple times
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')


#-- configurableAnalysis stuff -------------------------------------------------------
process.load("Workspace.ConfigurableAnalysis.configurableAnalysis_ForPattuple_cff")

#-- Calibration tag -----------------------------------------------------------
# Should match input file's tag
process.GlobalTag.globaltag = 'GR09_P_V6::All' #GR09_P_V6 for data, MC_31X_V9 for MC

#-- Missing ak5GenJets in 3.3.2 samples ---------------------------------------
# Comment this when you run on data
#from PhysicsTools.PatAlgos.tools.cmsswVersionTools import run33xOnReRecoMC
#run33xOnReRecoMC( process, "ak5GenJets" )

############################# START SUSYPAT specifics ####################################
from PhysicsTools.Configuration.SUSY_pattuple_cff import addDefaultSUSYPAT, removeMCDependence, getSUSY_pattuple_outputCommands
#addDefaultSUSYPAT(process,'HLT8E29') #second parameter is the name of the HLT menu
addDefaultSUSYPAT(process) #second parameter is the name of the HLT menu
# Uncomment next line when you run on data
removeMCDependence(process)
SUSY_pattuple_outputCommands = getSUSY_pattuple_outputCommands( process )
############################## END SUSYPAT specifics ####################################


#-- Output module configuration -----------------------------------------------
process.out.fileName = 'SUSYPAT.root'       # <-- CHANGE THIS TO SUIT YOUR NEEDS

# Custom settings
process.out.splitLevel = cms.untracked.int32(99)  # Turn on split level (smaller files???)
process.out.overrideInputFileSplitLevels = cms.untracked.bool(True)
process.out.dropMetaData = cms.untracked.string('DROPPED')   # Get rid of metadata related to dropped collections
process.out.outputCommands = cms.untracked.vstring('drop *', *SUSY_pattuple_outputCommands )

#-- Execution path ------------------------------------------------------------
# Full path
process.p = cms.Path( process.seqSUSYDefaultSequence+process.configurableAnalysis )

