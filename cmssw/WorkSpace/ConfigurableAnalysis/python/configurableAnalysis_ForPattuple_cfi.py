import FWCore.ParameterSet.Config as cms

#Main module to create the two trees and histograms (these histograms are not from the trees)
configurableAnalysis = cms.EDFilter("ConfigurableAnalysis",
    Selections = cms.PSet(
        filters = cms.PSet(
        ),
        selections = cms.PSet(
            minSelection = cms.PSet(
                filterOrder = cms.vstring(''),
                makeFinalPlots = cms.bool(True),
                makeSummaryTable = cms.bool(True),
                makeContentPlots = cms.bool(True),
                makeAllButOnePlots = cms.bool(True),
                ntuplize = cms.bool(True),
                nMonitor = cms.uint32(1000),
                makeCumulativePlots = cms.bool(True)
            )
        )
    ),
    Plotter = cms.PSet(
        TH1s = cms.PSet(

        ),
        ComponentName = cms.string('VariablePlotter'),
        TProfiles = cms.PSet(

        ),
        TH2s = cms.PSet(

        )
    ),
    Variables = cms.PSet(

    ),
    workAsASelector = cms.bool(True),
    flows = cms.vstring('minSelection'),
    Ntupler = cms.PSet(
        branchesPSet = cms.PSet(
    
            treeName = cms.string('eventB'),

									


			),
			ComponentName = cms.string('CompleteNTupler'),
			useTFileService = cms.bool(True) #, ## false for EDM; true for non EDM

	)
)




