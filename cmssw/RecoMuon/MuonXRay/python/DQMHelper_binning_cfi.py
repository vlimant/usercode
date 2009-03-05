import FWCore.ParameterSet.Config as cms

pTbinning = cms.PSet(
    nBins = cms.uint32(400),
    Max = cms.double(200.0),
    Min = cms.double(0.0)
)
d0binning = cms.PSet(
    nBins = cms.uint32(200),
    Max = cms.double(0.1),
    Min = cms.double(-0.1)
)
etabinning = cms.PSet(
    nBins = cms.uint32(120),
    Max = cms.double(3.0),
    Min = cms.double(-3.0)
)
Aetabinning = cms.PSet(
    nBins = cms.uint32(60),
    Max = cms.double(3.0),
    Min = cms.double(0.0)
)
phibinning = cms.PSet(
    nBins = cms.uint32(400),
    Max = cms.double(3.14159265359),
    Min = cms.double(-3.14159265359)
)

