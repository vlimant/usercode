import FWCore.ParameterSet.Config as cms

IDconverttoBinNum = cms.PSet(
        ranges = cms.VPSet(cms.PSet(
            pdgIDs = cms.vint32(0),
            label = cms.string('ID=0')
        ), 
            cms.PSet(
                pdgIDs = cms.vint32(211, -211),
                label = cms.string('#pi+/-')
            ), 
            cms.PSet(
                pdgIDs = cms.vint32(321, -321, 130, -130),
                label = cms.string('K')
            ), 
            cms.PSet(
                pdgIDs = cms.vint32(411, -411, 421, -421, 431, 
                    -431),
                label = cms.string('D')
            ), 
            cms.PSet(
                pdgIDs = cms.vint32(521, -521, 511, -511, 531, 
                    -531),
                label = cms.string('B')
            ), 
            cms.PSet(
                pdgIDs = cms.vint32(5122, -5122),
                label = cms.string('#Lambda_{b}')
            ), 
            cms.PSet(
                pdgIDs = cms.vint32(443, -443),
                label = cms.string('J/#Psi')
            ), 
            cms.PSet(
                pdgIDs = cms.vint32(553, -553, 100553, -100553, 200553, 
                    -200553, 300553, -300553),
                label = cms.string('#Upsilon(nS)')
            ), 
            cms.PSet(
                pdgIDs = cms.vint32(24, -24),
                label = cms.string('W^{+/-}')
            ), 
            cms.PSet(
                pdgIDs = cms.vint32(23),
                label = cms.string('Z^{0}')
            ), 
            cms.PSet(
                pdgIDs = cms.vint32(15, -15),
                label = cms.string('#tau^{+/-}')
            ), 
            cms.PSet(
                pdgIDs = cms.vint32(13, -13),
                label = cms.string('#mu^{+/-}')
            )),
        title = cms.string('a set of PDGids')
    )
