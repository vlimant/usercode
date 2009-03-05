import FWCore.ParameterSet.Config as cms

hltMuonTreeMaker = cms.EDAnalyzer("IsoMuAnalyzer",
    trackExtractorPSet = cms.PSet(
        Diff_z = cms.double(0.2),
        inputTrackCollection = cms.InputTag("hltPixelTracks"),
        BeamSpotLabel = cms.InputTag("hltOfflineBeamSpot"),
        ComponentName = cms.string('TrackExtractor'),
        DR_Max = cms.double(0.24),
        Diff_r = cms.double(0.1),
        Chi2Prob_Min = cms.double(-1.0),
        DR_Veto = cms.double(0.01),
        NHits_Min = cms.uint32(0),
        Chi2Ndof_Max = cms.double(1e+64),
        Pt_Min = cms.double(-1.0),
        DepositLabel = cms.untracked.string('PXLS'),
        BeamlineOption = cms.string('BeamSpotFromEvent')
    ),
    l3MuonLabel = cms.InputTag("hltL3Muons"),
    singleMuIsoTriggerName = cms.string('HLT_IsoMu11'),
    diMuIsoTriggerName = cms.string('HLT_DoubleIsoMu3'),
    l2MuonLabel = cms.InputTag("hltL2Muons","UpdatedAtVtx"),
    propagatorName = cms.string('SteppingHelixPropagatorAlong'),
    l1MuonLabel = cms.InputTag("hltL1extraParticles"),
    trackCutsPSet = cms.PSet(
        ConeSizes = cms.vdouble(0.24, 0.24, 0.24, 0.24, 0.24, 
            0.24, 0.24, 0.24, 0.24, 0.24, 
            0.24, 0.24, 0.24, 0.24, 0.24, 
            0.24, 0.24, 0.24, 0.24, 0.24, 
            0.24, 0.24, 0.24, 0.24, 0.24, 
            0.24),
        EtaBounds = cms.vdouble(0.0435, 0.1305, 0.2175, 0.3045, 0.3915, 
            0.4785, 0.5655, 0.6525, 0.7395, 0.8265, 
            0.9135, 1.0005, 1.0875, 1.1745, 1.2615, 
            1.3485, 1.4355, 1.5225, 1.6095, 1.6965, 
            1.785, 1.88, 1.9865, 2.1075, 2.247, 
            2.411),
        L3IsoTrackCutsName = cms.string('SimpleCuts'),
        Thresholds = cms.vdouble(1.1, 1.1, 1.1, 1.1, 1.2, 
            1.1, 1.2, 1.1, 1.2, 1.0, 
            1.1, 1.0, 1.0, 1.1, 1.0, 
            1.0, 1.1, 0.9, 1.1, 0.9, 
            1.1, 1.0, 1.0, 0.9, 0.8, 
            0.1)
    ),
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
    ),
    # InputTag l2MuonLabel = hltL2Muons
    # InputTag muSysTransientRecHitsLabel = MuonTransientTrackingRecHitBuilderESProducer
    ServiceParameters = cms.PSet(
        Propagators = cms.untracked.vstring('SteppingHelixPropagatorAny', 
            'SteppingHelixPropagatorAlong', 
            'SteppingHelixPropagatorOpposite', 
            'PropagatorWithMaterial', 
            'PropagatorWithMaterialOpposite', 
            'SmartPropagator', 
            'SmartPropagatorOpposite', 
            'SmartPropagatorAnyOpposite', 
            'SmartPropagatorAny', 
            'SmartPropagatorRK', 
            'SmartPropagatorAnyRK'),
        RPCLayers = cms.bool(True),
        UseMuonNavigation = cms.untracked.bool(True)
    ),
    triggerResults_ = cms.InputTag("TriggerResults","","HLT"),
    linkLabel = cms.InputTag("hltL3Muons"),
    calExtractorPSet = cms.PSet(
        DR_Veto_H = cms.double(0.1),
        Vertex_Constraint_Z = cms.bool(False),
        Threshold_H = cms.double(0.5),
        ConeSizes = cms.vdouble(0.24, 0.24, 0.24, 0.24, 0.24, 
            0.24, 0.24, 0.24, 0.24, 0.24, 
            0.24, 0.24, 0.24, 0.24, 0.24, 
            0.24, 0.24, 0.24, 0.24, 0.24, 
            0.24, 0.24, 0.24, 0.24, 0.24, 
            0.24),
        Thresholds = cms.vdouble(4.0, 3.7, 4.0, 3.5, 3.4, 
            3.4, 3.2, 3.4, 3.1, 2.9, 
            2.9, 2.7, 3.1, 3.0, 2.4, 
            2.1, 2.0, 2.3, 2.2, 2.4, 
            2.5, 2.5, 2.6, 2.9, 3.1, 
            2.9),
        DR_Max = cms.double(0.24),
        DR_Veto_E = cms.double(0.07),
        Vertex_Constraint_XY = cms.bool(False),
        Threshold_E = cms.double(0.2),
        Weight_E = cms.double(1.5),
        EtaBounds = cms.vdouble(0.0435, 0.1305, 0.2175, 0.3045, 0.3915, 
            0.4785, 0.5655, 0.6525, 0.7395, 0.8265, 
            0.9135, 1.0005, 1.0875, 1.1745, 1.2615, 
            1.3485, 1.4355, 1.5225, 1.6095, 1.6965, 
            1.785, 1.88, 1.9865, 2.1075, 2.247, 
            2.411),
        DepositLabel = cms.untracked.string('EcalPlusHcal'),
        calName = cms.string('CaloExtractor'),
        CaloTowerCollectionLabel = cms.InputTag("hltTowerMakerForMuons"),
        Weight_H = cms.double(1.0)
    ),
    singleMuNonIsoTriggerName = cms.string('HLT_Mu15'),
    diMuNonIsoTriggerName = cms.string('HLT_DoubleMu3'),
    l2MuonSeeds = cms.InputTag("hltL2MuonSeeds"),
    l3AssociatorName = cms.string('AssociatorByDeltaR1.0'),
    matrixPset = cms.PSet(
        action = cms.string('use'),
        atIP = cms.bool(True),
        errorMatrixValuesPSet = cms.PSet(
            pf3_V12 = cms.PSet(
                action = cms.string('scale'),            
                values = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0, 
                    1.0, 1.0, 1.0, 1.0, 1.0, 
                    1.0, 1.0)
            ),
            pf3_V13 = cms.PSet(
                action = cms.string('scale'),
                values = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0, 
                    1.0, 1.0, 1.0, 1.0, 1.0, 
                    1.0, 1.0)
            ),
            pf3_V11 = cms.PSet(
                action = cms.string('scale'),
                values = cms.vdouble(3.0, 3.0, 3.0, 5.0, 4.0, 
                    5.0, 10.0, 7.0, 10.0, 10.0, 
                    10.0, 10.0)
            ),
            pf3_V25 = cms.PSet(
                action = cms.string('scale'),
                values = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0, 
                    1.0, 1.0, 1.0, 1.0, 1.0, 
                    1.0, 1.0)
            ),
            pf3_V14 = cms.PSet(
                action = cms.string('scale'),
                values = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0, 
                    1.0, 1.0, 1.0, 1.0, 1.0, 
                    1.0, 1.0)
            ),
            pf3_V15 = cms.PSet(
                action = cms.string('scale'),
                values = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0, 
                    1.0, 1.0, 1.0, 1.0, 1.0, 
                    1.0, 1.0)
            ),
            pf3_V34 = cms.PSet(
                action = cms.string('scale'),
                values = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0, 
                    1.0, 1.0, 1.0, 1.0, 1.0, 
                    1.0, 1.0)
            ),
            yAxis = cms.vdouble(0.0, 1.0, 1.4, 2.5),
            pf3_V45 = cms.PSet(
                action = cms.string('scale'),
                values = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0, 
                    1.0, 1.0, 1.0, 1.0, 1.0, 
                    1.0, 1.0)
            ),
            pf3_V44 = cms.PSet(
                action = cms.string('scale'),
                values = cms.vdouble(3.0, 3.0, 3.0, 5.0, 4.0, 
                    5.0, 10.0, 7.0, 10.0, 10.0, 
                    10.0, 10.0)
            ),
            xAxis = cms.vdouble(0.0, 13.0, 30.0, 70.0, 500.0),
            pf3_V23 = cms.PSet(
                action = cms.string('scale'),
                values = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0, 
                    1.0, 1.0, 1.0, 1.0, 1.0, 
                    1.0, 1.0)
            ),
            pf3_V22 = cms.PSet(
                action = cms.string('scale'),
                values = cms.vdouble(3.0, 3.0, 3.0, 5.0, 4.0, 
                    5.0, 10.0, 7.0, 10.0, 10.0, 
                    10.0, 10.0)
            ),
            pf3_V55 = cms.PSet(
                action = cms.string('scale'),
                values = cms.vdouble(3.0, 3.0, 3.0, 5.0, 4.0, 
                    5.0, 10.0, 7.0, 10.0, 10.0, 
                    10.0, 10.0)
            ),
            zAxis = cms.vdouble(-3.14159, 3.14159),
            pf3_V35 = cms.PSet(
                action = cms.string('scale'),
                values = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0, 
                    1.0, 1.0, 1.0, 1.0, 1.0, 
                    1.0, 1.0)
            ),
            pf3_V33 = cms.PSet(
                action = cms.string('scale'),
                values = cms.vdouble(3.0, 3.0, 3.0, 5.0, 4.0, 
                    5.0, 10.0, 7.0, 10.0, 10.0, 
                    10.0, 10.0)
            ),
            pf3_V24 = cms.PSet(
                action = cms.string('scale'),
                values = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0, 
                    1.0, 1.0, 1.0, 1.0, 1.0, 
                    1.0, 1.0)
            )
        )
    ),
    trackingParticleLabel = cms.InputTag("tpMuon"),
    l2AssociatorName = cms.string('AssociatorByDeltaR1.0'),
    beamSpotLabel = cms.InputTag("hltOfflineBeamSpot"),


    TimerLabel = cms.InputTag("hltTimer"),
    TimingModules = cms.untracked.PSet(
    MuonL3IsoModules = cms.untracked.vstring('hltPixelTracks',
                                             'hltL3MuonIsolations'),
    TrackerRecModules = cms.untracked.vstring('hltSiPixelClusters',
                                              'hltSiPixelRecHits',
                                              'hltSiStripClusters',
                                              'hltSiStripRawToClustersFacility'),
    MuonLocalRecModules = cms.untracked.vstring('hltDt1DRecHits',
                                                'hltDt4DSegments',
                                                'hltRpcRecHits',
                                                'hltCsc2DRecHits',
                                                'hltCscSegments'),
    CaloDigiModules = cms.untracked.vstring('hltEcalPreshowerDigis',
                                            'hltEcalRegionalMuonsFEDs',
                                            'hltEcalRegionalMuonsDigis',
                                            'hltHcalDigis'),
    MuonL3RecModules = cms.untracked.vstring('hltL3TrajectorySeed',
                                             'hltL3TrackCandidateFromL2',
                                             'hltL3Muons',
                                             'hltL3MuonCandidates'),
    CaloRecModules = cms.untracked.vstring('hltEcalRegionalMuonsWeightUncalibRecHit',
                                           'hltEcalRegionalMuonsRecHitTmp',
                                           'hltEcalRegionalMuonsRecHit',
                                           'hltEcalPreshowerRecHit',
                                           'hltHbhereco',
                                           'hltHoreco',
                                           'hltHfreco',
                                           'hltTowerMakerForMuons'),
    MuonL2IsoModules = cms.untracked.vstring('hltL2MuonIsolations'),
    MuonDigiModules = cms.untracked.vstring('hltMuonCSCDigis',
                                            'hltMuonDTDigis',
                                            'hltMuonRPCDigis'),
    TrackerDigiModules = cms.untracked.vstring('hltSiPixelDigis'),
    MuonL2RecModules = cms.untracked.vstring('hltL2MuonSeeds',
                                             'hltL2Muons',
                                             'hltL2MuonCandidates'),
    )
)
