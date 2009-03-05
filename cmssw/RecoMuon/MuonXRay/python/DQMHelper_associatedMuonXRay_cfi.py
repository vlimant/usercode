import FWCore.ParameterSet.Config as cms
from RecoMuon.MuonXRay.DQMHelper_binning_cfi import *

DQMHelper_associatedMuonXRay = cms.PSet(
        H1s = cms.PSet(
            h_mu_eta = cms.PSet(
                xAxis = cms.PSet(
                    etabinning,
                    Label = cms.string('muon #eta^{reco}')
                ),
                title = cms.string('muon #eta^{reco}')
            ),
            h_mu_d0 = cms.PSet(
                xAxis = cms.PSet(
                    d0binning,
                    Label = cms.string('muon d_{0}^{reco} [cm]')
                ),
                title = cms.string('muon d_{0}^{reco}')
            ),
            h_mu_relDpT = cms.PSet(
                xAxis = cms.PSet(
                    nBins = cms.uint32(400),
                    Min = cms.double(-10.0),
                    Max = cms.double(10.0),
                    Label = cms.string('2*#Delta(L3Tk p_{T},L2 p_{T})/#Sigma(L3Tk p_{T},L2 p_{T})')
                ),
                title = cms.string('muon relative kink')
            ),
            num_reco = cms.PSet(
                nBins = cms.uint32(21),
                title = cms.string('number of reconstructed muons.'),
                Min = cms.double(-0.5),
                Max = cms.double(20.5),
                Label = cms.string('Nb. of reco #mu')
            ),
            h_mu_90pt = cms.PSet(
                xAxis = cms.PSet(
                    pTbinning,
                    Label = cms.string('muon 90% p_{T}^{reco} [GeV]')
                ),
                title = cms.string('muon 90% p_{T}^{reco}')
            ),
            h_mu_hits = cms.PSet(
                nBins = cms.uint32(81),
                title = cms.string('number of RecHit'),
                Min = cms.double(-0.5),
                Max = cms.double(80.5),
                Label = cms.string('Nb. of Hits')
            ),
            h_mu_sim_Aeta = cms.PSet(
                xAxis = cms.PSet(
                    Aetabinning,
                    Label = cms.string('muon |#eta^{sim}|')
                ),
                title = cms.string('muon |#eta^{sim}|')
            ),
            h_mu_pt = cms.PSet(
                xAxis = cms.PSet(
                    pTbinning,
                    Label = cms.string('muon p_{T}^{reco} [GeV]')
                ),
                title = cms.string('muon p_{T}^{reco}')
            ),
            h_mu_sim_phi = cms.PSet(
                xAxis = cms.PSet(
                    phibinning,
                    Label = cms.string('muon #varphi^{sim}')
                ),
                title = cms.string('muon #varphi^{sim}')
            ),
            h_mu_leading_pt = cms.PSet(
                xAxis = cms.PSet(
                    pTbinning,
                    Label = cms.string('leading muon p_{T}^{reco} [GeV]')
                ),
                title = cms.string('leading muon p_{T}^{reco}')
            ),
            h_mu_sim_pt = cms.PSet(
                xAxis = cms.PSet(
                    pTbinning,
                    Label = cms.string('muon p_{T}^{sim} [GeV]')
                ),
                title = cms.string('muon p_{T}^{sim}')
            ),
            h_mu_DpT = cms.PSet(
                xAxis = cms.PSet(
                    nBins = cms.uint32(400),
                    Min = cms.double(-100.0),
                    Max = cms.double(100.0),
                    Label = cms.string('#Delta(L3Tk p_{T},L2 p_{T}) [GeV]')
                ),
                title = cms.string('muon kink')
            ),
            h_mu_phi = cms.PSet(
                xAxis = cms.PSet(
                    phibinning,
                    Label = cms.string('muon #varphi^{reco}')
                ),
                title = cms.string('muon #varphi^{reco}')
            ),
            h_mu_Aeta = cms.PSet(
                xAxis = cms.PSet(
                    Aetabinning,
                    Label = cms.string('muon |#eta^{reco}|')
                ),
                title = cms.string('muon |#eta^{reco}|')
            ),
            num_reco_confirm = cms.PSet(
                nBins = cms.uint32(21),
                title = cms.string('number of reconstructed muons confirmed by trigger'),
                Min = cms.double(-0.5),
                Max = cms.double(20.5),
                Label = cms.string('Nb. of reco #mu making the trig.')
            ),
            h_mu_leading_phi = cms.PSet(
                xAxis = cms.PSet(
                    phibinning,
                    Label = cms.string('leading muon #varphi^{reco} [GeV]')
                ),
                title = cms.string('leading muon #varphi^{reco}')
            ),
            h_mu_leading_Aeta = cms.PSet(
                xAxis = cms.PSet(
                    Aetabinning,
                    Label = cms.string('leading muon |#eta^{reco}|')
                ),
                title = cms.string('leading muon |#eta^{reco}|')
            ),
            h_mu_leading_sim_pt = cms.PSet(
                xAxis = cms.PSet(
                    pTbinning,
                    Label = cms.string('leading muon p_{T}^{sim} [GeV]')
                ),
                title = cms.string('leading muon p_{T}^{sim}')
            ),
            h_mu_sim_eta = cms.PSet(
                xAxis = cms.PSet(
                    etabinning,
                    Label = cms.string('muon #eta^{sim}')
                ),
                title = cms.string('muon #eta^{sim}')
            ),
            h_mu_leading_d0 = cms.PSet(
                xAxis = cms.PSet(
                    d0binning,
                    Label = cms.string('leading muon d_{0}^{reco} [cm]')
                ),
                title = cms.string('leading muon d_{0}^{reco}')
            ),
            h_mu_leading_90pt = cms.PSet(
                xAxis = cms.PSet(
                    pTbinning,
                    Label = cms.string('leading muon 90% p_{T}^{reco} [GeV]')
                ),
                title = cms.string('leading muon 90% p_{T}^{reco}')
            ),
            h_mu_IsoDep = cms.PSet(
               xAxis = cms.PSet(
                   nBins = cms.uint32(100),
                   Max = cms.double(20.0),
                   Min = cms.double(0.0),
                   Label = cms.string('Energy deposition [GeV]')
               ),
               title = cms.string('Energy deposition [GeV]')
            ),
            h_mu_nIsoDep = cms.PSet(
               xAxis = cms.PSet(
                   nBins = cms.uint32(21),
                   Max = cms.double(20.5),
                   Min = cms.double(-0.5),
                   Label = cms.string('Number of object in deposit cone of muon')
               ),
               title = cms.string('Number of object in deposit cone of muon')
            ),
            h_mu_leading_IsoDep = cms.PSet(
               xAxis = cms.PSet(
                   nBins = cms.uint32(100),
                   Max = cms.double(20.0),
                   Min = cms.double(0.0),
                   Label = cms.string('Energy deposition around leading muon [GeV]')
               ),
               title = cms.string('Energy deposition around leading muon [GeV]')
            ),
            h_mu_leading_nIsoDep = cms.PSet(
               xAxis = cms.PSet(
                   nBins = cms.uint32(21),
                   Max = cms.double(20.5),
                   Min = cms.double(-0.5),
                   Label = cms.string('Number of object in deposit cone of leading muon')
               ),
               title = cms.string('Number of object in deposit cone of leading muon')
            )            
        ),
        H2s = cms.PSet(
            h_mu_d0_pT = cms.PSet(
                yAxis = cms.PSet(
                    d0binning,
                    Label = cms.string('muon d_{0}^{reco}')
                ),
                xAxis = cms.PSet(
                    pTbinning,
                    Label = cms.string('muon p_{T}^{reco} [GeV]')
                ),
                title = cms.string('muon d_{0}^{reco} versus p_{T}^{reco}')
            )
        )
    )
