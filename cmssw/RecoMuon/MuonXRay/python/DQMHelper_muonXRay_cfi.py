import FWCore.ParameterSet.Config as cms
from RecoMuon.MuonXRay.DQMHelper_binning_cfi import *

DQMHelper_muonXRay = cms.PSet(
        H1s = cms.PSet(
            h_leading_mu_sim_Aeta = cms.PSet(
                xAxis = cms.PSet(
                    Aetabinning,
                    Label = cms.string('leading muon  |#eta^{sim}|')
                ),
                title = cms.string('leading muon  |#eta^{sim}|')
            ),
            h_mu_sim_eta_tricky = cms.PSet(
                xAxis = cms.PSet(
                    etabinning,
                    Label = cms.string('muon |#eta|^{sim}')
                ),
                title = cms.string('|#eta|^{sim} of muon for tricky events')
            ),
            h_leading_mu_sim_vertex_position = cms.PSet(
                nBins = cms.uint32(200),
                title = cms.string('ledaing muon production vertex position'),
                Min = cms.double(0.0),
                Max = cms.double(400.0),
                Label = cms.string('d(origin, decay_{vertex}) [cm]')
            ),
            h_leading_mu_sim_phi = cms.PSet(
                xAxis = cms.PSet(
                    phibinning,
                    Label = cms.string('leading muon  #varphi^{sim}')
                ),
                title = cms.string('leading muon  #varphi^{sim}')
            ),
            h_mu_sim_phi = cms.PSet(
                xAxis = cms.PSet(
                    phibinning,
                    Label = cms.string('muon #varphi^{sim}')
                ),
                title = cms.string('muon #varphi^{sim}')
            ),
            h_mu_sim_pt = cms.PSet(
                xAxis = cms.PSet(
                    pTbinning,
                    Label = cms.string('muon p_{T}^{sim} [GeV]')
                ),
                title = cms.string('muon p_{T}^{sim}')
            ),
            h_num_sim = cms.PSet(
                nBins = cms.uint32(21),
                title = cms.string('number of sim muons.'),
                Min = cms.double(-0.5),
                Max = cms.double(20.5),
                Label = cms.string('Number of sim muons.')
            ),
            h_parent_id = cms.PSet(
                nBins = cms.uint32(16),
                title = cms.string('Name of parent particle from tight association'),
                Min = cms.double(-0.5),
                Max = cms.double(15.5),
                Label = cms.string('category.')
            ),
            h_mu_sim_pt_tricky = cms.PSet(
                xAxis = cms.PSet(
                    pTbinning,
                    Label = cms.string('muon p_{T}^{sim} [GeV]')
                ),
                title = cms.string('pt of sim muons for tricky events')
            ),
            h_leading_mu_sim_pt = cms.PSet(
                xAxis = cms.PSet(
                    pTbinning,
                    Label = cms.string('leading muon  p_{T}^{sim} [GeV]')
                ),
                title = cms.string('leading muon  p_{T}^{sim}')
            ),
            h_leading_mu_sim_pt_tricky = cms.PSet(
                xAxis = cms.PSet(
                    pTbinning,
                    Label = cms.string('leading muon  p_{T}^{sim} [GeV]')
                ),
                title = cms.string('pt of sim leading muons for tricky events')
            ),
            h_sim_vertex_position = cms.PSet(
                nBins = cms.uint32(200),
                title = cms.string('production vertex position'),
                Min = cms.double(0.0),
                Max = cms.double(400.0),
                Label = cms.string('d(origin, decay_{vertex}) [cm]')
            ),
            h_leading_mu_sim_eta_tricky = cms.PSet(
                xAxis = cms.PSet(
                    etabinning,
                    Label = cms.string('leading muon  |#eta|^{sim}')
                ),
                title = cms.string('|#eta|^{sim} of leading muon for tricky events')
            ),
            h_mu_sim_eta = cms.PSet(
                xAxis = cms.PSet(
                    etabinning,
                    Label = cms.string('muon #eta^{sim}')
                ),
                title = cms.string('muon #eta^{sim}')
            ),
            h_leading_mu_sim_eta = cms.PSet(
                xAxis = cms.PSet(
                    etabinning,
                    Label = cms.string('leading muon  #eta^{sim}')
                ),
                title = cms.string('leading muon  #eta^{sim}')
            ),
            h_mu_sim_Aeta = cms.PSet(
                xAxis = cms.PSet(
                    Aetabinning,
                    Label = cms.string('muon |#eta^{sim}|')
                ),
                title = cms.string('muon |#eta^{sim}|')
            )
        ),
        H2s = cms.PSet()
    )
