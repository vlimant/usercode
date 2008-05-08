
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "FWCore/Framework/interface/Event.h"
#include <vector>

class SUSYHelper
{
 public:

  float  SUSYEMETIndirectLepton(const edm::Event& iEvent) const;

  float SUSYEMETCorMuons(const pat::MET&  myMET,const std::vector< pat::Muon>& myMuons) const;
  float SUSYHT( const std::vector< pat::Jet >& myJets,const pat::MET&  myMET) const;
  float SUSYMETIso(const std::vector< pat::Jet >& myJets,const pat::MET&  myMET) const;
  float SUSYMeff( const std::vector< pat::Jet >& myJets,const pat::MET&  myMET) const;
  // calculated  EMFrack of event a la TDRII
  float SUSYETEMFrack( const std::vector< pat::Jet >& myJets) const;
  // calculated  Charge Frack of event a la TDRII
  float SUSYEChFrack( const std::vector< pat::Jet >& myJets) const;
  float SUSYEChFrack( const  pat::Jet& aJet) const;

  float SUSYAllGenRecoilMET(const std::vector< pat::Jet >& myJets,const std::vector< pat::Muon>& myMuons ,const pat::MET&  myMET) const;
  float SUSYMETRXY(const std::vector< pat::Jet >& myJets,const pat::MET&  myMET,unsigned int x,unsigned int y) const;
 
  // add muon momentum to MET from MET muons
  reco::Particle::Vector MuonCorrection( reco::Particle::Vector METmometum, const std::vector< pat::Muon>& myMuons) const;

  // loop over all jets, returns - summed momentum vector
  reco::Particle::Vector SUSYRecoilMET(const std::vector< pat::Jet >& myJets) const;
  // loop over jets within pt&eta range, returns - summed momentum vector
  // fake jets are eventually rejected
  reco::Particle::Vector SUSYRecoilMETCutted(const std::vector< pat::Jet >& myJets, float ptMin, float etaMax) const;
  // loop over over jets & muons,  returns - summed momentum vector
  reco::Particle::Vector SUSYAllRecoilMET(const std::vector< pat::Jet >& myJets,const std::vector< pat::Electron >& myElectrons,const std::vector< pat::Muon>& myMuons ,const std::vector< pat::Tau>&  myTau) const;
};
