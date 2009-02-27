#include "Workspace/ConfigurableAnalysis/interface/StringBasedNTupler.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Hemisphere.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include <DataFormats/BeamSpot/interface/BeamSpot.h>

#include "SimDataFormats/Track/interface/SimTrack.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"


//--------------------------------------------------------------------------------
//just define here a list of objects you would like to be able to have a branch of
//--------------------------------------------------------------------------------
#define ANOTHER_CLASS(C) if (class_==#C) return StringBranchHelper<C>(*this, iEvent)()

TreeBranch::value TreeBranch::branch(const edm::Event& iEvent){
  ANOTHER_CLASS(pat::Jet);
  else ANOTHER_CLASS(pat::Muon);
  else ANOTHER_CLASS(pat::Electron);
  else ANOTHER_CLASS(pat::MET);
  else ANOTHER_CLASS(pat::Tau);
  else ANOTHER_CLASS(pat::Hemisphere);
  else ANOTHER_CLASS(pat::Photon);
  else ANOTHER_CLASS(reco::Muon);
  else ANOTHER_CLASS(reco::Track);
  else ANOTHER_CLASS(SimTrack);
  else ANOTHER_CLASS(reco::GenParticle);
  else ANOTHER_CLASS(l1extra::L1ParticleMap);
  else ANOTHER_CLASS(reco::Vertex);
  else ANOTHER_CLASS(pat::GenericParticle);
  else {
    edm::LogError("TreeBranch")<<branchName()<<" failed to recognized class type: "<<class_;
    return TreeBranch::value(new std::vector<float>());
  }
}
#undef ANOTHER_CLASS
