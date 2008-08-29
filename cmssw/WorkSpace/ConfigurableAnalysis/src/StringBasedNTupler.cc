#include "Workspace/ConfigurableAnalysis/interface/StringBasedNTupler.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Hemisphere.h"

#include "SimDataFormats/Track/interface/SimTrack.h"

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
  else ANOTHER_CLASS(reco::Muon);
  else ANOTHER_CLASS(reco::Track);
  else ANOTHER_CLASS(SimTrack);
  else {
    edm::LogError("TreeBranch")<<branchName()<<" failed to recognized class type: "<<class_;
    return TreeBranch::value(new std::vector<double>());
  }
}
#undef ANOTHER_CLASS
