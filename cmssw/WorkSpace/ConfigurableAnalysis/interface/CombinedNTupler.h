#ifndef CombinedNTupler_H
#define CombinedNTupler_H

#include "Workspace/ConfigurableAnalysis/interface/NTupler.h"
#include "Workspace/ConfigurableAnalysis/interface/VariableNTupler.h"
#include "Workspace/ConfigurableAnalysis/interface/StringBasedNTupler.h"

class CombinedNTupler : public NTupler{
 public:
  CombinedNTupler(const edm::ParameterSet& iConfig){
    all_.push_back(new VariableNTupler(iConfig));
    all_.push_back(new StringBasedNTupler(iConfig));
  }

  //  uint registerleaves(edm::EDFilter * producer){
  uint registerleaves(edm::ProducerBase * producer){
    uint nLeaves=0;
    for (uint n=0;n!=all_.size();++n) nLeaves+=all_[n]->registerleaves(producer);
    return nLeaves;}

  void fill(edm::Event& iEvent){for (uint n=0;n!=all_.size();++n) all_[n]->fill(iEvent);}

 protected:
  std::vector<NTupler*> all_;
};

#endif
