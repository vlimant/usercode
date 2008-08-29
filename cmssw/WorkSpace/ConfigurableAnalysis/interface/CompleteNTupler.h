#include "PhysicsTools/UtilAlgos/interface/NTupler.h"

#include "Workspace/ConfigurableAnalysis/interface/StringBasedNTupler.h"
#include "Workspace/ConfigurableAnalysis/interface/VariableNTupler.h"

class CompleteNTupler : public NTupler {
 public:
  CompleteNTupler(const edm::ParameterSet& iConfig){
    sN = new StringBasedNTupler(iConfig);
    vN = new VariableNTupler(iConfig);
  }
  
  uint registerleaves(edm::ProducerBase * producer){
    uint nLeaves=0;
    nLeaves+=sN->registerleaves(producer);
    nLeaves+=vN->registerleaves(producer);
    return nLeaves;
  }
  void fill(edm::Event& iEvent){
    sN->fill(iEvent);
    vN->fill(iEvent);
  }

 private:
  StringBasedNTupler * sN;
  VariableNTupler * vN;  

};

