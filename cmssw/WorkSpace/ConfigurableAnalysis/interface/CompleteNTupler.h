#include "PhysicsTools/UtilAlgos/interface/NTupler.h"

#include "Workspace/ConfigurableAnalysis/interface/StringBasedNTupler.h"
#include "Workspace/ConfigurableAnalysis/interface/VariableNTupler.h"
#include "Workspace/ConfigurableAnalysis/interface/AdHocNTupler.h"

class CompleteNTupler : public NTupler {
 public:
  CompleteNTupler(const edm::ParameterSet& iConfig){
    sN = new StringBasedNTupler(iConfig);
    vN = new VariableNTupler(iConfig);
    if (iConfig.exists("AdHocNPSet"))
      if (!iConfig.getParameter<edm::ParameterSet>("AdHocNPSet").empty())
	aN = new AdHocNTupler(iConfig);
      else aN=0;
    else
      aN=0;
  }
  
  uint registerleaves(edm::ProducerBase * producer){
    uint nLeaves=0;
    nLeaves+=sN->registerleaves(producer);
    nLeaves+=vN->registerleaves(producer);
    if (aN)
      nLeaves+=aN->registerleaves(producer);
    return nLeaves;
  }
  void fill(edm::Event& iEvent){
    sN->fill(iEvent);
    vN->fill(iEvent);
    if (aN)
      aN->fill(iEvent);

    sN->callBack();
    vN->callBack();
    if (aN)
      aN->callBack();
  }

 private:
  StringBasedNTupler * sN;
  VariableNTupler * vN;  
  AdHocNTupler * aN;

};

