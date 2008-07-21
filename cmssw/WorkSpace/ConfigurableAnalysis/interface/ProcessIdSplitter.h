
#include "PhysicsTools/UtilAlgos/interface/CachingVariable.h"
#include "PhysicsTools/HepMCCandAlgos/interface/CSA07ProcessId.h"


class ProcessIdSplitter : public Splitter{
 public:
  ProcessIdSplitter(std::string  n,const edm::ParameterSet & iConfig) :
    Splitter("ProcessIdSplitter",n,iConfig){
    lumi_=iConfig.getParameter<double>("lumi");
    weightLabel_=iConfig.getParameter<std::string>("weightLabel");
    uint maxID = iConfig.getParameter<uint>("maxID");//70
    //fill the labels
    for (uint id=0;id!=maxID;++id){
      labels_.push_back(std::string(csa07::csa07ProcessName(id)));
      std::stringstream ss;
      ss<<"_processIDSplit_"<<id;
      short_labels_.push_back(ss.str());
    }
  }

    CachingVariable::evalType eval(const edm::Event & iEvent) const{
      int ID=csa07::csa07ProcessId(iEvent,lumi_,weightLabel_);
      return std::make_pair(true, (double)ID);
    }

 private:
  double lumi_;
  std::string weightLabel_;
};
