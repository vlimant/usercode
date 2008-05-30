#include "Workspace/ConfigurableAnalysis/interface/CachingVariable.h"
#include "Workspace/ConfigurableAnalysis/interface/VariableHelper.h"

//#define  VHS edm::Service<VariableHelperService>()->get()
#define VHS VariableHelperInstance::get()

const edm::Event & CachingVariable::event() const {return VHS.event();}
const edm::EventSetup  & CachingVariable::setup() const {return VHS.setup();}

CachingVariable::evalType VarSplitter::eval() const{
  double v=(VHS)(var_);
  if (v<slots_.front()){
    if (useUnderFlow_) return std::make_pair(true,0);
    else return std::make_pair(false,0);
  }
  if (v>=slots_.back()){
    if (useOverFlow_) return std::make_pair(true,(double)maxIndex());
    else return std::make_pair(false,0);
  }
  uint i=1;
  for (;i<slots_.size();++i)
    if (v<slots_[i]) break;

  if (useUnderFlow_) return std::make_pair(true,(double) i);
  //need to substract 1 because checking on upper edges
  else return std::make_pair(true,(double)i-1);
}

#include "PhysicsTools/HepMCCandAlgos/interface/CSA07ProcessId.h"

ProcessIdSplitter::ProcessIdSplitter(std::string  n,const edm::ParameterSet & iConfig) :
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

CachingVariable::evalType ProcessIdSplitter::eval() const{
  int ID=csa07::csa07ProcessId(event(),lumi_,weightLabel_);
  return std::make_pair(true, (double)ID);
}

CachingVariable::evalType Power::eval() const {
  double v=(VHS)(var_);
  double p=exp(power_*log(v));
  return std::make_pair(true,p);
}


CachingVariable::evalType CalculateMHT::eval() const {
   edm::Handle<std::vector<pat::Jet> > jH;
    event().getByLabel(jetSrc_,jH);
    if (jH.failedToGet()) 
      return std::make_pair(false,0);

    edm::Handle<std::vector<pat::MET> > mH;
    event().getByLabel(metSrc_,mH);
    if (mH.failedToGet())
      return std::make_pair(false,0);

    double MHT=helper_.SUSYHT(*jH,(*mH)[0]);
    return std::make_pair(true,MHT);
}
