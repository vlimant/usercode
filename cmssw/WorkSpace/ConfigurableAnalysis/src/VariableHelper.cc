#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "Workspace/ConfigurableAnalysis/interface/VariableHelper.h"
#include "Workspace/ConfigurableAnalysis/interface/CachingVariable.h"
#include "Workspace/ConfigurableAnalysis/interface/CachingVariableFactory.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

//VariableHelperService::VariableHelperService(const edm::ParameterSet & iConfig, edm::ActivityRegistry & r){instance_ = new VariableHelper(iConfig);}

VariableHelper * VariableHelperInstance::VariableHelperUniqueInstance_=0;

VariableHelper::VariableHelper(const edm::ParameterSet & iConfig){
  std::vector<std::string> psetNames;
  iConfig.getParameterSetNames(psetNames);
  for (uint i=0;i!=psetNames.size();++i){
    std::string & vname=psetNames[i];
    edm::ParameterSet vPset=iConfig.getParameter<edm::ParameterSet>(psetNames[i]);
    //std::string type=vPset.getParameter<std::string>("type");
    std::string type="helper";
    if (type=="helper"){
      std::string method=vPset.getParameter<std::string>("method");
      variables_[vname]=CachingVariableFactory::get()->create(method,vname,vPset);
    }
    /*    else if (type=="parser"){
	  std::string Class=vPset.getParameter<std::string>("class");
	  if (Class=="pat::Jet")
	  {
	  variables_[vname]=new ExpressionVariable<pat::Jet>(vname,vPset);
	  }
	  else if (Class=="pat::Muon")
	  {
	  variables_[vname]=new ExpressionVariable<pat::Muon>(vname,vPset);
	  }
	  else
	  {
	  //Class not recognized
	  throw;
	  }
	}*/
    else{
      //type not recognized
      throw;
    }
  }
}

void VariableHelper::update(const edm::Event & e, const edm::EventSetup & es) const
{
  ev_=&e;
  es_=&es;
}

/*
const CachingVariable* VariableHelper::variable(std::string & name) const{
std::map<std::string, CachingVariable*> ::const_iterator v=variables_.find(name);
if (v!=variables_.end())
return v->second;
else
{
//unknown request variable.
throw;
}
}
*/

const CachingVariable* VariableHelper::variable(std::string name) const{ 
  std::map<std::string, CachingVariable*> ::const_iterator v=variables_.find(name);
  if (v!=variables_.end())
    return v->second;
  else
    {
      edm::LogError("VariableHelper")<<"I don't know anything named: "<<name;
      return 0;
    }
}


double VariableHelper::operator() (std::string & name) const{  
  const CachingVariable* v = variable(name);
  return (*v)();
}

double VariableHelper::operator() (std::string name) const{
  const CachingVariable* v = variable(name);
  return (*v)();
}

