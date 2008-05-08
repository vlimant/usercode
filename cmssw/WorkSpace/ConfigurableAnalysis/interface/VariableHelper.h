#ifndef ConfigurableAnalysis_VariableHelper_H
#define ConfigurableAnalysis_VariableHelper_H

#include "FWCore/ServiceRegistry/interface/ActivityRegistry.h"
#include "Workspace/ConfigurableAnalysis/interface/CachingVariable.h"

namespace edm {
  class Event;
  class EventSetup;
}

class VariableHelper {
 public:
  VariableHelper(const edm::ParameterSet & iConfig);  
  typedef std::map<std::string, CachingVariable*>::const_iterator iterator;

  const CachingVariable* variable(std::string name)const ;
  double operator()(std::string & name) const;
  double operator()(std::string name) const;

  iterator begin() { return variables_.begin();}
  iterator end() { return variables_.end();}

  void update(const edm::Event & e, const edm::EventSetup & es) const;
  const edm::Event & event() const { if (!ev_) throw; return *ev_;}
  const edm::EventSetup & setup() const { if (!es_) throw; return *es_;}
 private:
  mutable const edm::Event * ev_;
  mutable const edm::EventSetup * es_;
  std::map<std::string, CachingVariable*> variables_;
};

class VariableHelperInstance {
 private:
  static VariableHelper * VariableHelperUniqueInstance_;
 public:
  static VariableHelper & init(const edm::ParameterSet & iConfig){
    if (VariableHelperUniqueInstance_) throw;
    else VariableHelperUniqueInstance_ = new VariableHelper(iConfig);
    return *VariableHelperUniqueInstance_;
  }
  static VariableHelper & get(){
    if (!VariableHelperUniqueInstance_) throw;
    else return (*VariableHelperUniqueInstance_);
  }
};

/*
class VariableHelperService {
 public:
  VariableHelperService(const edm::ParameterSet & iConfig, edm::ActivityRegistry &);
  const VariableHelper & get() const { return *instance_;}
 private:
  VariableHelper * instance_;
};
*/

#endif
