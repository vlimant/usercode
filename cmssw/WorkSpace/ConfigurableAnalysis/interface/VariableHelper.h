#ifndef ConfigurableAnalysis_VariableHelper_H
#define ConfigurableAnalysis_VariableHelper_H

#include "FWCore/ServiceRegistry/interface/ActivityRegistry.h"
#include "Workspace/ConfigurableAnalysis/interface/CachingVariable.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

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
  const edm::Event & event() const { if (!ev_){
      std::cerr<<"event has not been set. expect seg fault."<<std::endl;
      throw;}
    return *ev_;}
  const edm::EventSetup & setup() const { if (!es_){
      std::cerr<<"setup has not been set. expect seg fault."<<std::endl;
      throw;}
 return *es_;}
  void setHolder(std::string hn);
 private:
  mutable const edm::Event * ev_;
  mutable const edm::EventSetup * es_;
  std::map<std::string, CachingVariable*> variables_;
};

class VariableHelperInstance {
 private:
  static VariableHelper * SetVariableHelperUniqueInstance_;
  static std::map<std::string, VariableHelper* > multipleInstance_;

 public:

  static VariableHelper & init(std::string user, const edm::ParameterSet & iConfig){
    if (multipleInstance_.find(user)!=multipleInstance_.end()){
      std::cerr<<user<<" VariableHelper user already defined."<<std::endl;
      throw;}
    else SetVariableHelperUniqueInstance_ = new VariableHelper(iConfig);
    multipleInstance_[user] = SetVariableHelperUniqueInstance_;
    SetVariableHelperUniqueInstance_->setHolder(user);
    return (*SetVariableHelperUniqueInstance_);
  }
  
  static VariableHelper & get(){
    if (!SetVariableHelperUniqueInstance_)
      //    if (!VariableHelperUniqueInstance_ && !SetVariableHelperUniqueInstance_)
      {
	std::cerr<<" none of VariableHelperUniqueInstance_ or SetVariableHelperUniqueInstance_ is valid."<<std::endl;
	throw;
      }
    //    else {
    //      if (VariableHelperUniqueInstance_) return (*VariableHelperUniqueInstance_);
      else return (*SetVariableHelperUniqueInstance_);
    //    }
  }

  static VariableHelper & set(std::string user){
    std::map<std::string, VariableHelper* >::iterator f=multipleInstance_.find(user);
    if (f == multipleInstance_.end()){
      std::cerr<<user<<" VariableHelper user not defined."<<std::endl;
      throw;
    }
    else{
      SetVariableHelperUniqueInstance_ = (f->second);
      return (*SetVariableHelperUniqueInstance_);
    }
  }
};

#endif
