#ifndef VariableNtupler_NTupler_H
#define VariableNtupler_NTupler_H

#include "Workspace/ConfigurableAnalysis/interface/VariableHelper.h"
#include "Workspace/ConfigurableAnalysis/interface/UpdaterService.h"

#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include <FWCore/Framework/interface/ProducerBase.h>

#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
#include "TBranch.h"

#include "Workspace/ConfigurableAnalysis/interface/NTupler.h"

class VariableNTupler : public NTupler{
 public:
  VariableNTupler(const edm::ParameterSet& iConfig){
    edm::ParameterSet variablePSet=iConfig.getParameter<edm::ParameterSet>("variablesPSet");
    if (variablePSet.getParameter<bool>("allVariables"))
      {
	VariableHelper::iterator v=VariableHelperInstance::get().begin();
	VariableHelper::iterator v_end=VariableHelperInstance::get().end();
	for(;v!=v_end;++v){
	  leaves_[v->second->name()]=v->second;
	}
      }
    else{
      std::vector<std::string> leaves=variablePSet.getParameter<std::vector<std::string> >("leaves");
      for (uint i=0;i!=leaves.size();++i){
	leaves_[leaves[i]]= VariableHelperInstance::get().variable(leaves[i]);
      }
    }
    if (variablePSet.exists("useTFileService"))
      useTFileService_=variablePSet.getParameter<bool>("useTFileService");
    else
      useTFileService_=iConfig.getParameter<bool>("useTFileService");
  }
  
  //  uint registerleaves(edm::EDFilter * producer){
  uint registerleaves(edm::ProducerBase * producer){
    uint nLeaves=0;
    if (useTFileService_){
      //loop the leaves registered
      nLeaves=leaves_.size();
      // make arrays of pointer to the actual values
      dataHolder_=new double[nLeaves];
      iterator i=leaves_.begin();
      iterator i_end= leaves_.end();
      edm::Service<TFileService> fs;
      tree_=fs->make<TTree>("variable","VariableNTuple tree");
      uint iInDataHolder=0;
      for(;i!=i_end;++i,++iInDataHolder){
	tree_->Branch(i->first.c_str(), &(dataHolder_[iInDataHolder]), (i->first+"/D").c_str());
      }
    }else{
      //loop the leaves registered
      iterator i=leaves_.begin();
      iterator i_end= leaves_.end();
      for(;i!=i_end;++i){
	nLeaves++;
	producer->produces<double>(i->first).setBranchAlias(i->first);
      }
    }
    return nLeaves;
  }
  
  void fill(edm::Event& iEvent){
    //protection against stupid users :-(
    //    if (!edm::Service<UpdaterService>()->checkOnce("VariableNTupler::fill")) return;
    
    if (useTFileService_){
      //fill the data holder
      iterator i=leaves_.begin();
      iterator i_end=leaves_.end();
      uint iInDataHolder=0;
      for(;i!=i_end;++i,++iInDataHolder){
	dataHolder_[iInDataHolder]=(*i->second)();
      }
      //fill into root;
      tree_->Fill();
    }else{
      //other leaves
      iterator i=leaves_.begin();
      iterator i_end=leaves_.end();
      for(;i!=i_end;++i){
	std::auto_ptr<double> leafValue(new double((*i->second)()));
	iEvent.put(leafValue, i->first);
      }
    }
  }
 protected:
  typedef std::map<std::string, const CachingVariable *>::iterator iterator;
  std::map<std::string, const CachingVariable *> leaves_;

  double * dataHolder_;
};


#endif
