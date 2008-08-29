#ifndef VariableNtupler_NTupler_H
#define VariableNtupler_NTupler_H

#include "PhysicsTools/UtilAlgos/interface/VariableHelper.h"
#include "PhysicsTools/UtilAlgos/interface/UpdaterService.h"

#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include <FWCore/Framework/interface/ProducerBase.h>

#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"

#include "PhysicsTools/UtilAlgos/interface/NTupler.h"

class VariableNTupler : public NTupler{
 public:
  VariableNTupler(const edm::ParameterSet& iConfig){
    ownTheTree_=false;
    edm::ParameterSet variablePSet=iConfig.getParameter<edm::ParameterSet>("variablesPSet");
    if (variablePSet.getParameter<bool>("allVariables"))
      {
	VariableHelper::iterator v=edm::Service<VariableHelperService>()->get().begin();
	VariableHelper::iterator v_end=edm::Service<VariableHelperService>()->get().end();
	for(;v!=v_end;++v){
	  leaves_[v->second->name()]=v->second;
	}
      }
    else{
      std::vector<std::string> leaves=variablePSet.getParameter<std::vector<std::string> >("leaves");
      for (uint i=0;i!=leaves.size();++i){
	leaves_[leaves[i]]= edm::Service<VariableHelperService>()->get().variable(leaves[i]);
      }
    }
    if (variablePSet.exists("useTFileService"))
      useTFileService_=variablePSet.getParameter<bool>("useTFileService");
    else
      useTFileService_=iConfig.getParameter<bool>("useTFileService");

    if (useTFileService_){
      if (variablePSet.exists("treeName"))
	treeName_=variablePSet.getParameter<std::string>("treeName");
      else
	treeName_=iConfig.getParameter<std::string>("treeName");
    }
  }
  
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
      //      tree_=dynamic_cast<TTree*>(fs->file().FindObjectAny(treeName_.c_str()));
      //      if (!tree_){
      //      std::cout<<"VariableNTupler owns its tree"<<std::endl;
	ownTheTree_=true;
	tree_=fs->make<TTree>(treeName_.c_str(),"VariableNTuple tree");
	//	fs->file().Add(tree_);
	//      }
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
    //    if (!edm::Service<UpdaterService>()->checkOnce("VariableNTupler::fill")) return;
    //    std::cout<<"I am trying to fill the tree VariableNTupler "<< useTFileService_ <<" " <<ownTheTree_<< std::endl;
    if (useTFileService_){
      //fill the data holder
      iterator i=leaves_.begin();
      iterator i_end=leaves_.end();
      uint iInDataHolder=0;
      //      std::cout<<"looping: "<<leaves_.size()<<std::endl;
      for(;i!=i_end;++i,++iInDataHolder){
	//	std::cout<<"index: "<<iInDataHolder<<std::endl;
	dataHolder_[iInDataHolder]=(*i->second)(iEvent);
	//	std::cout<<"succeeded"<<std::endl;
      }
      //fill into root;
      if (ownTheTree_) {
	//	std::cout<<"I am filling the tree VariableNTupler"<<std::endl;
	tree_->Fill();
      }
    }else{
      //other leaves
      iterator i=leaves_.begin();
      iterator i_end=leaves_.end();
      for(;i!=i_end;++i){
	std::auto_ptr<double> leafValue(new double((*i->second)(iEvent)));
	iEvent.put(leafValue, i->first);
      }
    }
  }
 protected:
  typedef std::map<std::string, const CachingVariable *>::iterator iterator;
  std::map<std::string, const CachingVariable *> leaves_;

  bool ownTheTree_;
  std::string treeName_;
  double * dataHolder_;
};


#endif
