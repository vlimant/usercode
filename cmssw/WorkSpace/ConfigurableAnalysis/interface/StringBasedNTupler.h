#ifndef StringBasedNTupler_NTupler_H
#define StringBasedNTupler_NTupler_H

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

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "PhysicsTools/Utilities/interface/StringObjectFunction.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


class TreeBranch {
 public:
  TreeBranch(): class_(""),expr_(""),maxIndexName_(""),branchAlias_("") {}
    TreeBranch(std::string C, edm::InputTag S, std::string E, std::string Mi, std::string Ba) :
      class_(C),src_(S),expr_(E),maxIndexName_(Mi),branchAlias_(Ba){}
    
  const std::string & className() const { return class_;}
  const edm::InputTag & src() const { return src_;}
  const std::string & expr() const { return expr_;}
  const std::string & maxIndexName() const { return maxIndexName_;}
  const std::string branchName()const{ 
	TString name(branchAlias_);
	name.ReplaceAll("_","0");
	return std::string(name.Data());}
  const std::string & branchAlias()const{ return branchAlias_;}

  typedef std::auto_ptr<std::vector<double> > value;
  value branch(const edm::Event& iEvent);

  std::vector<double>** dataHolderPtrAdress() { return &dataHolderPtr_;}
  std::vector<double>* dataHolderPtr() { return dataHolderPtr_;}
  void assignDataHolderPtr(std::vector<double> * data) { dataHolderPtr_=data;}
 private:
  std::string class_;
  edm::InputTag src_;
  std::string expr_;
  std::string maxIndexName_;
  std::string branchAlias_;

  std::vector<double> * dataHolderPtr_;
};

template <typename Object, typename Collection=std::vector<Object> >
class StringBranchHelper {
public:
  typedef TreeBranch::value value;
  value operator()() { return value_;}

  StringBranchHelper(const TreeBranch & B, const edm::Event& iEvent)
    {
      //    grab the collection
      edm::Handle<Collection> oH;
      iEvent.getByLabel(B.src(), oH);
      //empty vector if product not found
      if (oH.failedToGet()){
	edm::LogError("StringBranchHelper")<<"cannot open: "<<B.src();
	std::auto_ptr<std::vector<double> > ret(new std::vector<double>());
      }
      else{
	//parser for the object expression
	StringObjectFunction<Object> expr(B.expr());
	value_.reset(new std::vector<double>(oH->size()));

	//actually fill the vector of values
	uint i_end=oH->size();
	for (uint i=0;i!=i_end;++i){
	  try {
	    (*value_)[i]=(expr)((*oH)[i]); 
	  }catch(...){
	    edm::LogError("StringBranchHelper")<<"could not evaluate expression: "<<B.expr()<<" on class: "<<B.className(); }
	}
      }
    }
 private:
  value value_;
};



class StringBasedNTupler : public NTupler {
 public:
  StringBasedNTupler(const edm::ParameterSet& iConfig){
    edm::ParameterSet branchesPSet = iConfig.getParameter<edm::ParameterSet>("branchesPSet");
    std::vector<std::string> branches;
    branchesPSet.getParameterSetNames(branches);
    for (uint b=0;b!=branches.size();++b){
      edm::ParameterSet bPSet = branchesPSet.getParameter<edm::ParameterSet>(branches[b]);
      std::string className=bPSet.getParameter<std::string>("class");
      edm::InputTag src=bPSet.getParameter<edm::InputTag>("src");
      edm::ParameterSet leavesPSet=bPSet.getParameter<edm::ParameterSet>("leaves");
      std::vector<std::string> leaves=leavesPSet.getParameterNamesForType<std::string>();
      std::string maxName="N"+branches[b];
      for (uint l=0;l!=leaves.size();++l){
	std::string leave_expr=leavesPSet.getParameter<std::string>(leaves[l]);
	std::string branchAlias=branches[b]+"_"+leaves[l];
	branches_[maxName].push_back(TreeBranch(className, src, leave_expr, maxName, branchAlias));
      }//loop the provided leaves
    }//loop the provided branches
    if (branchesPSet.exists("useTFileService"))
      useTFileService_=branchesPSet.getParameter<bool>("useTFileService");         
    else
      useTFileService_=iConfig.getParameter<bool>("useTFileService");
    treeName_=branchesPSet.getParameter<std::string>("treeName");
  }

  //  uint registerleaves(edm::EDFilter * producer){
  uint registerleaves(edm::ProducerBase * producer){
    uint nLeaves=0;

    if (useTFileService_){
      edm::Service<TFileService> fs;
      tree_=fs->make<TTree>(treeName_.c_str(),"StringBasedNTupler tree");
      //reserve memory for the indexes      
      indexDataHolder_ = new uint[branches_.size()];
      // loop the automated leafer
      Branches::iterator iB=branches_.begin();
      Branches::iterator iB_end=branches_.end();
      uint indexOfIndexInDataHolder=0;
      for(;iB!=iB_end;++iB,++indexOfIndexInDataHolder){
	//branch for the index
	//---	std::cout<<"happy so far 1"<<std::endl;
	tree_->Branch(iB->first.c_str(), &(indexDataHolder_[indexOfIndexInDataHolder]),(iB->first+"/i").c_str());
	//---	std::cout<<"happy so far 2"<<std::endl;
	std::vector<TreeBranch>::iterator iL=iB->second.begin();
	std::vector<TreeBranch>::iterator iL_end=iB->second.end();
	for(;iL!=iL_end;++iL){
	  TreeBranch & b=*iL;
	  //---	  std::cout<<"happy so far 3"<<std::endl;
	  //branch for the leave
	  tree_->Branch(b.branchAlias().c_str(),"std::vector<double>",iL->dataHolderPtrAdress());
	  //---	  std::cout<<"happy so far 4"<<std::endl;
	}
      }
    }
    else{
      // loop the automated leafer
      Branches::iterator iB=branches_.begin();
      Branches::iterator iB_end=branches_.end();
      for(;iB!=iB_end;++iB){
	//the index. should produce it only once
	nLeaves++;
	producer->produces<uint>(iB->first).setBranchAlias(iB->first);
	std::vector<TreeBranch>::iterator iL=iB->second.begin();
	std::vector<TreeBranch>::iterator iL_end=iB->second.end();
	for(;iL!=iL_end;++iL){
	  TreeBranch & b=*iL;
	  //the values
	  producer->produces<std::vector<double> >(b.branchName()).setBranchAlias(b.branchAlias());
	}
      }
    }
    return nLeaves;
  }

  void fill(edm::Event& iEvent){
    //protection against stupid users :-(
    //    if (!edm::Service<UpdaterService>()->checkOnce("StringBasedNTupler::fill")) return;
    //well if you do that, you cannot have two ntupler of the same type in the same job...

    if (useTFileService_){
      // loop the automated leafer
      Branches::iterator iB=branches_.begin();
      Branches::iterator iB_end=branches_.end();
      uint indexOfIndexInDataHolder=0;
      for(;iB!=iB_end;++iB,++indexOfIndexInDataHolder){
	std::vector<TreeBranch>::iterator iL=iB->second.begin();
	std::vector<TreeBranch>::iterator iL_end=iB->second.end();
	uint maxS=0;
	for(;iL!=iL_end;++iL){
	  TreeBranch & b=*iL;
	  //---	  std::cout<<"happy so far 5"<<std::endl;
	  std::auto_ptr<std::vector<double> > branch(b.branch(iEvent));
	  //---	  std::cout<<"happy so far 6"<<std::endl;
	  if (branch->size()>maxS) maxS=branch->size();
	  //---	  std::cout<<"happy so far 7"<<std::endl;
	  b.assignDataHolderPtr(branch.release());
	  //---	  std::cout<<"happy so far 8"<<std::endl;
	}
	indexDataHolder_[indexOfIndexInDataHolder]=maxS;
      }
      tree_->Fill();

      //de-allocate memory now: allocated in branch(...) and released to the pointer.
      for(;iB!=iB_end;++iB,++indexOfIndexInDataHolder){
	std::vector<TreeBranch>::iterator iL=iB->second.begin();
	std::vector<TreeBranch>::iterator iL_end=iB->second.end();
	for(;iL!=iL_end;++iL){
          TreeBranch & b=*iL;
	  delete b.dataHolderPtr();
	}}
    }else{
      // loop the automated leafer
      Branches::iterator iB=branches_.begin();
      Branches::iterator iB_end=branches_.end();
      for(;iB!=iB_end;++iB){
	std::vector<TreeBranch>::iterator iL=iB->second.begin();
	std::vector<TreeBranch>::iterator iL_end=iB->second.end();
	uint maxS=0;
	for(;iL!=iL_end;++iL){
	  TreeBranch & b=*iL;
	  std::auto_ptr<std::vector<double> > branch(b.branch(iEvent));
	  if (branch->size()>maxS) maxS=branch->size();
	  iEvent.put(branch, b.branchName());
	}
	//index should be put only once per branch
	//FIXME. max or min of branch size?
	std::auto_ptr<uint> maxN(new uint(maxS));
	iEvent.put(maxN, iB->first);
    }
    }
  }

  protected:
  typedef std::map<std::string, std::vector<TreeBranch> > Branches;
  Branches branches_;

  std::string treeName_;
  uint * indexDataHolder_;
};


#endif
