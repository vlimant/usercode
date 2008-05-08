#ifndef NTupler_H
#define NTupler_H

#include "Workspace/ConfigurableAnalysis/interface/VariableHelper.h"
#include "Workspace/ConfigurableAnalysis/interface/UpdaterService.h"

#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/EDFilter.h"


/*
 * Description:
 * placeholder for common ntuplizer tools
 *
 */

class NTupler {
 public:
  NTupler(){}
  virtual ~NTupler(){}

  virtual uint registerleaves(edm::EDFilter * producer) =0;
  virtual void fill(edm::Event& iEvent)=0;
};


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

 private:
  std::string class_;
  edm::InputTag src_;
  std::string expr_;
  std::string maxIndexName_;
  std::string branchAlias_;
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
	  //	  try {
	    (*value_)[i]=(expr)((*oH)[i]); 
	    //	  }catch(...){	    edm::LogError("StringBranchHelper")<<"could not evaluate expression: "<<B.expr()<<" on class: "<<B.className(); }
	}
      }
    }
 private:
  value value_;
};

class VariableNTupler : public NTupler{
 public:
  VariableNTupler(const edm::ParameterSet& iConfig){
    if (iConfig.getParameter<bool>("allVariables"))
      {
	VariableHelper::iterator v=VariableHelperInstance::get().begin();
	VariableHelper::iterator v_end=VariableHelperInstance::get().end();
	for(;v!=v_end;++v){
	  leaves_[v->second->name()]=v->second;
	}
      }
    else{
      std::vector<std::string> leaves=iConfig.getParameter<std::vector<std::string> >("leaves");
      for (uint i=0;i!=leaves.size();++i){
	leaves_[leaves[i]]= VariableHelperInstance::get().variable(leaves[i]);
      }
    }
  }
  
  uint registerleaves(edm::EDFilter * producer){
    //loop the leaves registered
    uint nLeaves=0;
    iterator i=leaves_.begin();
    iterator i_end= leaves_.end();
    for(;i!=i_end;++i){
      nLeaves++;
      producer->produces<double>(i->first).setBranchAlias(i->first);
    }
    return nLeaves;
  }

  void fill(edm::Event& iEvent){
    //protection against stupid users :-(
    if (!edm::Service<UpdaterService>()->checkOnce("VariableNTupler::fill")) return;
    
    //other leaves
    iterator i=leaves_.begin();
    iterator i_end=leaves_.end();
    for(;i!=i_end;++i){
      std::auto_ptr<double> leafValue(new double((*i->second)()));
      iEvent.put(leafValue, i->first);
    }
  }
 protected:
  typedef std::map<std::string, const CachingVariable *>::iterator iterator;
  std::map<std::string, const CachingVariable *> leaves_;
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

  }

  uint registerleaves(edm::EDFilter * producer){
    uint nLeaves=0;
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
    return nLeaves;
  }

  void fill(edm::Event& iEvent){
    //protection against stupid users :-(
    if (!edm::Service<UpdaterService>()->checkOnce("StringBasedNTupler::fill")) return;

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
  protected:
  typedef std::map<std::string, std::vector<TreeBranch> > Branches;
  Branches branches_;


};



class CombinedNTupler : public NTupler{
 public:
  CombinedNTupler(const edm::ParameterSet& iConfig){
    all_.push_back(new VariableNTupler(iConfig));
    all_.push_back(new StringBasedNTupler(iConfig));
  }
  uint registerleaves(edm::EDFilter * producer){
    uint nLeaves=0;
    for (uint n=0;n!=all_.size();++n) nLeaves+=all_[n]->registerleaves(producer);
    return nLeaves;}

  void fill(edm::Event& iEvent){for (uint n=0;n!=all_.size();++n) all_[n]->fill(iEvent);}

 protected:
  std::vector<NTupler*> all_;
};


#endif
