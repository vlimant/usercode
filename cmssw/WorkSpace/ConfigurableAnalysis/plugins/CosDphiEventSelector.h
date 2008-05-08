#ifndef Workspace_CosDphiEventSelector_h_
#define Workspace_CosDphiEventSelector_h_

// system include files
#include <memory>

// user include files
#include "Workspace/EventSelectors/interface/SusyEventSelector.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include <vector>
#include <string>

template< typename LHS, typename RHS>
class CosDphiEventSelector : public SusyEventSelector {
public:
  CosDphiEventSelector (const edm::ParameterSet& pset) :
    srcLhs_(pset.getParameter<edm::InputTag>("srcLhs")),
    indexLhs_(pset.getParameter<uint>("indexLhs")),
    srcRhs_(pset.getParameter<edm::InputTag>("srcRhs")),
    indexRhs_(pset.getParameter<uint>("indexRhs")),
    minCosDphi_(pset.getParameter<double>("minCosDPhi")),
    maxCosDphi_(pset.getParameter<double>("maxCosDPhi"))
      {
	std::stringstream ss;
	description_.push_back("cos delta phi selection, satisfying: ");
	ss<<minCosDphi_<<" < Cos(DeltaPhi( Obj1, Oj2 )) < "<<maxCosDphi_;
	description_.push_back(ss.str());	ss.str("");
	ss<<"with Obj1 at index: "<<indexLhs_<<" of: "<<srcLhs_;
	description_.push_back(ss.str());	ss.str("");
	ss<<"with Obj2 at index: "<<indexRhs_<<" of: "<<srcRhs_;
	description_.push_back(ss.str());	ss.str("");
      }
    

  virtual bool select (const edm::Event& e) const{
    edm::Handle<std::vector<LHS> > lhsH;
    e.getByLabel(srcLhs_, lhsH);
    if (lhsH.failedToGet()){
      edm::LogError("CosDphiEventSelector")<<name()<<" could not get a collection with label: "<<srcLhs_;
      return false;
    }
    if (indexLhs_>=lhsH->size()){
      edm::LogError("CosDphiEventSelector")<<name()<<" tries to access index: "<<indexLhs_<<" of: "<<srcLhs_<<" with: "<<lhsH->size()<<" entries.";
      return false;}
    const LHS & lhs = (*lhsH)[indexLhs_];
    
    edm::Handle<std::vector<RHS> > rhsH;
    e.getByLabel(srcRhs_, rhsH);
    if (rhsH.failedToGet()){
      edm::LogError("CosDphiEventSelector")<<name()<<" could not get a collection with label: "<<srcLhs_;
      return false;
    }
    if (indexRhs_>=rhsH->size()){
      edm::LogError("CosDphiEventSelector")<<name()<<" tries to access index: "<<indexRhs_<<" of: "<<srcRhs_<<" with: "<<rhsH->size()<<" entries.";
      return false;}
    const RHS & rhs = (*rhsH)[indexRhs_];

    double cdphi = cos(lhs.phi()-rhs.phi());

    if (cdphi>minCosDphi_ && cdphi<maxCosDphi_) return true;
    return false;
    }
  virtual ~CosDphiEventSelector () {}
private:
  edm::InputTag srcLhs_;
  uint indexLhs_;
  edm::InputTag srcRhs_;
  uint indexRhs_;
  double minCosDphi_;
  double maxCosDphi_;
};

#endif
