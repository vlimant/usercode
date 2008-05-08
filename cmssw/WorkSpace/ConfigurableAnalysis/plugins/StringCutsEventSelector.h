#ifndef _StringCutsEventSelector_H
#define _StringCutsEventSelector_H

#include "Workspace/EventSelectors/interface/SusyEventSelector.h"
#include "PhysicsTools/Utilities/interface/StringCutObjectSelector.h"

template<typename Object>
class  StringCutsEventSelector : public SusyEventSelector {
 public:
  StringCutsEventSelector(const edm::ParameterSet& pset) :
    SusyEventSelector(pset),
    src_(pset.getParameter<edm::InputTag>("src"))
      {
	std::vector<std::string> selection=pset.getParameter<std::vector<std::string > >("cut");
	std::stringstream ss;
	ss<<"string cut based selection on collection: "<<src_;
	description_.push_back(ss.str());	    ss.str("");
	description_.push_back("selection cuts are:");
	for (uint i=0;i!=selection.size();i++)
	  if (selection[i]!="-"){
	    f_.push_back( new StringCutObjectSelector<Object>(selection[i]));
	    ss<<"["<<i<<"]: "<<selection[i];
	    description_.push_back(ss.str());           ss.str("");
	  }
	  else
	    {
	      f_.push_back(0);
	      ss<<"["<<i<<"]: no selection";
	      description_.push_back(ss.str());           ss.str("");
	    }
      }
    
    bool select (const edm::Event& e) const{
      edm::Handle<std::vector<Object> > oH;
      e.getByLabel(src_, oH);
      uint i=0;
      if (oH->size()<f_.size()) return false;
      for (;i!=f_.size();i++)
	{  
	  if (!f_[i]) continue;
	  const Object & o = (*oH)[i];
	  if (!(*f_[i])(o)) return false;
	}
      return true;
    }
    
 private:
    edm::InputTag src_;
    std::vector<StringCutObjectSelector<Object> *> f_;
};

#endif
