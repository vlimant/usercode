#ifndef _StringCutEventSelector
#define _StringCutEventSelector

#include "Workspace/EventSelectors/interface/SusyEventSelector.h"
#include "PhysicsTools/Utilities/interface/StringCutObjectSelector.h"
#include "Workspace/ConfigurableAnalysis/interface/InputTagDistributor.h"

template<typename Object>
class  StringCutEventSelector : public SusyEventSelector {
 public:
  StringCutEventSelector(const edm::ParameterSet& pset) :
    SusyEventSelector(pset),
    //    src_(pset.getParameter<edm::InputTag>("src")),
    src_(InputTagDistributor::retrieve("src",pset)),
    f_(pset.getParameter<std::string>("cut")),
    //put this guy to 0 to do the check on "all" object in the collection
    nFirst_(pset.getParameter<uint>("nFirst"))
      {
	  std::stringstream ss;
	  ss<<"string cut based selection on collection: "<<src_;
	  description_.push_back(ss.str());
	  ss.str("");
	  description_.push_back(std::string("selection cut is: ")+pset.getParameter<std::string>("cut"));
      }
    
    bool select (const edm::Event& e) const{
      edm::Handle<std::vector<Object> > oH;
      e.getByLabel(src_, oH);
      //reject events if not enough object in collection
      //      if ((nFirst_!=0) && (oH->size()<nFirst_)) return false;
      uint i=0;
      for (;i!=oH->size();i++)
	{
	  //stop doing the check if reaching too far in the collection
	  if ((nFirst_!=0) && (i>=nFirst_)) break;
	  const Object & o = (*oH)[i];
	  if (!f_(o)) return false;
	}
      return true;
    }
    
 private:
    edm::InputTag src_;
    StringCutObjectSelector<Object> f_;
    uint nFirst_;
};

#endif
