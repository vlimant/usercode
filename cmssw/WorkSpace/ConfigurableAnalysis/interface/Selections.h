#ifndef Selections_H
#define Selections_H

#include "Workspace/EventSelectors/interface/EventSelectorFactory.h"

class Filter {
 public:
  Filter(const edm::ParameterSet& iConfig);
  Filter(std::string name, edm::ParameterSet& iConfig) : 
    name_(name),inverted_(false), selector_(0)
  {
    if (!iConfig.empty()){
      const std::string d("name");
      iConfig.addUntrackedParameter<std::string>(d,name);
      std::string componentName = iConfig.getParameter<std::string>("selector");
      selector_ = EventSelectorFactory::get()->create(componentName, iConfig);
    }
  }
 
  const std::string & name() {return name_;}
  bool accept(edm::Event& iEvent) {
    bool decision=false;
    if (selector_)
      decision=selector_->select(iEvent);
    else
      decision=true;
    if (inverted_) return !decision;
    else return decision;
  }
  void setInverted() {    inverted_=true; }

 private:
  std::string name_;
  bool inverted_;//too allow !filter
  SusyEventSelector * selector_;
};

//forward declaration for friendship
class Selections;

class Selection {
 public:
  typedef std::vector<Filter*>::iterator iterator;
  friend class Selections;

  Selection(std::string name, const edm::ParameterSet& iConfig) :
    name_(name), 
    ntuplize_(iConfig.getParameter<bool>("ntuplize")),
    makeFinalPlots_(iConfig.getParameter<bool>("makeFinalPlots")),
    makeCumulativePlots_(iConfig.getParameter<bool>("makeCumulativePlots")),
    makeAllButOnePlots_(iConfig.getParameter<bool>("makeAllButOnePlots")),
    nSeen_(0),
    makeSummaryTable_(iConfig.getParameter<bool>("makeSummaryTable"))
  {
  }

  const std::string & name() {return name_;}
  iterator begin() { return filters_.begin();}
  iterator end() { return filters_.end();}

  std::map<std::string, bool> accept(edm::Event& iEvent){
    nSeen_++;
    std::map<std::string, bool> ret;
    bool global=true;
    for (iterator filter=begin(); filter!=end();++filter){
      const std::string & fName=(*filter)->name();
      Count & count=counts_[fName];
      count.nSeen_++;
      bool decision=(*filter)->accept(iEvent);
      ret[fName]=decision;
      if (decision) count.nPass_++;
      global=global && decision;
      if (global) count.nCumulative_++;
    }
    return ret;
  }

  //print to LogVerbatim("Selections|<name()>")
  void print(){
    if (!makeSummaryTable_) return;
    const std::string category ="Selections|"+name();
    edm::LogInfo(category)<<"   Summary table for selection: "<<name()<<" with: "<<nSeen_<<" events run.";
    std::cout<<"   Summary table for selection: "<<name()<<" with: "<<nSeen_<<" events run."<<std::endl;
    if (nSeen_==0) return;
    for (iterator filter=begin(); filter!=end();++filter){
      const std::string & fName=(*filter)->name();
      const Count & count=counts_[fName];
      edm::LogVerbatim(category)<<fName<<" has: "<<count.nPass_<<" passed events. "<<(float)(count.nPass_/count.nSeen_)*100.<<" [%]";
      std::cout<<fName<<" has: "<<count.nPass_<<" passed events. "<<(float)(count.nPass_/count.nSeen_)*100.<<" [%]"<<std::endl;
      std::cout<<fName<<" has: "<<count.nCumulative_<<" cumulative passed events. "<<(float)(count.nCumulative_/count.nSeen_)*100.<<" [%]"<<std::endl;
    }
    edm::LogVerbatim(category)<<"-------------------------------------";
    std::cout<<"-------------------------------------"<<std::endl;
  };


  bool ntuplize() {return ntuplize_;}
  bool makeFinalPlots() { return makeFinalPlots_;}
  bool makeCumulativePlots() { return makeCumulativePlots_;}
  bool makeAllButOnePlots() { return makeAllButOnePlots_;}
  bool makeSummaryTable() { return makeSummaryTable_;}

 private:
  std::string name_;
  std::vector<Filter*> filters_;
  //some options
  bool ntuplize_;
  bool makeFinalPlots_;
  bool makeCumulativePlots_;
  bool makeAllButOnePlots_;

  uint nSeen_;
  struct Count{
    uint nPass_;
    uint nSeen_;
    uint nCumulative_;
  };
  std::map<std::string, Count> counts_;
  bool makeSummaryTable_;
};

class Selections {
 public:
  typedef std::vector<Selection>::iterator iterator;

  Selections(const edm::ParameterSet& iConfig) : 
    filtersPSet_(iConfig.getParameter<edm::ParameterSet>("filters")),
    selectionPSet_(iConfig.getParameter<edm::ParameterSet>("selections"))
  {
    //FIXME. what about nested filters
    //make all configured filters
    std::vector<std::string> filterNames;
    uint nF=filtersPSet_.getParameterSetNames(filterNames);
    for (uint iF=0;iF!=nF;iF++){
      edm::ParameterSet pset = filtersPSet_.getParameter<edm::ParameterSet>(filterNames[iF]);
      filters_.insert(std::make_pair(filterNames[iF],Filter(filterNames[iF],pset)));
    }

    //parse all configured selections
    std::vector<std::string> selectionNames;
    std::map<std::string, std::vector<std::string> > selectionFilters;
    uint nS=selectionPSet_.getParameterSetNames(selectionNames);
    for (uint iS=0;iS!=nS;iS++){
      edm::ParameterSet pset=selectionPSet_.getParameter<edm::ParameterSet>(selectionNames[iS]);
      selections_.push_back(Selection(selectionNames[iS],pset));
      //      selections_.insert(std::make_pair(selectionNames[iS],Selection(selectionNames[iS],pset)));
      //keep track of list of filters for this selection for further dependency resolution
      selectionFilters[selectionNames[iS]]=pset.getParameter<std::vector<std::string> >("filterOrder");
    }


    //watch out of recursive dependency
    uint nestedDepth=0; //FIXME not taken care of

    //resolving dependencies
    for (std::map<std::string, std::vector<std::string> >::iterator sIt= selectionFilters.begin();sIt!=selectionFilters.end();++sIt)
      {
	//parse the vector of filterNames
	for (std::vector<std::string>::iterator fOrS=sIt->second.begin();fOrS!=sIt->second.end();++fOrS)
	  {
	    if (filters_.find(*fOrS)==filters_.end())
	      {
		//not a know filter names uncountered
		// look for a selection name
		std::map<std::string, std::vector<std::string> >::iterator s=selectionFilters.find(*fOrS);
		if (s==selectionFilters.end()){
		  //error. 
		  edm::LogError("SelectionHelper")<<"unresolved filter/selection name: "<<*fOrS;
		}
		else{
		  //remove the occurence
		  std::vector<std::string>::iterator newLoc=sIt->second.erase(fOrS);
		  //insert the list of filters corresponding to this selection in there
		  sIt->second.insert(newLoc,s->second.begin(),s->second.end());
		  //decrement selection iterator to come back to it
		  sIt--;
		  break;
		}
	      }
	      
	  }//loop over the string in "filterOrder"
      }//loop over all defined Selection

    //finally, configure the Selections
    //loop the selections instanciated
    //    for (std::map<std::string, Selection>::iterator sIt=selections_.begin();sIt!=selections_.end();++sIt)
    //      const std::string & sName=sIt->first;
    //Selection & selection =sIt->second;
    for (std::vector<Selection>::iterator sIt=selections_.begin();sIt!=selections_.end();++sIt){
      const std::string & sName=sIt->name();    
      Selection & selection =*sIt;

      //parse the vector of filterNames
      std::vector<std::string> & listOfFilters=selectionFilters[sName];
      for (std::vector<std::string>::iterator fIt=listOfFilters.begin();fIt!=listOfFilters.end();++fIt)
	{
	  std::map<std::string, Filter>::iterator filterInstance=filters_.find(*fIt);
	  if (filterInstance==filters_.end()){
	    //error
	    edm::LogError("Selections")<<"cannot resolve: "<<*fIt;
	  }
	  else{
	    //actually increment the filter
	    selection.filters_.push_back(&filterInstance->second);
	  }
	}
    }


  }

  iterator begin() {return selections_.begin(); }
  iterator end() { return selections_.end();}

  //print each selection 
  void print(){ for (std::vector<Selection>::iterator sIt=selections_.begin();sIt!=selections_.end();++sIt) sIt->print();}
    
 private:
  edm::ParameterSet filtersPSet_;
  std::map<std::string, Filter> filters_;

  edm::ParameterSet selectionPSet_;
  //  std::map<std::string, Selection> selections_;
  std::vector<Selection> selections_;
};


/*
 * Description: 
 * class to allows user defined collections from the event (pT ordered jets...)
*/
/*
class Retriever {
 public:
  typedef std::vector<CandidateRefToBase> Objects;
  typedef std::map<std::string, Objects> AvailableObjects;

  Retriever(const edm::ParameterSet& iConfig);
  Objects & get(std::string id){
  //should first retrieve it if not there yet
  return objects_[id]; }
  void clear();

 private:
  AvailableObjects objects_;
};
*/

/*
 * Description:
 * placeholder for common plotting tools
 *
 */
class Plotter{
 public:
  Plotter(const edm::ParameterSet& iConfig){};

  //set the directory on which to put plot next
  void setDir(std::string dir){};
  
  void fill(std::string selectionName, edm::Event& iEvent){};

  //writeout the root file
  void write(){};
  
 private:
  
};


/*
 * Description:
 * placeholder for common ntuplizer tools
 *
 */
class Ntupler{
 public:
  Ntupler(const edm::ParameterSet& iConfig){};

  void registerleaves(edm::EDProducer * producer, std::string selectionName){
    //produces<double>("leafName");
  }
  void fill(std::string selectionName, const edm::Event& iEvent){};
 private:
};

#endif
