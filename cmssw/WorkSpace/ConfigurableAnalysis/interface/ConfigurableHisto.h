#ifndef ConfigurableAnalysis_ConfigurableHisto_H
#define ConfigurableAnalysis_ConfigurableHisto_H

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/UtilAlgos/interface/TFileDirectory.h"

#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "THStack.h"

#include "Workspace/ConfigurableAnalysis/interface/VariableHelper.h"
#include "Workspace/ConfigurableAnalysis/interface/CachingVariable.h"


class ConfigurableAxis {
 public:
  ConfigurableAxis(){}
  ConfigurableAxis(edm::ParameterSet par) :
    nBin_(0),Min_(0),Max_(0),Label_(""){
    Label_=par.getParameter<std::string>("Label");

    if (par.exists("nBins")){
      nBin_=par.getParameter<uint>("nBins");
      Min_=par.getParameter<double>("Min");
      Max_=par.getParameter<double>("Max");
    }else{
      if (par.exists("vBins"))
	vBins_=par.getParameter<std::vector<double> >("vBins");
      else{
	Min_=par.getParameter<double>("Min");
	Max_=par.getParameter<double>("Max");
      }
    }
  }
  bool variableSize(){return vBins_.size()!=0;}
  uint nBin(){if (variableSize()) return vBins_.size()-1;else return nBin_;}
  double Min(){if (variableSize()) return vBins_.front(); else return Min_;}
  double Max(){if (variableSize()) return vBins_.back(); else return Max_;}
  const std::string & Label(){ return Label_;}
  const double * xBins(){ if (vBins_.size()!=0) return &(vBins_.front()); else return 0;}

 private:
  std::vector<double> vBins_;
  uint nBin_;
  double Min_;
  double Max_;
  std::string Label_;
};


class ConfigurableHisto {
 public:
  enum HType { h1 ,h2, prof };
  ConfigurableHisto(HType t, std::string name, edm::ParameterSet & iConfig) :
    type_(t),h_(0),name_(name), conf_(iConfig),x_(0),y_(0),z_(0),w_(0){}

  virtual ~ConfigurableHisto(){}

  virtual ConfigurableHisto * clone() const { return new ConfigurableHisto(*this);}

  virtual void book(TFileDirectory* dir){
    std::string title=conf_.getParameter<std::string>("title");
    edm::ParameterSet xAxisPSet=conf_.getParameter<edm::ParameterSet>("xAxis");
    ConfigurableAxis xAxis(xAxisPSet);
    //    x_=edm::Service<VariableHelperService>()->variable(xAxisPSet.getParameter<std::string>("var"));
    x_=VariableHelperInstance::get().variable(xAxisPSet.getParameter<std::string>("var"));

    std::string yLabel="";    
    bool yVBin=false;
    ConfigurableAxis yAxis;
    if (conf_.exists("yAxis")){
      edm::ParameterSet yAxisPSet=conf_.getParameter<edm::ParameterSet>("yAxis");
      yAxis=ConfigurableAxis(yAxisPSet);
      yLabel=yAxis.Label();
      //at least TH2 or TProfile
      if (yAxisPSet.exists("var"))
	//y_=edm::Service<VariableHelperService>()->variable(yAxisPSet.getParameter<std::string>("var"));
	y_=VariableHelperInstance::get().variable(yAxisPSet.getParameter<std::string>("var"));
      yVBin=yAxis.variableSize();
    }
    
    if (conf_.exists("zAxis")){
      throw;
    }

    bool xVBin=xAxis.variableSize();


    if (type()==h1){
      //make TH1F
      if (xVBin)
	h_=dir->make<TH1F>(name_.c_str(),title.c_str(),
			   xAxis.nBin(), xAxis.xBins());
      else
	h_=dir->make<TH1F>(name_.c_str(),title.c_str(),
			   xAxis.nBin(),xAxis.Min(),xAxis.Max());
    }
    else if (type()==h2){
      //make TH2F
      if (xVBin){
	if (yVBin)
	  h_=dir->make<TH2F>(name_.c_str(),title.c_str(),
			     xAxis.nBin(),xAxis.xBins(),
			     yAxis.nBin(),yAxis.xBins());
	else
	  h_=dir->make<TH2F>(name_.c_str(),title.c_str(),
			     xAxis.nBin(),xAxis.xBins(),
			     yAxis.nBin(),yAxis.Min(),yAxis.Max());
      }else{
	if (yVBin)
	  h_=dir->make<TH2F>(name_.c_str(),title.c_str(),
			     xAxis.nBin(),xAxis.Min(),xAxis.Max(),
			     yAxis.nBin(),yAxis.xBins());
	else
	  h_=dir->make<TH2F>(name_.c_str(),title.c_str(),
			     xAxis.nBin(),xAxis.Min(),xAxis.Max(),
			     yAxis.nBin(),yAxis.Min(),yAxis.Max());
      }
    }
    else if (type()==prof){
      //make TProfile
      std::string pFopt="";
      if (conf_.exists("Option"))
	pFopt=conf_.getParameter<std::string>("Option");
      
      if (xVBin)
      h_=dir->make<TProfile>(name_.c_str(),title.c_str(),
			     xAxis.nBin(), xAxis.xBins(),
			     yAxis.Min(),yAxis.Max(),
			     pFopt.c_str());
      else
	h_=dir->make<TProfile>(name_.c_str(),title.c_str(),
			       xAxis.nBin(),xAxis.Min(),xAxis.Max(),
			       yAxis.Min(),yAxis.Max(),
			       pFopt.c_str());
      
    }
    else {
      edm::LogError("ConfigurableHisto")<<"cannot book: "<<name_<<"\n"<<conf_.dump();
      throw;
    }

    //cosmetics
    h_->GetXaxis()->SetTitle(xAxis.Label().c_str());
    h_->SetYTitle(yLabel.c_str());
    
    if (conf_.exists("weight"))
      {
	w_=VariableHelperInstance::get().variable(conf_.getParameter<std::string>("weight"));
      }
  }
  
  virtual void fill(const edm::Event &) {
    if (!h_)      return;

    double weight=1.0;
    if (w_){
      if (w_->compute())
	weight=(*w_)();
      else{
	edm::LogInfo("ConfigurableHisto")<<"could not compute the weight for: "<<name_
					 <<" with config:\n"<<conf_.dump()
					 <<" default to 1.0";
      }
    }
    
    switch(type_){
    case h1:
      if (!h_) throw;
      if (x_->compute()) h_->Fill((*x_)(),weight);
      else{
	edm::LogInfo("ConfigurableHisto")<<"could not fill: "<<name_
					 <<" with config:\n"<<conf_.dump();
      }
      break;
    case prof:
      TProfile * p=dynamic_cast<TProfile*>(h_);
      if (!p) throw;
      if (x_->compute() && y_->compute()) p->Fill((*x_)(),(*y_)(),weight);
      else{
	edm::LogInfo("ConfigurableHisto")<<"could not fill: "<<name_
					 <<" with config:\n"<<conf_.dump();
      }
      break;
    case h2:
      TH2 * h2=dynamic_cast<TH2*>(h_);
      if (!h2) throw;
      if (x_->compute() && y_->compute()) h2->Fill((*x_)(),(*y_)(),weight);
      else{
	edm::LogInfo("ConfigurableHisto")<<"could not fill: "<<name_
                                         <<" with config:\n"<<conf_.dump();
      }
      break;
    }
  }
  const HType & type() { return type_;}  

  void complete() {}
  TH1 * h() {return h_;}

 protected:
  ConfigurableHisto(const ConfigurableHisto & master){
    type_=master.type_;
    h_=0; //no histogram attached in copy constructor
    name_=master.name_;
    conf_=master.conf_;
    x_=master.x_;
    y_=master.y_;
    z_=master.z_;
    w_=master.w_;
  }
  HType type_;
  TH1 * h_;
  std::string name_;
  edm::ParameterSet conf_;
  
  const CachingVariable * x_;
  const CachingVariable * y_;
  const CachingVariable * z_;
  const CachingVariable * w_;
};

class SplittingConfigurableHisto : public ConfigurableHisto {
 public:
  SplittingConfigurableHisto(HType t, std::string name, edm::ParameterSet & pset) :
    ConfigurableHisto(t,name,pset) , splitter_(0) {
    //---    std::cout<<"getting IN the constructor of SplittingConfigurableHisto for "<< name_<<std::endl;
    std::string title=pset.getParameter<std::string>("title");

    //allow for many splitters ...
    if (pset.exists("splitters")){
      //you want more than one splitter
      std::vector<std::string> splitters = pset.getParameter<std::vector<std::string> >("splitters");
      //---      std::cout<<splitters.size()<<" splitters"<<std::endl;
      for (uint s=0;s!=splitters.size();++s){
	//---	std::cout<<"trying with "<<splitters[s]<<std::endl;
	const CachingVariable * v=VariableHelperInstance::get().variable(splitters[s]);
	const Splitter * splitter = dynamic_cast<const Splitter*>(v);
	if (!splitter){
	  edm::LogError("SplittingConfigurableHisto")<<"for: "<<name_<<" the splitting variable: "<<splitters[s]<<" is not a Splitter";
	  continue;
	}
	//---	std::cout<<"ok for:"<<name_<<" and "<<splitters[s]<<std::endl;

	//insert in the map
	std::vector<ConfigurableHisto*> & insertedHisto=subHistoMap_[splitter];
	//now configure the histograms
	uint mSlots=splitter->maxSlots();
	for (uint i=0;i!=mSlots;++i){
	  //---	  std::cout<<" slot: "<<i<<std::endl;
	  const std::string & slabel=splitter->shortLabel(i);
	  const std::string & label=splitter->label(i);
	  edm::ParameterSet mPset=pset;
	  edm::Entry e("string",title+" for "+label,true);
	  mPset.insert(true,"title",e);
	  //---	  std::cout<<"new for: "<<name_<<" and "<<slabel<<std::endl;
	  insertedHisto.push_back(new ConfigurableHisto(t,name+slabel ,mPset));
	}//loop on slots
      }//loop on splitters

    }//if splitters exists
    else{
      //single splitter
      //get the splitting variable
      const CachingVariable * v=VariableHelperInstance::get().variable(pset.getParameter<std::string>("splitter"));
      //---    std::cout<<"trying to cast"<<std::endl;
      splitter_ = dynamic_cast<const Splitter*>(v);
      if (!splitter_){
	//      std::cout<<"for: "<<name_<<" the splitting variable: "<<v->name()<<" is not a Splitter"<<std::endl;
	edm::LogError("SplittingConfigurableHisto")<<"for: "<<name_<<" the splitting variable: "<<v->name()<<" is not a Splitter";
      }
      else{
	//---      std::cout<<"configure slots"<<std::endl;
	//configure the splitted plots
	uint mSlots=splitter_->maxSlots();
	for (uint i=0;i!=mSlots;i++){
	  const std::string & slabel=splitter_->shortLabel(i);
	  const std::string & label=splitter_->label(i);
	  edm::ParameterSet mPset=pset;
	  edm::Entry e("string",title+" for "+label,true);
	  mPset.insert(true,"title",e);
	  subHistos_.push_back(new ConfigurableHisto(t,name+slabel, mPset));
	}
      }
    }
      //---    std::cout<<"getting OUT the constructor of SplittingConfigurableHisto"<<std::endl;
  }//end of ctor
    
  void book(TFileDirectory* dir){
    //book the base histogram
    ConfigurableHisto::book(dir);

    if (subHistoMap_.size()!=0){
      SubHistoMap::iterator i=subHistoMap_.begin();
      SubHistoMap::iterator i_end=subHistoMap_.end();
      for (;i!=i_end;++i){for (uint h=0;h!=i->second.size();++h){ 
	  i->second[h]->book(dir);}
	//book the THStack
	std::string sName= name_+"_"+i->first->name();
	std::string sTitle="Stack histogram of "+name_+" for splitter "+i->first->name();
	subHistoStacks_[i->first]= dir->make<THStack>(sName.c_str(),sTitle.c_str());
      }
    }else{
      for (uint h=0;h!=subHistos_.size();h++){subHistos_[h]->book(dir);}
      //book a THStack
      std::string sName= name_+"_"+splitter_->name();
      std::string sTitle="Stack histogram of "+name_+" for splitter "+splitter_->name();
      stack_ = dir->make<THStack>(sName.c_str(),sTitle.c_str());
    }
    
  }

  ConfigurableHisto * clone() const {    return new SplittingConfigurableHisto(*this);  }

  void fill(const edm::Event & e) {
    //fill the base histogram
    ConfigurableHisto::fill(e);
    
    if (subHistoMap_.size()!=0){
      SubHistoMap::iterator i=subHistoMap_.begin();
      SubHistoMap::iterator i_end=subHistoMap_.end();
      for (;i!=i_end;++i){
	const Splitter * splitter=i->first;
	if (!splitter) continue;
	if (!splitter->compute()) continue;
	uint slot=(uint)  (*splitter)();
	//fill in the proper slot
	i->second[slot]->fill(e);
      }
    }
    else{
      //fill the component histograms
      //    std::cout<<"splitter is: "<<splitter_<<std::endl;
      if (!splitter_) return;
      if (!splitter_->compute()){
	//      std::cout<<"variable does not compute."<<std::endl;
	return;}
      uint slot=(uint) (*splitter_)();
      //    std::cout<<"filling in slot: "<<slot<<std::endl;
      subHistos_[slot]->fill(e);
    }
  }

  void complete(){
    if (subHistoMap_.size()!=0){
      //fill up the stacks
      SubHistoMap::iterator i=subHistoMap_.begin();
      SubHistoMap::iterator i_end=subHistoMap_.end();
      for (;i!=i_end;++i){
	for (uint h=0;h!=i->second.size();h++){
	  //	  if (i->second[h]->h()->Integral==0) continue;// do not add empty histograms. NO, because it will be tough to merge two THStack
	  subHistoStacks_[i->first]->Add(i->second[h]->h(), i->first->label(h).c_str());
	}
      }

    }else{
      //fill up the only stack
      for (uint i=0;i!=subHistos_.size();i++){	stack_->Add(subHistos_[i]->h(), splitter_->label(i).c_str());      }
    }


  }
 private:
  SplittingConfigurableHisto(const SplittingConfigurableHisto & master) : ConfigurableHisto(master){
    splitter_ = master.splitter_;
    if (master.subHistoMap_.size()!=0){
      SubHistoMap::const_iterator i=master.subHistoMap_.begin();
      SubHistoMap::const_iterator i_end=master.subHistoMap_.end();
      for (;i!=i_end;++i){
	const std::vector<ConfigurableHisto*> & masterHistos=i->second;
	std::vector<ConfigurableHisto*> & clonedHistos=subHistoMap_[i->first];
	for (uint i=0;i!=masterHistos.size();i++){clonedHistos.push_back(masterHistos[i]->clone());}
      }
    }else{
      for (uint i=0;i!=master.subHistos_.size();i++){subHistos_.push_back(master.subHistos_[i]->clone());}
    }
  }

  typedef std::map<const Splitter *, std::vector<ConfigurableHisto* > > SubHistoMap;
  typedef std::map<const Splitter *, THStack *> SubHistoStacks;
  SubHistoStacks subHistoStacks_;
  SubHistoMap subHistoMap_;
  
  const Splitter * splitter_;
  std::vector<ConfigurableHisto* > subHistos_;
  // do you want to have a stack already made from the splitted histos?
  THStack * stack_;
};

#endif
