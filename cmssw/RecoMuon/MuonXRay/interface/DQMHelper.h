#ifndef RecoMuon_MuonXRay_DQMHelper
#define RecoMuon_MuonXRay_DQMHelper

#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"


template <class H>
H* MEA1(MonitorElement* me ){
  //  H * h =dynamic_cast<H*>(&(**((MonitorElementRootH1 *)me)));
  H * h = dynamic_cast<H*>(me->getTH1());
  if (!h) {edm::LogError("MonitorElementAdaptor")<<"in MEA, MonitorElement "<<me->getName()<<" does not cast to H1 types.";}
  return h;
}
template <class H>
H* MEA2(MonitorElement* me ){
  //  H * h =dynamic_cast<H*>(&(**((MonitorElementRootH2 *)me)));
  H * h = dynamic_cast<H*>(me->getTH1());
  if (!h) {edm::LogError("MonitorElementAdaptor")<<"in MEA, MonitorElement "<<me->getName()<<" does not cast to H2 types.";}
  return h;
}


#include "RecoMuon/MuonXRay/interface/ConfigurableHisto.h"




class DQMHelper {
public:
  DQMHelper(edm::ParameterSet  par):
    dbe_(edm::Service<DQMStore>().operator->()){
    //    dump_=dbe_->book1D("dump","value for non-register histograms",1,0,1);

    if (par.exists("directory")){
      std::string dir=par.getParameter<std::string>("directory");
      //replace @module_label by module_label in directory structure
      //can this grab the module label of where it si located ... ?
      //    std::string module_label=par.getParameter<std::string>("@module_label");
      //    edm::LogWarning("DQMHelper")<<"find module_label="<<module_label;
      //and set it
      dbe_->setCurrentFolder(dir);
    }

    edm::ParameterSet H1s = par.getParameter<edm::ParameterSet>("H1s");
    std::vector<std::string> pnames=H1s.getParameterNames();
    for (uint ip=0;ip!=pnames.size();++ip){
      edm::ParameterSet ph1=H1s.getParameter<edm::ParameterSet>(pnames[ip]);
      ph1.addParameter<std::string>("name",pnames[ip]);
      ConfigurableH1 h1(ph1,dbe_);
      H[h1.name()]=h1.me();
    }

    edm::ParameterSet H2s = par.getParameter<edm::ParameterSet>("H2s");
    pnames=H2s.getParameterNames();
    for (uint ip=0;ip!=pnames.size();++ip){
      edm::ParameterSet ph2=H2s.getParameter<edm::ParameterSet>(pnames[ip]);
      ph2.addParameter<std::string>("name",pnames[ip]);
      ConfigurableH2 h2(ph2,dbe_);
      H[h2.name()]=h2.me();
    }

    /*
       edm::ParameterSet stackH1s = par.getParameter<edm::ParameterSet>("sH1s");
       for (uint ip=0;ip!=stackH1s.size();++ip){
       ConfigurableStackH1 sh1(stackH1s[ip]);
    */
  }

  DQMHelper():
    dbe_(edm::Service<DQMStore>().operator->()){
    //dump_=dbe_->book1D("dump","value for non-register histograms",1,0,1);
  }
  DQMStore* dbe() {return dbe_;}
  ~DQMHelper(){
    //remove the dump histogram.
    //    if (dump_){      dbe_->removeElement("dump");}
  }
private:
  //DQM interface
  DQMStore* dbe_;
  typedef std::vector<MonitorElement* > byCategory;
  typedef std::map<std::string, byCategory > histosByCategory;
  typedef std::map<std::string, MonitorElement* > histos;
public:
  histos H;
  //  MonitorElement * dump_;
  template <class H1> H1* h1(std::string name) {  return MEA1<H1>(e(name));}
  template <class H2> H2* h2(std::string name) {  return MEA2<H2>(e(name));}
  
  MonitorElement * e(std::string name){
    histos::iterator f=H.find(name);
    if (f!=H.end()){  
      return f->second;
    }else{
      edm::LogError("DQMHelper")<<name<<" is not a registerd histogram. send to dump. or throw";
      //--a bit harch don't you think?      
      throw;
      //      return dump_;
    }
  }
  histosByCategory HBC;

  TH1* book1D(uint bin, std::string identifier,
	      std::string name, std::string title,
	      int nB, double xl, double xh,  
	      std::string Xl ="", std::string Yl =""){
    MonitorElement * me=dbe_->book1D(name, title, nB, xl, xh);
    HBC[identifier][bin]=me;
    TH1* h = MEA1<TH1>(me);
    h->SetXTitle(Xl.c_str());
    h->SetYTitle(Yl.c_str());
    return h;}

  TH1* book1D(std::string name, std::string title,
	      int nB, double xl, double xh,  
	      std::string Xl ="", std::string Yl =""){
    MonitorElement * me=dbe_->book1D(name, title, nB, xl, xh);
    H[name]=me;
    TH1* h = MEA1<TH1>(me);
    h->SetXTitle(Xl.c_str());
    h->SetYTitle(Yl.c_str());
    return h;}

  TH2* book2D(std::string name, std::string title,
	      int nBx, double xl, double xh, 
	      int nBy, double yl, double yh, 
	      std::string Xl="", std::string Yl=""){
    MonitorElement * me=dbe_->book2D(name, title, nBx, xl, xh, nBy, yl, yh);
    H[name]=me;
    TH2* h = MEA2<TH2>(me);
    h->SetXTitle(Xl.c_str());
    h->SetYTitle(Yl.c_str());
    return h;}

};


#endif
