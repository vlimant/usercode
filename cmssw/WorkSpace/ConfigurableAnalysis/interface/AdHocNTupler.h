#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerAlgorithm.h"

using namespace std;

class AdHocNTupler : public NTupler {
 public:
  AdHocNTupler (const edm::ParameterSet& iConfig){
    edm::ParameterSet adHocPSet = iConfig.getParameter<edm::ParameterSet>("AdHocNPSet");

    if (adHocPSet.exists("useTFileService"))
      useTFileService_=adHocPSet.getParameter<bool>("useTFileService");         
    else
      useTFileService_=iConfig.getParameter<bool>("useTFileService");

    if (useTFileService_){
      if (adHocPSet.exists("treeName")){
	treeName_=adHocPSet.getParameter<std::string>("treeName");
	ownTheTree_=true;
      }
      else{
	treeName_=iConfig.getParameter<std::string>("treeName");
	ownTheTree_=false;
      }
    }
    
    trigger_prescalevalue = new std::vector<float>;
    trigger_name = new std::vector<std::string>;
    trigger_decision = new std::vector<float>;
    trigger_lastfiltername = new std::vector<std::string>;
    triggerobject_pt = new std::vector<std::vector<float> >;
    triggerobject_px = new std::vector<std::vector<float> >;
    triggerobject_py = new std::vector<std::vector<float> >;
    triggerobject_pz = new std::vector<std::vector<float> >;
    triggerobject_et = new std::vector<std::vector<float> >;
    triggerobject_energy = new std::vector<std::vector<float> >;
    triggerobject_phi = new std::vector<std::vector<float> >;
    triggerobject_eta = new std::vector<std::vector<float> >;
    triggerobject_collectionname = new std::vector<std::vector<TString> >;
    L1trigger_bit = new std::vector<float>;
    L1trigger_techTrigger = new std::vector<float>;
    L1trigger_prescalevalue = new std::vector<float>;
    L1trigger_name = new std::vector<std::string>;
    L1trigger_alias = new std::vector<std::string>;
    L1trigger_decision = new std::vector<float>;
    L1trigger_decision_nomask = new std::vector<float>;
  }

  ~AdHocNTupler(){
    delete trigger_prescalevalue;
    delete trigger_name;
    delete trigger_decision;
    delete trigger_lastfiltername;
    delete triggerobject_pt;
    delete triggerobject_px;
    delete triggerobject_py;
    delete triggerobject_pz;
    delete triggerobject_et;
    delete triggerobject_energy;
    delete triggerobject_phi;
    delete triggerobject_eta;
    delete triggerobject_collectionname;
    delete L1trigger_bit;
    delete L1trigger_techTrigger;
    delete L1trigger_prescalevalue;
    delete L1trigger_name;
    delete L1trigger_alias;
    delete L1trigger_decision;
    delete L1trigger_decision_nomask;

  }

  uint registerleaves(edm::ProducerBase * producer){
    uint nLeaves=0;
    if (useTFileService_){
      edm::Service<TFileService> fs;      
      if (ownTheTree_){
	ownTheTree_=true;
	tree_=fs->make<TTree>(treeName_.c_str(),"StringBasedNTupler tree");
      }else{
	TObject * object = fs->file().Get(treeName_.c_str());
	if (!object){
	  ownTheTree_=true;
	  tree_=fs->make<TTree>(treeName_.c_str(),"StringBasedNTupler tree");
	}
	tree_=dynamic_cast<TTree*>(object);
	if (!tree_){
	  ownTheTree_=true;
	  tree_=fs->make<TTree>(treeName_.c_str(),"StringBasedNTupler tree");
	}
      }
      
      //register the leaves by hand
      tree_->Branch("trigger_prescalevalue",&trigger_prescalevalue);
      tree_->Branch("trigger_name",&trigger_name);
      tree_->Branch("trigger_decision",&trigger_decision);
      tree_->Branch("trigger_lastfiltername",&trigger_lastfiltername);
      tree_->Branch("triggerobject_pt",&triggerobject_pt);
      tree_->Branch("triggerobject_px",&triggerobject_px);
      tree_->Branch("triggerobject_py",&triggerobject_py);
      tree_->Branch("triggerobject_pz",&triggerobject_pz);
      tree_->Branch("triggerobject_et",&triggerobject_et);
      tree_->Branch("triggerobject_energy",&triggerobject_energy);
      tree_->Branch("triggerobject_phi",&triggerobject_phi);
      tree_->Branch("triggerobject_eta",&triggerobject_eta);
      tree_->Branch("triggerobject_collectionname",&triggerobject_collectionname);
      tree_->Branch("L1trigger_bit",&L1trigger_bit);
      tree_->Branch("L1trigger_techTrigger",&L1trigger_techTrigger);
      tree_->Branch("L1trigger_prescalevalue",&L1trigger_prescalevalue);
      tree_->Branch("L1trigger_name",&L1trigger_name);
      tree_->Branch("L1trigger_alias",&L1trigger_alias);
      tree_->Branch("L1trigger_decision",&L1trigger_decision);
      tree_->Branch("L1trigger_decision_nomask",&L1trigger_decision_nomask);

    }

    else{
      //EDM COMPLIANT PART
      //      producer->produce<ACertainCollection>(ACertainInstanceName);
    }


    return nLeaves;
  }

  void fill(edm::Event& iEvent){
    //open the collection that you want
    //retrieve the objects
    //fill the variable for tree filling 


    edm::Handle< pat::TriggerEvent > triggerevent;
    iEvent.getByLabel("patTriggerEvent",triggerevent);  
    //    std::cout<<"The trigger HLT table is"<<triggerevent->nameHltTable()<<std::endl;


    edm::Handle< std::vector<pat::TriggerPath> > triggerpaths;
    iEvent.getByLabel("patTrigger",triggerpaths);  
    for( std::vector<pat::TriggerPath>::const_iterator tp=triggerpaths->begin(); tp!=triggerpaths->end(); ++tp ){
      double prescalevalue = tp->prescale(); 
      std::string name = tp->name();
      float decision = tp->wasAccept();
      (*trigger_prescalevalue).push_back(prescalevalue);
      (*trigger_name).push_back(name);
      (*trigger_decision).push_back(decision);

      std::vector<TString> collection_names;
      std::vector<float> pt_vector;
      std::vector<float> px_vector;
      std::vector<float> py_vector;
      std::vector<float> pz_vector;
      std::vector<float> et_vector;
      std::vector<float> energy_vector;
      std::vector<float> phi_vector;
      std::vector<float> eta_vector;

      /*
      cout<<""<<endl;
      cout<<"Trigger names is: "<<name<<endl;
      cout<<"Trigger decision is: "<<decision<<endl;
      */

      edm::RefVector< pat::TriggerFilterCollection> toFilt = triggerevent->pathFilters(name);
      edm::RefVector<pat::TriggerFilterCollection>::const_iterator toFilt_it=toFilt.end();
       if(toFilt.size()>0){
	const pat::TriggerFilter *triggerfilter = (*(--toFilt_it)).get();
	(*trigger_lastfiltername).push_back(triggerfilter->label());
	//cout<<"The trigger filter is: "<<triggerfilter->label()<<endl;
	edm::RefVector< pat::TriggerObjectCollection > tocoll = triggerevent->filterObjects(triggerfilter->label());
      
	for( edm::RefVector<pat::TriggerObjectCollection>::const_iterator to=tocoll.begin(); to!=tocoll.end(); ++to ){
	  const pat::TriggerObject *triggerObject = (*to).get();
	  double pt = triggerObject->pt(); 
	  double px = triggerObject->px(); 
	  double py = triggerObject->py(); 
	  double pz = triggerObject->pz(); 
          double et = triggerObject->et();
	  double energy = triggerObject->energy(); 
          double phi = triggerObject->phi();
          double eta = triggerObject->eta();
	  TString collname(triggerObject->collection());
	  //cout<<"The trigger collname is: "<<collname<<endl;
	  //cout<<"The trigger objectpt is: "<<pt<<endl;
	  collection_names.push_back(collname);
	  pt_vector.push_back(pt);
	  px_vector.push_back(px);
	  py_vector.push_back(py);
	  pz_vector.push_back(pz);
          et_vector.push_back(et);
	  energy_vector.push_back(energy);
          phi_vector.push_back(phi);
          eta_vector.push_back(eta);
	}
	
	(*triggerobject_collectionname).push_back(collection_names);
	(*triggerobject_pt).push_back(pt_vector);
	(*triggerobject_px).push_back(px_vector);
	(*triggerobject_py).push_back(py_vector);
	(*triggerobject_pz).push_back(pz_vector);
        (*triggerobject_et).push_back(et_vector);
	(*triggerobject_energy).push_back(energy_vector);
        (*triggerobject_phi).push_back(phi_vector);
        (*triggerobject_eta).push_back(eta_vector);
       }//end of if statement requiring that reVector be greater than 0
       else{
	 (*trigger_lastfiltername).push_back("none");
	 (*triggerobject_collectionname).push_back(collection_names);
	 (*triggerobject_pt).push_back(pt_vector);
	 (*triggerobject_px).push_back(px_vector);
	 (*triggerobject_py).push_back(py_vector);
	 (*triggerobject_pz).push_back(pz_vector);
         (*triggerobject_et).push_back(et_vector);
	 (*triggerobject_energy).push_back(energy_vector);
         (*triggerobject_phi).push_back(phi_vector);
         (*triggerobject_eta).push_back(eta_vector);
       }
      collection_names.clear();
      pt_vector.clear();
      px_vector.clear();
      py_vector.clear();
      pz_vector.clear();
      et_vector.clear();
      energy_vector.clear();
      phi_vector.clear();
      eta_vector.clear();
    }
    


    edm::Handle< std::vector<pat::TriggerAlgorithm> > triggeralgos;
    iEvent.getByLabel("patTrigger",triggeralgos);
    for( std::vector<pat::TriggerAlgorithm>::const_iterator ta=triggeralgos->begin(); ta!=triggeralgos->end(); ++ta ){
      float prescalevalue = ta->prescale();
      std::string name = ta->name();
      std::string alias = ta->alias();
      float bit = ta->bit();
      float techTrig = ta->techTrigger();
      float decision = ta->decision();
      float decision_nomask = ta->decisionBeforeMask();
      (*L1trigger_prescalevalue).push_back(prescalevalue);
      (*L1trigger_name).push_back(name);
      (*L1trigger_alias).push_back(alias);
      (*L1trigger_bit).push_back(bit);
      (*L1trigger_techTrigger).push_back(techTrig);
      (*L1trigger_decision).push_back(decision);
      (*L1trigger_decision_nomask).push_back(decision_nomask);
    }



    //fill the tree    
    if (ownTheTree_){ tree_->Fill(); }
    (*trigger_prescalevalue).clear();
    (*trigger_name).clear();
    (*trigger_decision).clear();
    (*trigger_lastfiltername).clear();
    (*triggerobject_pt).clear();
    (*triggerobject_px).clear();
    (*triggerobject_py).clear();
    (*triggerobject_pz).clear();
    (*triggerobject_et).clear();
    (*triggerobject_energy).clear();
    (*triggerobject_phi).clear();
    (*triggerobject_eta).clear();
    (*triggerobject_collectionname).clear();
    (*L1trigger_bit).clear();
    (*L1trigger_techTrigger).clear();
    (*L1trigger_prescalevalue).clear();
    (*L1trigger_name).clear();
    (*L1trigger_alias).clear();
    (*L1trigger_decision).clear();
    (*L1trigger_decision_nomask).clear();

  }

  void callBack(){
    //clean up whatever memory was allocated
  }

 private:
  bool ownTheTree_;
  std::string treeName_;
  bool useTFileService_;

  std::vector<float> * trigger_prescalevalue;
  std::vector<std::string> * trigger_name;
  std::vector<float> * trigger_decision;
  std::vector<std::string> * trigger_lastfiltername;
  std::vector<std::vector<float> > * triggerobject_pt;
  std::vector<std::vector<float> > * triggerobject_px;
  std::vector<std::vector<float> > * triggerobject_py;
  std::vector<std::vector<float> > * triggerobject_pz;
  std::vector<std::vector<float> > * triggerobject_et;
  std::vector<std::vector<float> > * triggerobject_energy;
  std::vector<std::vector<float> > * triggerobject_phi;
  std::vector<std::vector<float> > * triggerobject_eta;
  std::vector<std::vector<TString> > * triggerobject_collectionname;
  std::vector<float> * L1trigger_bit;
  std::vector<float> * L1trigger_techTrigger;
  std::vector<float> * L1trigger_prescalevalue;
  std::vector<std::string> * L1trigger_name;
  std::vector<std::string> * L1trigger_alias;
  std::vector<float> * L1trigger_decision;
  std::vector<float> * L1trigger_decision_nomask;

};
