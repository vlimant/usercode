#include "DataFormats/PatCandidates/interface/Electron.h"

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
    
    e_leading_pt_ = new double; 

  }

  ~AdHocNTupler(){
    delete e_leading_pt_;
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
      tree_->Branch("e_leading_pt",e_leading_pt_,"e_leading_pt/d");


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

    edm::Handle< std::vector<pat::Electron> > electrons;
    iEvent.getByLabel("cleanPatElectrons",electrons);

    double highest_pt_e = 0;
    for( std::vector<pat::Electron>::const_iterator elec=electrons->begin(); elec!=electrons->end(); ++elec ){
      double pt = elec->pt();
      if(pt>highest_pt_e){highest_pt_e=pt;}
    }


    *e_leading_pt_= highest_pt_e;

    //fill the tree    
    if (ownTheTree_){ tree_->Fill(); }
  }

  void callBack(){
    //clean up whatever memory was allocated
  }

 private:
  bool ownTheTree_;
  std::string treeName_;
  bool useTFileService_;

  double * e_leading_pt_;

};
