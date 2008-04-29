// -*- C++ -*-
//
// Package:    MuonXRay
// Class:      MuonXRay
// 
/**\class MuonXRay MuonXRay.cc Workspace/MuonXRay/src/MuonXRay.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Finn Rebassoo
//         Created:  Tue Jan  8 17:56:28 CST 2008
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include <SimDataFormats/HepMCProduct/interface/HepMCProduct.h>
#include "HepMC/GenEvent.h"
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TVector3.h"
#include "TMath.h"

//
// class decleration
//

#include "RecoMuon/MuonXRay/interface/IDconverttoBinNum.h"
#include "RecoMuon/MuonXRay/interface/MotherSearch.h"
#include "RecoMuon/MuonXRay/interface/DQMHelper.h"

class MuonXRay : public edm::EDAnalyzer {
   public:
      explicit MuonXRay(const edm::ParameterSet&);
      ~MuonXRay();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
  std::string theCategory;
  DQMHelper h;
  IDconverttoBinNum wantMotherBin;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MuonXRay::MuonXRay(const edm::ParameterSet& iConfig):
  theCategory("MuonXRay"),
  wantMotherBin(iConfig.getParameter<edm::ParameterSet>("IDconverttoBinNum"))
{

  //set current folder
  h.dbe()->setCurrentFolder("MuonXRay/"+iConfig.getParameter<std::string>("@module_label"));
  //now that the directory is set
  h = DQMHelper(iConfig.getParameter<edm::ParameterSet>("DQMHelper"));

  wantMotherBin.splitH1ByCategory("h_mu_sim_pt",h,"%s coming from %s");
  wantMotherBin.splitH1ByCategory("h_mu_sim_eta",h,"%s coming from %s");
  wantMotherBin.splitH1ByCategory("h_mu_sim_Aeta",h,"%s coming from %s");
  wantMotherBin.splitH1ByCategory("h_sim_vertex_position",h,"%s coming from %s");           

  /*
  //  h.book1D("h_generated_pt_no_sim_vertex","pt of gen muons with no sim vertex",400,0,200);  
  //++  h.book1D("h_parent_id","Name of parent particle from tight association",16,-0.5,15.5,"parent");

  //++  h.book1D("h_mu_sim_pt","pt of sim muons",400,0,200, "muon p_{T}^{sim} [GeV]");
  //++  h.book1D("h_mu_sim_eta","#eta^{sim} of muon",50,-2.5,2.5, "muon #eta^{sim}");
  //++  h.book1D("h_sim_vertex_position","Sim vertex position",200,0,400,"muon production vertex distance [cm]");
  
  //++  h.HBC["h_mu_sim_pt_assoc_ID"].resize(wantMotherBin.size());
  //++  h.HBC["h_mu_sim_eta_assoc_ID"].resize(wantMotherBin.size());
  //++  h.HBC["h_sim_vertex_position_assoc_ID"].resize(wantMotherBin.size());
  std::string htitle,hname,binname;
  for(int iut = 0;iut<wantMotherBin.size();++iut)
  {
  binname = wantMotherBin.GetBinName(iut);
  
  //++      htitle = Form("p_{T}^{sim} of muon coming from %s",binname.c_str());
  //++      hname = Form("h_mu_sim_pt_assoc_ID_%i",iut);
  //++      h.book1D(iut,"h_mu_sim_pt_assoc_ID",hname,htitle,400,0,200,"muon p_{T}^{sim} [GeV]");
  
  //++      htitle = Form("|#eta|^{sim} of muon coming from %s",binname.c_str());
  //++      hname = Form("h_mu_sim_eta_assoc_ID_%i",iut);
  //++      h.book1D(iut,"h_mu_sim_eta_assoc_ID",hname,htitle,25,0,2.5,"muon |#eta|^{sim}");
  
  //++      htitle = Form("production vertex position distance from %s",binname.c_str());
  //++      hname = Form("h_sim_vertex_position_assoc_ID_%i",iut);
  //++      h.book1D(iut,"h_sim_vertex_position_assoc_ID",hname,htitle,200,0,400,"muon production vertex distance [cm]");
  }
  */

  //++  h.book1D("h_mu_sim_pt_tricky","pt of sim muons for tricky events",400,0,200, "muon p_{T}^{sim}","");
  //++  h.book1D("h_mu_sim_eta_tricky","|#eta|^{sim} of muon for tricky events",25,0,2.5, "muon |#eta|^{sim}","");

}


MuonXRay::~MuonXRay()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//



// ------------ method called to for each event  ------------
void
MuonXRay::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Handle<SimTrackContainer> SimTk;
   Handle<SimVertexContainer> SimVtx;
   Handle<HepMCProduct> hepmc;

   iEvent.getByType(hepmc);
   iEvent.getByLabel("g4SimHits",SimVtx);
   iEvent.getByLabel("g4SimHits",SimTk);
 
   std::map<int,int> thetrackIdDone;

   double eventWeight=1;
   /*
     uint nMuon=0;
     for(std::vector<SimTrack>::const_iterator isimtk = SimTk->begin();isimtk!=SimTk->end();++isimtk)
     {
     LogDebug(theCategory)<<"I am here 1";
     if(isimtk->type()==13||isimtk->type()==-13) {
     if (thetrackIdDone.find(isimtk->trackId()) != thetrackIdDone.end()){ continue;}
     }
     }
     //there wil be nMuon since in the event. could use this as the event-weight
     if (nMuon!=0){eventWeight=1./(double)nMuon;}
   */

   uint countMuons=0;
   thetrackIdDone.clear();
   for(std::vector<SimTrack>::const_iterator isimtk = SimTk->begin();isimtk!=SimTk->end();++isimtk)
     {
       LogDebug(theCategory)<<"I am here 1";
       if(isimtk->type()==13||isimtk->type()==-13)
	 {
	   //make sure that two muons with the same track id don't appear
	   if (thetrackIdDone.find(isimtk->trackId()) != thetrackIdDone.end()){edm::LogError(theCategory)<<"This muon has the same track Id as another muon"; continue;}
	   thetrackIdDone[isimtk->trackId()]=1;
	   countMuons++;

	   //calculate motherhood
	   MotherSearch mother(&*isimtk, SimTk, SimVtx, hepmc);
	
	   //use sim or gen information about mother
	   if (mother.SimIsValid()){
	     h.e("h_mu_sim_pt")->Fill(isimtk->momentum().perp(),eventWeight);
	     h.e("h_mu_sim_eta")->Fill(isimtk->momentum().eta(),eventWeight);
	     h.e("h_mu_sim_Aeta")->Fill(fabs(isimtk->momentum().eta()),eventWeight);
	     h.e("h_sim_vertex_position")->Fill(mother.Sim_vertex->position().v().mag(),eventWeight);

	     int motherBinNumber = wantMotherBin.GetBinNum(mother.Sim_mother->type());
	     h.e("h_parent_id")->Fill(motherBinNumber,eventWeight);
	     h.HBC["h_mu_sim_pt_assoc_ID"][motherBinNumber]->Fill(mother.simtrack->momentum().perp(),eventWeight);
	     h.HBC["h_mu_sim_eta_assoc_ID"][motherBinNumber]->Fill(mother.simtrack->momentum().eta(),eventWeight);
	     h.HBC["h_mu_sim_Aeta_assoc_ID"][motherBinNumber]->Fill(fabs(mother.simtrack->momentum().eta()),eventWeight);
	     h.HBC["h_sim_vertex_position_assoc_ID"][motherBinNumber]->Fill(mother.Sim_vertex->position().v().mag(),eventWeight);
	   }//use Sim
	   else if (mother.GenIsValid()){
	     h.e("h_mu_sim_pt")->Fill(isimtk->momentum().perp(),eventWeight);
	     h.e("h_mu_sim_eta")->Fill(isimtk->momentum().eta(),eventWeight);
	     h.e("h_mu_sim_Aeta")->Fill(fabs(isimtk->momentum().eta()),eventWeight);
	     h.e("h_sim_vertex_position")->Fill(mother.Gen_vertex->point3d().mag(),eventWeight);
	     
	     int motherBinNumber = wantMotherBin.GetBinNum(mother.Gen_mother->pdg_id());
	     h.e("h_parent_id")->Fill(motherBinNumber,eventWeight);
	     h.HBC["h_mu_sim_pt_assoc_ID"][motherBinNumber]->Fill(mother.gentrack->momentum().perp(),eventWeight);
	     h.HBC["h_mu_sim_eta_assoc_ID"][motherBinNumber]->Fill(mother.gentrack->momentum().eta(),eventWeight);
	     h.HBC["h_mu_sim_Aeta_assoc_ID"][motherBinNumber]->Fill(fabs(mother.gentrack->momentum().eta()),eventWeight);
	     h.HBC["h_sim_vertex_position_assoc_ID"][motherBinNumber]->Fill(mother.Gen_vertex->point3d().mag(),eventWeight);
	   }//both gen mother and vertex valid
	   else{
	     edm::LogWarning(theCategory)<<"tricky event. the muon has no parent whatsoever";
	     h.e("h_mu_sim_pt_tricky")->Fill(isimtk->momentum().perp(),eventWeight);
	     h.e("h_mu_sim_eta_tricky")->Fill(isimtk->momentum().eta(),eventWeight);
	   }
	 }//Simtrack is a muon
     }//loop over simtracks
   h.e("h_num_sim")->Fill(countMuons);
}


// ------------ method called once each job just before starting event loop  ------------
void 
MuonXRay::beginJob(const edm::EventSetup&)
{
  TH1* hist = MEA1<TH1>(h.e("h_parent_id"));
  for(int counting = 0;counting<wantMotherBin.size();counting++)
    {
      std::string binname = wantMotherBin.GetBinName(counting);
      hist->GetXaxis()->SetBinLabel(counting+1,binname.c_str());
    }
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonXRay::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonXRay);
