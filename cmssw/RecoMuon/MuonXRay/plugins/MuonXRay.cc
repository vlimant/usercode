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
// $Id: MuonXRay.cc,v 1.3 2008/08/27 06:56:36 vlimant Exp $
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
  wantMotherBin.splitH1ByCategory("h_mu_sim_phi",h,"%s coming from %s");
  wantMotherBin.splitH1ByCategory("h_mu_sim_Aeta",h,"%s coming from %s");
  wantMotherBin.splitH1ByCategory("h_sim_vertex_position",h,"%s coming from %s");           

  wantMotherBin.splitH1ByCategory("h_leading_mu_sim_pt",h,"%s coming from %s");
  wantMotherBin.splitH1ByCategory("h_leading_mu_sim_eta",h,"%s coming from %s");
  wantMotherBin.splitH1ByCategory("h_leading_mu_sim_phi",h,"%s coming from %s");
  wantMotherBin.splitH1ByCategory("h_leading_mu_sim_Aeta",h,"%s coming from %s");
  wantMotherBin.splitH1ByCategory("h_leading_mu_sim_vertex_position",h,"%s coming from %s");           


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

   double leadingPt=0;
   std::vector<SimTrack>::const_iterator leadingSimTrack = SimTk->end();
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

	   if (isimtk->momentum().pt() > leadingPt){
	     leadingPt = isimtk->momentum().pt();
	     leadingSimTrack = isimtk;
	   }

	   //calculate motherhood
	   MotherSearch mother(&*isimtk, SimTk, SimVtx, hepmc);
	
	   //use sim or gen information about mother
	   if (mother.SimIsValid()){
	     h.e("h_mu_sim_pt")->Fill(isimtk->momentum().pt(),eventWeight);
	     h.e("h_mu_sim_eta")->Fill(isimtk->momentum().eta(),eventWeight);
	     h.e("h_mu_sim_phi")->Fill(isimtk->momentum().phi(),eventWeight);
	     h.e("h_mu_sim_Aeta")->Fill(fabs(isimtk->momentum().eta()),eventWeight);
	     h.e("h_sim_vertex_position")->Fill(mother.Sim_vertex->position().mag(),eventWeight);

	     int motherBinNumber = wantMotherBin.GetBinNum(mother.Sim_mother->type());
	     h.e("h_parent_id")->Fill(motherBinNumber,eventWeight);
	     h.HBC["h_mu_sim_pt_assoc_ID"][motherBinNumber]->Fill(mother.simtrack->momentum().pt(),eventWeight);
	     h.HBC["h_mu_sim_eta_assoc_ID"][motherBinNumber]->Fill(mother.simtrack->momentum().eta(),eventWeight);
	     h.HBC["h_mu_sim_phi_assoc_ID"][motherBinNumber]->Fill(mother.simtrack->momentum().phi(),eventWeight);
	     h.HBC["h_mu_sim_Aeta_assoc_ID"][motherBinNumber]->Fill(fabs(mother.simtrack->momentum().eta()),eventWeight);
	     h.HBC["h_sim_vertex_position_assoc_ID"][motherBinNumber]->Fill(mother.Sim_vertex->position().mag(),eventWeight);
	   }//use Sim
	   else if (mother.GenIsValid()){
	     h.e("h_mu_sim_pt")->Fill(isimtk->momentum().pt(),eventWeight);
	     h.e("h_mu_sim_eta")->Fill(isimtk->momentum().eta(),eventWeight);
	     h.e("h_mu_sim_phi")->Fill(isimtk->momentum().phi(),eventWeight);
	     h.e("h_mu_sim_Aeta")->Fill(fabs(isimtk->momentum().eta()),eventWeight);
	     h.e("h_sim_vertex_position")->Fill(mother.Gen_vertex->point3d().mag(),eventWeight);
	     
	     int motherBinNumber = wantMotherBin.GetBinNum(mother.Gen_mother->pdg_id());
	     h.e("h_parent_id")->Fill(motherBinNumber,eventWeight);
	     h.HBC["h_mu_sim_pt_assoc_ID"][motherBinNumber]->Fill(mother.gentrack->momentum().perp(),eventWeight);
	     h.HBC["h_mu_sim_eta_assoc_ID"][motherBinNumber]->Fill(mother.gentrack->momentum().eta(),eventWeight);
	     h.HBC["h_mu_sim_phi_assoc_ID"][motherBinNumber]->Fill(mother.gentrack->momentum().phi(),eventWeight);
	     h.HBC["h_mu_sim_Aeta_assoc_ID"][motherBinNumber]->Fill(fabs(mother.gentrack->momentum().eta()),eventWeight);
	     h.HBC["h_sim_vertex_position_assoc_ID"][motherBinNumber]->Fill(mother.Gen_vertex->point3d().mag(),eventWeight);
	   }//both gen mother and vertex valid
	   else{
	     edm::LogWarning(theCategory)<<"tricky event. the muon has no parent whatsoever";
	     h.e("h_mu_sim_pt_tricky")->Fill(isimtk->momentum().pt(),eventWeight);
	     h.e("h_mu_sim_eta_tricky")->Fill(isimtk->momentum().eta(),eventWeight);
	   }
	 }//Simtrack is a muon
     }//loop over simtracks
   
   //fill things for the leading track now.
   if (leadingSimTrack!=SimTk->end()){
     
	   //calculate motherhood
	   MotherSearch mother(&*leadingSimTrack, SimTk, SimVtx, hepmc);
	
	   //use sim or gen information about mother
	   if (mother.SimIsValid()){
	     h.e("h_leading_mu_sim_pt")->Fill(leadingSimTrack->momentum().pt(),eventWeight);
	     h.e("h_leading_mu_sim_eta")->Fill(leadingSimTrack->momentum().eta(),eventWeight);
	     h.e("h_leading_mu_sim_phi")->Fill(leadingSimTrack->momentum().phi(),eventWeight);
	     h.e("h_leading_mu_sim_Aeta")->Fill(fabs(leadingSimTrack->momentum().eta()),eventWeight);
	     h.e("h_leading_mu_sim_vertex_position")->Fill(mother.Sim_vertex->position().mag(),eventWeight);

	     int motherBinNumber = wantMotherBin.GetBinNum(mother.Sim_mother->type());
	     h.e("h_parent_id")->Fill(motherBinNumber,eventWeight);
	     h.HBC["h_leading_mu_sim_pt_assoc_ID"][motherBinNumber]->Fill(mother.simtrack->momentum().pt(),eventWeight);
	     h.HBC["h_leading_mu_sim_eta_assoc_ID"][motherBinNumber]->Fill(mother.simtrack->momentum().eta(),eventWeight);
	     h.HBC["h_leading_mu_sim_phi_assoc_ID"][motherBinNumber]->Fill(mother.simtrack->momentum().phi(),eventWeight);
	     h.HBC["h_leading_mu_sim_Aeta_assoc_ID"][motherBinNumber]->Fill(fabs(mother.simtrack->momentum().eta()),eventWeight);
	     h.HBC["h_leading_mu_sim_vertex_position_assoc_ID"][motherBinNumber]->Fill(mother.Sim_vertex->position().mag(),eventWeight);
	   }//use Sim
	   else if (mother.GenIsValid()){
	     h.e("h_leading_mu_sim_pt")->Fill(leadingSimTrack->momentum().pt(),eventWeight);
	     h.e("h_leading_mu_sim_eta")->Fill(leadingSimTrack->momentum().eta(),eventWeight);
	     h.e("h_leading_mu_sim_phi")->Fill(leadingSimTrack->momentum().phi(),eventWeight);
	     h.e("h_leading_mu_sim_Aeta")->Fill(fabs(leadingSimTrack->momentum().eta()),eventWeight);
	     h.e("h_leading_mu_sim_vertex_position")->Fill(mother.Gen_vertex->point3d().mag(),eventWeight);
	     
	     int motherBinNumber = wantMotherBin.GetBinNum(mother.Gen_mother->pdg_id());
	     h.e("h_parent_id")->Fill(motherBinNumber,eventWeight);
	     h.HBC["h_leading_mu_sim_pt_assoc_ID"][motherBinNumber]->Fill(mother.gentrack->momentum().perp(),eventWeight);
	     h.HBC["h_leading_mu_sim_eta_assoc_ID"][motherBinNumber]->Fill(mother.gentrack->momentum().eta(),eventWeight);
	     h.HBC["h_leading_mu_sim_phi_assoc_ID"][motherBinNumber]->Fill(mother.gentrack->momentum().phi(),eventWeight);
	     h.HBC["h_leading_mu_sim_Aeta_assoc_ID"][motherBinNumber]->Fill(fabs(mother.gentrack->momentum().eta()),eventWeight);
	     h.HBC["h_leading_mu_sim_vertex_position_assoc_ID"][motherBinNumber]->Fill(mother.Gen_vertex->point3d().mag(),eventWeight);
	   }//both gen mother and vertex valid
	   else{
	     edm::LogWarning(theCategory)<<"tricky event. the muon has no parent whatsoever";
	     h.e("h_leading_mu_sim_pt_tricky")->Fill(leadingSimTrack->momentum().pt(),eventWeight);
	     h.e("h_leading_mu_sim_eta_tricky")->Fill(leadingSimTrack->momentum().eta(),eventWeight);
	   }


   }

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
