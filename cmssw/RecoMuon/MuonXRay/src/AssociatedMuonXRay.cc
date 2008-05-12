// -*- C++ -*-
//
// Package:    AssociatedMuonXRay
// Class:      AssociatedMuonXRay
// 
/**\class AssociatedMuonXRay AssociatedMuonXRay.cc Workspace/AssociatedMuonXRay/src/AssociatedMuonXRay.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Finn Rebassoo
//         Created:  Wed Oct  3 16:40:27 CDT 2007
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
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include <string>
#include <cstdlib>

#include <DataFormats/TrackReco/interface/Track.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <TrackingTools/PatternTools/interface/TSCPBuilderNoMaterial.h>


#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByChi2.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByPosition.h"

#include <MagneticField/Records/interface/IdealMagneticFieldRecord.h>
#include <MagneticField/Engine/interface/MagneticField.h>

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "HepMC/GenEvent.h"
#include "TLorentzVector.h"

#include <DataFormats/TrackReco/interface/Track.h>

#include "DataFormats/MuonReco/interface/MuonTrackLinks.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "HepPDT/ParticleDataTable.hh"

#include "RecoMuon/MuonXRay/interface/IDconverttoBinNum.h"
#include "RecoMuon/MuonXRay/interface/MotherSearch.h"  
#include "RecoMuon/MuonXRay/interface/DQMHelper.h"  

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

//#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/HLTFilterObject.h"

#include <DataFormats/MuonReco/interface/MuonTrackLinks.h>

//
// class declaration
//

class AssociatedMuonXRay : public edm::EDAnalyzer {
public:
  explicit AssociatedMuonXRay(const edm::ParameterSet&);
  ~AssociatedMuonXRay();
  

private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  // ----------member data ---------------------------
  std::string theCategory;
  edm::InputTag thetrackLabel;
  edm::InputTag trackingParticleLabel;

  std::string theAssocLabel;
  edm::ESHandle<TrackAssociatorBase> theAssociator;

  bool theConfirm;
  edm::InputTag theConfirmLabel;

  bool theNeedLink;
  edm::InputTag theLinkLabel;

  edm::ESHandle<MagneticField> field;
  edm::ESHandle<HepPDT::ParticleDataTable> fPDGTable;

  IDconverttoBinNum wantMotherBin;
  DQMHelper h;
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

double pt90(const reco::TrackRef & tk, const edm::Event & ev){
  std::string prov=ev.getProvenance(tk.id()).moduleLabel();
  std::string instance = ev.getProvenance(tk.id()).productInstanceName();

  double nsigma_Pt=0;
  //  std::cout<<" looking at a track from: "<<prov<<":"<<instance<<std::endl;

  if (prov.find("L2")!= std::string::npos){ nsigma_Pt=3.9;}
  else if (prov.find("L3")!= std::string::npos){
    if (instance=="") nsigma_Pt=2.2;
    else nsigma_Pt=1.64; //90% for sure
  }
  else { edm::LogError("p90")<<"provenance: "<<prov<<" instance not recognized.";}

  double pt = tk->pt();
  double err0 = tk->error(0);
  double abspar0 = fabs(tk->parameter(0));
  double ptLx = pt;
  if (abspar0>0) ptLx += nsigma_Pt*err0/abspar0*pt;

  //  std::cout<<"initial pt: "<<pt<<", corrected pt: "<<ptLx<<std::endl;

  return ptLx;
}




AssociatedMuonXRay::AssociatedMuonXRay(const edm::ParameterSet& iConfig):
  theCategory("AssociatedMuonXRay"),
  wantMotherBin(iConfig.getParameter<edm::ParameterSet>("IDconverttoBinNum"))
{

  thetrackLabel = iConfig.getParameter<edm::InputTag>("trackLabel");
  trackingParticleLabel = iConfig.getParameter<edm::InputTag>("trackingParticleLabel");

  theAssocLabel = iConfig.getParameter<std::string>("associatorName");

  theConfirm = iConfig.exists("confirmLabel");
  if (theConfirm){
    theConfirmLabel = iConfig.getParameter<edm::InputTag>("confirmLabel");
    theNeedLink = iConfig.exists("linkLabel");
    if (theNeedLink)
      theLinkLabel = iConfig.getParameter<edm::InputTag>("linkLabel");
  }


  h.dbe()->setCurrentFolder("AssociatedMuonXRay/"+iConfig.getParameter<std::string>("@module_label"));
  //now that the directory is set
  h = DQMHelper(iConfig.getParameter<edm::ParameterSet>("DQMHelper"));

  wantMotherBin.splitH1ByCategory("h_mu_pt",h,"%s coming from %s");
  wantMotherBin.splitH1ByCategory("h_mu_90pt",h,"%s coming from %s");
  wantMotherBin.splitH1ByCategory("h_mu_leading_pt",h,"%s coming from %s");
  wantMotherBin.splitH1ByCategory("h_mu_leading_90pt",h,"%s coming from %s");
  wantMotherBin.splitH1ByCategory("h_mu_hits",h,"%s coming from %s");
  wantMotherBin.splitH1ByCategory("h_mu_leading_Aeta",h,"%s coming from %s");
  wantMotherBin.splitH1ByCategory("h_mu_eta",h,"%s coming from %s");
  wantMotherBin.splitH1ByCategory("h_mu_Aeta",h,"%s coming from %s");
  wantMotherBin.splitH1ByCategory("h_mu_d0",h,"%s coming from %s");

  wantMotherBin.splitH1ByCategory("h_mu_sim_pt",h,"%s coming from %s");
  wantMotherBin.splitH1ByCategory("h_mu_leading_sim_pt",h,"%s coming from %s");
  wantMotherBin.splitH1ByCategory("h_mu_sim_eta",h,"%s coming from %s");
  wantMotherBin.splitH1ByCategory("h_mu_sim_Aeta",h,"%s coming from %s");

  h.HBC["h_part_sim_pt_assoc_ID"].resize(wantMotherBin.size());
  h.HBC["h_part_pt_assoc_ID"].resize(wantMotherBin.size());
  h.HBC["quality_assoc_ID"].resize(wantMotherBin.size());

    //by categories of mother
  std::string htitle,hname,binname;
  for(int iut = 0;iut<wantMotherBin.size();iut++)
    {
      binname = wantMotherBin.GetBinName(iut);
      //whatever is associated
      htitle = Form("p_{T}^{sim} of associated %s",binname.c_str());
      hname = Form("h_part_sim_pt_assoc_ID_%i",iut);
      h.book1D(iut,"h_part_sim_pt_assoc_ID",hname,htitle,400,0,200,"particle p_{T}^{sim} [GeV]");

      htitle = Form("muon p_{T}^{reco} when associated to %s",binname.c_str());
      hname = Form("h_part_pt_assoc_ID_%i",iut);
      h.book1D(iut,"h_part_pt_assoc_ID",hname,htitle,400,0,200,"muon p_{T}^{reco} [GeV]");

      htitle = Form("association quality of reco::Track to %s",binname.c_str());
      hname = Form("quality_assoc_ID_%i",iut);
      h.book1D(iut,"quality_assoc_ID",hname,htitle,400,-2,2,"quality");
    }

}


AssociatedMuonXRay::~AssociatedMuonXRay()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

	    
//
// member functions
//
	    
	    
	    
// ------------ method called to for each event  ------------
	    
void
AssociatedMuonXRay::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
   
   //get the mag field
   iSetup.get<IdealMagneticFieldRecord>().get(field);
   
   //open a collection of track
   edm::Handle<reco::TrackCollection> tracks;
   iEvent.getByLabel(thetrackLabel,tracks);

   //open a collection of Tracking Particles
   edm::Handle<TrackingParticleCollection> TPtracks;
   iEvent.getByLabel(trackingParticleLabel,TPtracks);

   //get a hold off simulated information
   Handle<SimTrackContainer> SimTk;
   Handle<SimVertexContainer> SimVtx;
   //get a hold on generated information
   Handle<HepMCProduct> hepmc;

   iEvent.getByType(hepmc);
   iEvent.getByLabel("g4SimHits",SimVtx);
   iEvent.getByLabel("g4SimHits",SimTk);
	    
   //get the associators
   iSetup.get<TrackAssociatorRecord>().get(theAssocLabel,theAssociator);
   
   //associate
   reco::RecoToSimCollection recSimColl = theAssociator->associateRecoToSim(tracks,TPtracks, &iEvent);

   LogDebug(theCategory)<<"I have found: "<<recSimColl.size()<<" associations in total.";   
   //loop over tracks
   uint iT =0;
   uint iTmax = tracks->size();
   h.e("num_reco")->Fill(iTmax);

   //trigger collections
   edm::Handle<reco::HLTFilterObjectWithRefs> triggered;
   edm::Handle<reco::MuonTrackLinksCollection> links;
   if (theConfirm){
     iEvent.getByLabel(theConfirmLabel, triggered);
     if (theNeedLink) iEvent.getByLabel(theLinkLabel, links);
   }

   double highest_pt_reco = -1;
   reco::TrackRef leading_track;
   const SimTrack * leading_simtrack=0, *leading_simMother=0;;
   const HepMC::GenParticle * leading_gentrack=0, *leading_genMother=0;

   uint countConfirm=0;
   for(;iT!=iTmax;++iT){
     reco::TrackRef ref(tracks, iT);

     //confirm that it has triggered.
     bool isConfirmed=false;
     if (theConfirm){
       /*
               //cannot do any of that easilty in 167...
	       std::vector<reco::RecoChargedCandidateRef> cands;
	       triggered->getObjects(trigger::TriggerMuon,cands);
	       std::vector<reco::RecoChargedCandidateRef>::iterator it= cands.begin();
	       std::vector<reco::RecoChargedCandidateRef>::iterator end= cands.end();
	       for (;ti!=end;++ti){
	       if (ti->get<reco::TrackRef>() == ref){
	       isConfirmed=true; break;
	       }
	       }
       */
       reco::HLTFilterObjectWithRefs::const_iterator ti = triggered->begin();
       reco::HLTFilterObjectWithRefs::const_iterator end = triggered->end();
       reco::TrackRef refToCompareTo=ref;
       if (theNeedLink){
	 //find link back to tracker track
	 reco::MuonTrackLinksCollection::const_iterator il= links->begin();
	 reco::MuonTrackLinksCollection::const_iterator end= links->end();
	 for (;il!=end;++il){
	   if (il->trackerTrack() == refToCompareTo){
	     refToCompareTo=il->globalTrack(); break;}
	 }
       }
       
       for (;ti!=end;++ti){
	 const reco::RecoChargedCandidate * cand = dynamic_cast<const reco::RecoChargedCandidate*>(&(*ti));
	 if (!cand) { edm::LogError(theCategory)<<" casting is failing."; continue;}
	 //	 edm::LogInfo(theCategory)<<"comparing references from: "<<cand->track().id()<<" and: "<<refToCompareTo.id();
	 if (cand->track() == refToCompareTo){
	   isConfirmed=true; break; }
       }

     }
     if (theConfirm && !isConfirmed) continue;
     countConfirm++;

     bool isLeading=false;
     if(ref->pt()>highest_pt_reco) {
       highest_pt_reco= ref->pt();
       leading_track = ref;
       isLeading=true;}

     h.e("h_mu_pt")->Fill(ref->pt());
     h.e("h_mu_90pt")->Fill(pt90(ref,iEvent));
     h.e("h_mu_hits")->Fill(ref->recHitsSize());
     h.e("h_mu_eta")->Fill(ref->eta());
     h.e("h_mu_Aeta")->Fill(fabs(ref->eta()));
     h.e("h_mu_d0")->Fill(ref->d0());
     h.e("h_mu_d0_pT")->Fill(ref->pt(),ref->d0());
     
     reco::RecoToSimCollection::const_iterator findRefIt = recSimColl.find(ref);
     if (findRefIt != recSimColl.end()){
       const std::vector<std::pair<TrackingParticleRef,double> > & tp = findRefIt->val;
       const TrackingParticleRef & trp = tp.begin()->first;

       int particle_ID = trp->pdgId();
       int myBin = wantMotherBin.GetBinNum(particle_ID);
       
       h.HBC["quality_assoc_ID"][myBin]->Fill(tp.begin()->second);       
       h.HBC["h_part_sim_pt_assoc_ID"][myBin]->Fill(trp->pt());
       h.HBC["h_part_pt_assoc_ID"][myBin]->Fill(ref->pt());
       
       if(abs(particle_ID) == 13){
	 h.e("h_mu_sim_pt")->Fill(trp->pt());
	 h.e("h_mu_sim_Aeta")->Fill(fabs(trp->eta()));
	 h.e("h_mu_sim_eta")->Fill(trp->eta());
	 
	 //---------------------- MOTHERHOOD --------------------------------
	 //find the parent of tracking particle
	 for(TrackingParticle::g4t_iterator isimtk = trp->g4Track_begin();isimtk!=trp->g4Track_end();isimtk++)
	   {
	     LogDebug(theCategory)<<"I am here 1";      
	     if(isimtk->type()==13||isimtk->type()==-13)
	       {
		 //calculate mother hood
		 MotherSearch mother(&*isimtk, SimTk, SimVtx, hepmc);

		 //FIXME, use reco::Particle mother.mother();
		 double pt,eta;
		 int parentID;
		 int motherBinNumber;
		 if (mother.IsValid()){
		   if (mother.SimIsValid()){
		     parentID= mother.Sim_mother->type();
		     motherBinNumber = wantMotherBin.GetBinNum(mother.Sim_mother->type());
		     pt=mother.simtrack->momentum().perp();
		     eta=mother.simtrack->momentum().eta();
		     if (isLeading){
		       leading_simtrack = mother.simtrack;
		       leading_gentrack = 0;
		       leading_genMother=0;
		       leading_simMother=mother.Sim_mother;}
		   }
		   else {
		     parentID= mother.Gen_mother->pdg_id();
		     motherBinNumber = wantMotherBin.GetBinNum(mother.Gen_mother->pdg_id());
		     pt=mother.gentrack->momentum().perp();
		     eta=mother.gentrack->momentum().eta();
		     if (isLeading){
		       leading_simtrack = 0;
		       leading_gentrack = mother.gentrack;
		       leading_genMother= mother.Gen_mother;
		       leading_simMother=0;}
		   }

		   //fill histograms
		   h.HBC["h_mu_pt_assoc_ID"][motherBinNumber]->Fill(ref->pt());
		   h.HBC["h_mu_90pt_assoc_ID"][motherBinNumber]->Fill(pt90(ref,iEvent));
		   h.HBC["h_mu_hits_assoc_ID"][motherBinNumber]->Fill(ref->recHitsSize());
		   h.HBC["h_mu_eta_assoc_ID"][motherBinNumber]->Fill(ref->eta());
		   h.HBC["h_mu_Aeta_assoc_ID"][motherBinNumber]->Fill(fabs(ref->eta()));
		   h.HBC["h_mu_d0_assoc_ID"][motherBinNumber]->Fill(ref->d0());

		   h.HBC["h_mu_sim_pt_assoc_ID"][motherBinNumber]->Fill(pt);
		   h.HBC["h_mu_sim_eta_assoc_ID"][motherBinNumber]->Fill(eta);
		   h.HBC["h_mu_sim_Aeta_assoc_ID"][motherBinNumber]->Fill(fabs(eta));

		   //do it once per tracking particle once it succeed
		   break;
		 }
		 else{edm::LogError(theCategory)<<"tricky muon from TrackingParticle.";}
		
			
	       }//sim track is a muon
	     else{  edm::LogError(theCategory)<<"the sim track attached to the tracking particle is not a muon.";}
	   }//loop over SimTrack of tracking particle
       }//muon associated
       else{
	 //a reco muon is associated to something else than a muon
	 edm::LogError(theCategory)<<"a reconstructed muon is associated to: "<<particle_ID;
       }
     }//track has an association
     else{
       //this track was not associated.
       edm::LogError(theCategory)<<"a reconstructed muon is not associated.";
     }
     
   }//loop over tracks


   h.e("num_reco_confirm")->Fill(countConfirm);

   //make a plot for the leading track
   if (!leading_track.isNull()){
     h.e("h_mu_leading_pt")->Fill(leading_track->pt());   
     h.e("h_mu_leading_90pt")->Fill(pt90(leading_track,iEvent));   
     h.e("h_mu_leading_Aeta")->Fill(fabs(leading_track->eta()));
     if (leading_simtrack){
       h.e("h_mu_leading_sim_pt")->Fill(leading_simtrack->momentum().perp());
       int motherBinNumber = wantMotherBin.GetBinNum(leading_simMother->type());
       h.HBC["h_mu_leading_sim_pt_assoc_ID"][motherBinNumber]->Fill(leading_simtrack->momentum().perp());
       h.HBC["h_mu_leading_pt_assoc_ID"][motherBinNumber]->Fill(leading_track->pt());
       h.HBC["h_mu_leading_90pt_assoc_ID"][motherBinNumber]->Fill(pt90(leading_track,iEvent));
     }
     else if(leading_gentrack){
       h.e("h_mu_leading_sim_pt")->Fill(leading_gentrack->momentum().perp());
       int motherBinNumber = wantMotherBin.GetBinNum(leading_genMother->pdg_id());
       h.HBC["h_mu_leading_sim_pt_assoc_ID"][motherBinNumber]->Fill(leading_gentrack->momentum().perp());
       h.HBC["h_mu_leading_pt_assoc_ID"][motherBinNumber]->Fill(leading_track->pt());
       h.HBC["h_mu_leading_90pt_assoc_ID"][motherBinNumber]->Fill(pt90(leading_track,iEvent));
     }
   }
      
   /*
   //associate sim to reco
   reco::SimToRecoCollection simRecoColl;
   simRecoColl = theAssociator->associateSimToReco(tracks,TPtracks, &iEvent);
   //loop over map
   for(TrackingParticleCollection::size_type i=0;i<TPtracks->size();++i)
     {
       TrackingParticleRef tp(TPtracks,i);
       std::vector<std::pair<reco::TrackRef,double> > recoMuVec;
       reco::TrackRef recoMu;
       if(simRecoColl.find(tp) !=simRecoColl.end())
	 {recoMuVec = simRecoColl[tp];
	   if(recoMuVec.size()!=0){
	     recoMu = recoMuVec.begin()->first;
	     //	     double simptMu = tp->pt();
	     //	       double simptEta = tp->eta();
	     //	       double recoptMu = recoMu->pt();
	     //	       double recoetaMu = recoMu->eta();
	     //--	     sim_pt_sim_to_reco->Fill(tp->pt());
	     //--	     sim_eta_sim_to_reco->Fill(tp->eta());
	     //--	     reco_pt_sim_to_reco->Fill(recoMu->pt());
	     //--	     reco_eta_sim_to_reco->Fill(recoMu->eta());
	     if(tp->pdgId()==13||tp->pdgId()==-13)
	       {
		 //--		 sim_pt_sim_to_reco_mu->Fill(tp->pt());
		 //--		 sim_eta_sim_to_reco_mu->Fill(tp->eta());
	       }
	     else{
	       //--	       sim_pt_sim_to_reco_not_mu->Fill(tp->pt());
	       //--	       sim_eta_sim_to_reco_not_mu->Fill(tp->eta());
	     }
	     
	   }
	 }
     }
   */
}


// ------------ method called once each job just before starting event loop  ------------
void 
AssociatedMuonXRay::beginJob(const edm::EventSetup&)
{

}

// ------------ method called once each job just after ending the event loop  ------------
void 
AssociatedMuonXRay::endJob() {
  /*

  int NXbins = h_muon_pt_reco_assoc->GetNbinsX();
  for(int count = 0;count<NXbins;count++)
    { double total_integral_value = h_muon_pt_reco_assoc->Integral(0,NXbins);
      double integral_value = h_muon_pt_reco_assoc->Integral(0,count+1);
      h_muon_pt_reco_assoc_cumul->SetBinContent(count+1,integral_value/total_integral_value);
    }

  int NXbins_non_assoc = h_muon_pt_reco->GetNbinsX();
  for(int count = 0;count<NXbins_non_assoc;count++)
    {double total_integral_value = h_muon_pt_reco_assoc->Integral(0,NXbins_non_assoc);
      double integral_value = h_muon_pt_reco_assoc->Integral(0,count+1);
      h_muon_pt_reco_cumul->SetBinContent(count+1,integral_value/total_integral_value);
    }
  */



}

//define this as a plug-in
DEFINE_FWK_MODULE(AssociatedMuonXRay);