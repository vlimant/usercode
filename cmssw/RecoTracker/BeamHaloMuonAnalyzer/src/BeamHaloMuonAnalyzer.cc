// -*- C++ -*-
//
// Package:    BeamHaloMuonAnalyzer
// Class:      BeamHaloMuonAnalyzer
// 
/**\class BeamHaloMuonAnalyzer BeamHaloMuonAnalyzer.cc RecoTracker/BeamHaloMuonAnalyzer/src/BeamHaloMuonAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jean-Roch Vlimant
//         Created:  Tue Apr  3 19:51:45 CDT 2007
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

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <RecoMuon/L3MuonAnalyzer/interface/Analyzer.h>
#include <SimDataFormats/Track/interface/SimTrackContainer.h>
#include <SimDataFormats/Vertex/interface/SimVertexContainer.h>
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include <TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h>
#include <DataFormats/GeometrySurface/interface/Plane.h>
#include <DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DCollection.h>

#include "TFile.h"
#include "TTree.h"
#include "RecoMuon/L3MuonAnalyzer/interface/DumpClass.h"

//
// class decleration
//

class BeamHaloMuonAnalyzer : public edm::EDAnalyzer {
   public:
      explicit BeamHaloMuonAnalyzer(const edm::ParameterSet&);
      ~BeamHaloMuonAnalyzer();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
  edm::InputTag seedGeneratorLabel;
  edm::InputTag trackProducerLabel;
  std::string _category;

  IntrusiveAnalyzer  _analyzer;

  std::string _fileName;
  TFile * _file;
  TTree * _tree;
  EventDump * _dump;
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
BeamHaloMuonAnalyzer::BeamHaloMuonAnalyzer(const edm::ParameterSet& iConfig)
{
  _category = "BeamHaloMuonAnalyzer";
  seedGeneratorLabel = iConfig.getParameter<edm::InputTag>("seedGeneratorLabel");
  trackProducerLabel = iConfig.getParameter<edm::InputTag>("trackProducerLabel");

  //create a roottree
  _fileName=iConfig.getParameter<std::string>("@module_label")+".root";

}


BeamHaloMuonAnalyzer::~BeamHaloMuonAnalyzer()
{
 
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
BeamHaloMuonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   //setup the analyzer
   _analyzer.setUp(iEvent,iSetup);

   //initialize the tree
   _dump->clean();
   bool fillMe=false;
   _dump->run = (int)iEvent.id().run();
   _dump->evt = (int)iEvent.id().event();
   
   //retreive sim tracks
   Handle<SimTrackContainer> simTracks;
   iEvent.getByType<SimTrackContainer>(simTracks);
   Handle<SimVertexContainer> simVertex;
   iEvent.getByType<SimVertexContainer>(simVertex);
   
   //retreive seeds from reconstruction
   Handle<TrajectorySeedCollection> seeds;
   iEvent.getByLabel(seedGeneratorLabel, seeds);
   
   //retreive tracks from reconstruction
   Handle<reco::TrackCollection> tracks;
   iEvent.getByLabel(trackProducerLabel, tracks);

   //build a plane right in the middle of the detector
   const GlobalPoint center(0,0,0);
   //axis are the same as the global frame
   Surface::RotationType rotation;
   Plane::PlanePointer medianPlane = Plane::build(center,rotation);

   //show what is in the simtracks
   edm::LogInfo(_category)<<"simulated track content ("<< simTracks->size()<<") SimTrack";
   //   std::vector<SimTrackRef> mc_muons=_analyzer.particles(13,0);
   for (SimTrackContainer::const_iterator simIt = simTracks->begin(); simIt!=simTracks->end();++simIt){
     if (abs(simIt->type())!=13) continue; //muon only

     //create a candidate
     RoadCandidate * candidate = _dump->add();
     fillMe=true;

     const HepLorentzVector& mom= simIt->momentum();
     GlobalVector momentum(mom.x(),mom.y(),mom.z());
     GlobalPoint position;
     if (simIt->vertIndex()!=-1){
       const HepLorentzVector& pos = (*simVertex)[simIt->vertIndex()].position();
       position=  GlobalPoint(pos.x(),pos.y(),pos.z());}

     //check if you can get a muon trajectorystate
     TrajectoryStateOnSurface muonState;
     TrajectoryStateOnSurface medianState,Tk_medianState;

     //get the list of tracker simthits.
     std::map< double ,TrackPSimHitRef > mc_hit_Htof=_analyzer.MapRefCount(simIt->trackId(),true,false,true,"HighTof");
     std::map< double ,TrackPSimHitRef > mc_hit_Ltof=_analyzer.MapRefCount(simIt->trackId(),true,false,true,"LowTof");
     
	     
     if (mc_hit_Htof.size()!=0) {
       //convert Htof PSimHit
       muonState = _analyzer.convert(*mc_hit_Htof.begin()->second);
       if (muonState.isValid()){edm::LogVerbatim(_category)<<"this track has a state on the tracker system (Htof):\n"<<muonState; }
     }
     else if (mc_hit_Ltof.size()!=0) {
       //convert Ltof PSimHit
       muonState = _analyzer.convert(*mc_hit_Ltof.begin()->second);
       if (muonState.isValid()){edm::LogVerbatim(_category)<<"this track has a state on the tracker system (Ltof):\n"<<muonState; }
     }
     else {
       //get a PSimHit in the muon system
       muonState = _analyzer.convertFromMuonSystem(simIt->trackId());     
       if (muonState.isValid()){edm::LogVerbatim(_category)<<"this track has a state on the muon system:\n"<<muonState; }
     }

     if (muonState.isValid()) {
       //you have a muon state
       //get it at the median plane
       medianState = _analyzer.prop(anyDirection)->propagate(*muonState.freeState(),*medianPlane);
	 
       if (medianState.isValid()) {
	 //it reached the median plane. tell me where it is
	 double cosine_angle=medianState.globalMomentum().unit().z(); 
	 const double PI=(3.1415926);
	 double angle = acos(cosine_angle)*180/PI;
	 double r=medianState.globalPosition().perp();
	   
	 edm::LogVerbatim(_category)<<"this muon can reach the median plane: \n"
				    <<"cosine(angle): "<<cosine_angle<<"\n"
				    <<"angle: "<<angle<<"\n"
				    <<"180-angle: "<<180-angle<<"\n"
				    <<"pz: "<<medianState.globalMomentum().z()<<"\n"
				    <<"pT: "<<medianState.globalMomentum().perp()<<"\n"
				    <<"radius: "<<r<<"\n"
				    <<"state:\n"<<medianState;
	 //set the tree content
	 candidate->hasSim=true;
	 candidate->sim<<medianState;
	 
	 //check whether there is a rough track match
	 const double dMag_default = 10000 ; //cm
	 double dMag=dMag_default;
	 const reco::Track * best=0;
	 for (reco::TrackCollection::const_iterator Tkit = tracks->begin(); Tkit!= tracks->end();++Tkit){
	   const reco::Track & muon = *Tkit;
	   
	   //build a transient track from this guy
	   reco::TransientTrack Tmuon = _analyzer.TTbuilder().build(&muon);
	   
	   TrajectoryStateOnSurface tmp_Tk_medianState = _analyzer.prop(anyDirection)->propagate(*Tmuon.outermostMeasurementState().freeState(),
											     *medianPlane);
	   if (!tmp_Tk_medianState.isValid()) continue;
	   if ((tmp_Tk_medianState.globalPosition() - medianState.globalPosition()).mag()< dMag){
	     Tk_medianState = tmp_Tk_medianState;
	     dMag = (tmp_Tk_medianState.globalPosition() - medianState.globalPosition()).mag();
	     best= &(*Tkit);}
	 }//loop over the reco tracks

	 if (Tk_medianState.isValid()){
	   double Tk_cosine_angle=Tk_medianState.globalMomentum().unit().z();
	   double Tk_angle = acos(Tk_cosine_angle)*180/PI;
	   double Tk_r = Tk_medianState.globalPosition().perp();
	  
	   //set the tree content
	   candidate->hasBest=true;
	   candidate->best.updated<<Tk_medianState;
	   trackingRecHit_iterator rh = best->recHitsBegin();
	   for (; rh !=best->recHitsEnd();++rh){
	     if (!(*rh)->isValid()){continue;}
	     std::pair<IntrusiveAnalyzer::RecHit_match_PSimHit, TrajectoryStateOnSurface> matchHit;
	     //match in Htof (last argument to false)
	     matchHit = _analyzer.match(**rh,false,simIt->trackId(),false);
	     //try Ltof if fails in Htof
	     if (!matchHit.second.isValid()) matchHit = _analyzer.match(**rh,false,simIt->trackId(),true);
	     DressedHit * hit = candidate->best.add();
	     hit->matchcode = matchHit.first;
	     TransientTrackingRecHit::RecHitPointer Trh = _analyzer.TTRHbuilder().build(&(**rh));
	     if (!matchHit.second.isValid()){
		 if( !Trh->isValid()){continue;}
	       hit->simhit<<(*Trh);}
	     else{hit->simhit<<matchHit.second;}
	     if (Trh->isValid()){hit->rechit<<(*Trh);}
	   }//loop over own rechits

	   edm::LogVerbatim(_category)<<"probable match to a reconstructed track:\n"
				      <<"distance in the median plan: "<<dMag<<"\n"
				      <<"cosine(angle): "<<Tk_cosine_angle<<"\n"
				      <<"angle: "<<Tk_angle<<"\n"
				      <<"180-angle: "<<180-Tk_angle<<"\n"
				      <<"pz: "<<Tk_medianState.globalMomentum().z()<<"\n"
				      <<"pT: "<<Tk_medianState.globalMomentum().perp()<<"\n"
				      <<"radius: "<<Tk_r<<"\n"
				      <<"with ("<< best->recHitsSize()<<") rechits \n"
				      <<"state:\n"<<Tk_medianState;
	 }//dMag match     

       }//median plane reached
       else{ edm::LogVerbatim(_category)<<"this track has not reached the median plane";}
     }//muPSimhit valid
     else{ edm::LogVerbatim(_category)<<"this track has not reached the either muon/tracker system";}


     edm::LogVerbatim(_category)<<"Id: "<<simIt->trackId()<<"\n"
				<<"type: "<<simIt->type()<<"\n"
				<<"momentum: "<<momentum<<"\n"
				<<"vertex: "<<position<<"\n"
				<<"with ("<<mc_hit_Ltof.size()<<") PSimHit (Ltof) in the tracker\n"
				<<"with ("<<mc_hit_Htof.size()<<") PSimHit (Htof) in the tracker";
     
     uint ips=0;
     for (std::map< double ,TrackPSimHitRef >::iterator hit=mc_hit_Ltof.begin();hit!=mc_hit_Ltof.end();++hit){
       //set the tree content
       DressedHit * simhit = candidate->addSimHit();
       simhit->simhit<<_analyzer.convert(*hit->second);
       TransientTrackingRecHit::ConstRecHitPointer rechit= _analyzer.match(*hit->second);
       if (rechit){simhit->rechit<<(*rechit);}
       if (Tk_medianState.isValid()){
	 const BoundPlane & surface = _analyzer.tracker()->idToDet(DetId(hit->second->detUnitId()))->surface();
	 TrajectoryStateOnSurface atSurface = _analyzer.prop(anyDirection)->propagate(*Tk_medianState.freeState(),surface);
	 if (atSurface.isValid()) {simhit->road<<atSurface;}}

       
       const GeomDet * gd =_analyzer.tracker()->idToDet(DetId(hit->second->detUnitId()));
       edm::LogVerbatim(_category)<<ips++<<"] PSimHit on Id: "<<hit->second->detUnitId() <<"\n"
				  <<"entry position: "<<gd->surface().toGlobal(hit->second->entryPoint())<<"\n"
				  <<"entry momentum: "<<gd->surface().toGlobal(hit->second->momentumAtEntry())<<"\n"
				  <<"at distance: "<<hit->first<<"\n"
				  <<"exit position:  "<<gd->surface().toGlobal(hit->second->exitPoint())<<"\n"
				  <<"time of flight: "<<hit->second->tof()<<"\n"
				  <<"from process type: "<<hit->second->processType()<<"\n"
				  <<"energy loss: "<<hit->second->energyLoss()<<"\n";}
     
     for (std::map< double ,TrackPSimHitRef >::iterator hit=mc_hit_Htof.begin();hit!=mc_hit_Htof.end();++hit){
       //set the tree content
       DressedHit * simhit = candidate->addSimHit();
       simhit->simhit<<_analyzer.convert(*hit->second);
       TransientTrackingRecHit::ConstRecHitPointer rechit= _analyzer.match(*hit->second);
       if (rechit){simhit->rechit<<(*rechit);}
       if (Tk_medianState.isValid()){
	 const BoundPlane & surface = _analyzer.tracker()->idToDet(DetId(hit->second->detUnitId()))->surface();
	 TrajectoryStateOnSurface atSurface = _analyzer.prop(anyDirection)->propagate(*Tk_medianState.freeState(),surface);
	 if (atSurface.isValid()) {simhit->road<<atSurface;}}

       const GeomDet * gd =_analyzer.tracker()->idToDet(DetId(hit->second->detUnitId()));
       edm::LogVerbatim(_category)<<ips++<<"] PSimHit on Id: "<<hit->second->detUnitId() <<"\n"
				  <<"entry position: "<<gd->surface().toGlobal(hit->second->entryPoint())<<"\n"
				  <<"entry momentum: "<<gd->surface().toGlobal(hit->second->momentumAtEntry())<<"\n"
				  <<"at distance: "<<hit->first<<"\n"
				  <<"exit position:  "<<gd->surface().toGlobal(hit->second->exitPoint())<<"\n"
				  <<"time of flight: "<<hit->second->tof()<<"\n"
				  <<"from process type: "<<hit->second->processType()<<"\n"
				  <<"energy loss: "<<hit->second->energyLoss()<<"\n";}



   }//loop over simtracks

   //dump the content of the reconstructed hits
   //-------retreive strip rechits
   edm::Handle<SiStripRecHit2DCollection> rphirecHits; 
   const std::string recHitProducer="siStripMatchedRecHits";
   iEvent.getByLabel(recHitProducer, "rphiRecHit", rphirecHits);
   
   edm::LogVerbatim(_category)<<"there are ("<< rphirecHits->size()<<") r-phi recHits in the event";
   uint irh=0;
   for(SiStripRecHit2DCollection::const_iterator hit_it = rphirecHits->begin(); hit_it!=rphirecHits->end();++hit_it){

     TransientTrackingRecHit::RecHitPointer Trh =_analyzer.TTRHbuilder().build(&(*hit_it));
     if (Trh->isValid()){
       if (Trh->det()){
	 edm::LogVerbatim(_category)<<"       "<<++irh<<") rechits on det: "<<Trh->det()->geographicalId().rawId()
				    <<"               at: "<<Trh->globalPosition();}
       else{
	 edm::LogVerbatim(_category)<<"       "<<++irh<<") rechits on a layer: "
				    <<"               at: "<<Trh->globalPosition();}}
     else{
       if (Trh->det()){
	 edm::LogVerbatim(_category)<<"       "<<++irh<<") rechits on det: "<<Trh->det()->geographicalId().rawId()
				    <<"               invalid...";}
       else{
	 edm::LogVerbatim(_category)<<"       "<<++irh<<") rechits on a layer: "
				    <<"               invalid...";}}
   }//loop over the rechit in the event

   //what are the reconstructerd tracks ?
   for (reco::TrackCollection::const_iterator Tkit = tracks->begin(); Tkit!= tracks->end();++Tkit){
     const reco::Track & muon = *Tkit;

     //build a transient track from this guy
     reco::TransientTrack Tmuon = _analyzer.TTbuilder().build(&muon);

     TrajectoryStateOnSurface medianState = _analyzer.prop(anyDirection)->propagate(*Tmuon.outermostMeasurementState().freeState(),
										    *medianPlane);

     FreeTrajectoryState cIPFTS = _analyzer.transformer().initialFreeState(muon,_analyzer.field().product());
     edm::LogVerbatim(_category)<<"Tk track reconstructed with initial state:\n"<<cIPFTS
				<<"\n outer measurement state:\n"<<Tmuon.outermostMeasurementState()
				<<"\n inner measurement state:\n"<<Tmuon.innermostMeasurementState();
     if (medianState.isValid()){edm::LogVerbatim(_category)<<"\n state at median plane:\n"<<medianState;}
     edm::LogVerbatim(_category)<<"\nwith ("<<muon.recHitsSize()<<") rechits \n";
     
     trackingRecHit_iterator rh_it = muon.recHitsBegin();
     uint irh=0;
     for (; rh_it !=muon.recHitsEnd();++rh_it){
       TransientTrackingRecHit::RecHitPointer Trh =_analyzer.TTRHbuilder().build(&(**rh_it));
       if (Trh->isValid()){
	 if (Trh->det()){
	   edm::LogVerbatim(_category)<<"       "<<++irh<<") rechits on det: "<<Trh->det()->geographicalId().rawId()
					     <<"               at: "<<Trh->globalPosition();}
	 else{
	   edm::LogVerbatim(_category)<<"       "<<++irh<<") rechits on a layer: "
					     <<"               at: "<<Trh->globalPosition();}}
       else{
	 if (Trh->det()){
	   edm::LogVerbatim(_category)<<"       "<<++irh<<") rechits on det: "<<Trh->det()->geographicalId().rawId()
					     <<"               invalid...";}
	 else{
	   edm::LogVerbatim(_category)<<"       "<<++irh<<") rechits on a layer: "
					     <<"               invalid...";}}
     }


   }//loop over the reconstructed tracks
   
   if (fillMe){ _tree->Fill();}

}


// ------------ method called once each job just before starting event loop  ------------
void 
BeamHaloMuonAnalyzer::beginJob(const edm::EventSetup&)
{
  edm::LogInfo(_category)<<"opening a dump root file: "<<_fileName;
  _file = new TFile(_fileName.c_str(),"recreate"); 
  //  _file = TFile::Open(_fileName.c_str(),"recreate");
  _file->cd();
  _dump = new EventDump;
  _tree = new TTree("dump","event dump to monitor beam halo muon tracking");
  _tree->Branch("aDump","EventDump",&_dump);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
BeamHaloMuonAnalyzer::endJob() {
  edm::LogInfo(_category)<<"writing a dump root file with: "<<_tree->GetEntriesFast()<<" entries(fast)";
  _file->cd();
  _tree->Write();
  _file->Close();
}

//define this as a plug-in
DEFINE_FWK_MODULE(BeamHaloMuonAnalyzer);
