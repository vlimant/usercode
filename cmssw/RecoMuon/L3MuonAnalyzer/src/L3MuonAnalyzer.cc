// -*- C++ -*-
//
// Package:    L3MuonAnalyzer
// Class:      L3MuonAnalyzer
// 
/**\class L3MuonAnalyzer L3MuonAnalyzer.cc RecoMuon/L3MuonAnalyzer/src/L3MuonAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jean-Roch Vlimant
//         Created:  Thu Mar 22 17:02:13 CDT 2007
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

#include "RecoMuon/L3MuonAnalyzer/interface/Analyzer.h"

#include "DataFormats/TrackReco/interface/TrackToManyTrajectorySeedAssociation.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedToManyTrackCandidateAssociation.h"
#include "DataFormats/TrackReco/interface/TrackTrackCandidateAssociation.h"
#include "DataFormats/TrackReco/interface/TrackToManyTrackAssociation.h"

#include "RecoMuon/L3MuonAnalyzer/interface/DumpClass.h"

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

//
// class decleration
//

class L3MuonAnalyzer : public edm::EDAnalyzer {
   public:
      explicit L3MuonAnalyzer(const edm::ParameterSet&);
      ~L3MuonAnalyzer();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
  edm::InputTag _standAloneLabel;
  edm::InputTag _seedGeneratorLabel;
  edm::InputTag _trackCandMakerLabel;
  edm::InputTag _trackProducerLabel;
  edm::InputTag _associationMakerLabel;
  std::string _category;

  IntrusiveAnalyzer _analyzer;

  bool _skipPixel;

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
L3MuonAnalyzer::L3MuonAnalyzer(const edm::ParameterSet& iConfig)
{
  _category = "L3MuonAnalyzer";
   //now do what ever initialization is needed

  _standAloneLabel = iConfig.getParameter<edm::InputTag>("standAloneLabel");
  _seedGeneratorLabel = iConfig.getParameter<edm::InputTag>("seedGeneratorLabel");
  _trackCandMakerLabel = iConfig.getParameter<edm::InputTag>("trackCandMakerLabel");
  _trackProducerLabel = iConfig.getParameter<edm::InputTag>("trackProducerLabel");
  _associationMakerLabel = iConfig.getParameter<edm::InputTag>("associationMakerLabel");
  
  
  //create a roottree
  _fileName=iConfig.getParameter<std::string>("@module_label")+".root";
  
  _skipPixel=!iConfig.getParameter<bool>("considerPixel");
}


L3MuonAnalyzer::~L3MuonAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)


}


//
// member functions
//

// ------------ method called to for each event  ------------
void
L3MuonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   _analyzer.setUp(iEvent,iSetup);

   _dump->clean();
   bool fillMe=false;
   _dump->run = (int)iEvent.id().run();
   _dump->evt = (int)iEvent.id().event();

   edm::Handle<reco::TrackCollection> STAtrack;
   edm::Handle<reco::TrackToManyTrajectorySeedAssociationCollection> STAtrack_to_Tkseed;
   edm::Handle<reco::TrajectorySeedToManyTrackCandidateAssociationCollection> Tkseed_to_Tkcand;
   edm::Handle<reco::TrackTrackCandidateAssociationCollection> Tktrack_to_Tkcand;
   edm::Handle<reco::TrackToManyTrackAssociationCollection> STAtrack_to_Tktrack;

   //retreive the inputs
   iEvent.getByLabel(_standAloneLabel,STAtrack);
   iEvent.getByLabel(_seedGeneratorLabel,STAtrack_to_Tkseed);
   iEvent.getByLabel(_trackCandMakerLabel,Tkseed_to_Tkcand);
   iEvent.getByLabel(_trackProducerLabel,Tktrack_to_Tkcand);
   iEvent.getByLabel(_associationMakerLabel,STAtrack_to_Tktrack);

   edm::LogInfo(_category)<<" ++++++++++++++++++++++++++REPORT++++++++++++++++++++++++++++++++ ";
   edm::LogVerbatim(_category)<<" ============ STA to Tkt ================= ";
   for (reco::TrackToManyTrackAssociationCollection::const_iterator STA_to_Tk_it=  STAtrack_to_Tktrack->begin(); STA_to_Tk_it!=STAtrack_to_Tktrack->end() ;++STA_to_Tk_it) {
     const reco::Track & muon = *STA_to_Tk_it->key;

     //build a transient track from this guy
     reco::TransientTrack Tmuon = _analyzer.TTbuilder().build(&muon);

     FreeTrajectoryState cIPFTS = _analyzer.transformer().initialFreeState(muon,_analyzer.field().product());
     edm::LogVerbatim(_category)<<"("<<STA_to_Tk_it->val.size()<<") Tk track associated with this STA track\n"
       				       <<"STA track: \n"<<cIPFTS<<"------";

     uint it=0;
     for (edm::RefVector<reco::TrackCollection>::const_iterator Tkit = STA_to_Tk_it->val.begin(); Tkit!=STA_to_Tk_it->val.end(); ++Tkit){
       const reco::Track & tk = **Tkit;
       FreeTrajectoryState tcIPFTS = _analyzer.transformer().initialFreeState(tk,_analyzer.field().product());
       edm::LogVerbatim(_category)<<++it<<"]   Tk track with ("<< tk.recHitsSize()<<") rechits \n"<<tcIPFTS<<"------";
       trackingRecHit_iterator rh_it = tk.recHitsBegin();
       uint irh=0;
       for (; rh_it !=tk.recHitsEnd();++rh_it){
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
     }
     //is it matched to a sim track ?
     const SimTrack * match = &(*_analyzer.match(Tmuon));//spatial matching
     if (match){
       //only match to muon
       FreeTrajectoryState matchFts = _analyzer.convert(match);
       edm::LogVerbatim(_category)<<"The STA track is matched to simulated muon:\n"
					 <<matchFts;
       //do something more about it.
     }
     else{
       edm::LogVerbatim(_category)<<"this muons's a not a muon";
       //tell me what is this STA muon
     }

     _analyzer.setDebug(false);
     //gather the rechits in the road 
     double chi2 =40;
     IntrusiveAnalyzer::RecHitWithStateContainer roadcontent = _analyzer.roadContent(Tmuon,chi2,true,true);//with pixel and strips

     edm::LogVerbatim(_category)<<roadcontent.size()<<" rechits within chi2<"<<chi2<<" of the muon STA track";
     uint rh=0;
     for (IntrusiveAnalyzer::RecHitWithStateContainer::const_iterator iTThit = roadcontent.begin();iTThit != roadcontent.end();iTThit++){
       if (iTThit->first->isValid()){
	 edm::LogVerbatim(_category)<<++rh<<"] rechit on det: "<<iTThit->first->det()->geographicalId().rawId()
					   <<"    at: "<<iTThit->first->globalPosition()<<"\n"
					   <<"    with state position: "<<iTThit->second.globalPosition()
					   <<"    with state momentum: "<<iTThit->second.globalMomentum();}
       else{
	 edm::LogVerbatim(_category)<<++rh<<"] rechit on det: "<<iTThit->first->det()->geographicalId().rawId()
					   <<"    invalid...\n"
					   <<"    with state position: "<<iTThit->second.globalPosition()
					   <<"    with state momentum: "<<iTThit->second.globalMomentum();}
     }
     _analyzer.setDebug(false);

     edm::LogVerbatim(_category)<<"\n------------------\n";

     



     /*     //fill the tree according to those information
	    RoadCandidate * candidate = _dump->add();
	    candidate->hasSeed=true;
	    LogDebug(_category)<<O(Tmuon.innermostMeasurementState());
	    candidate->seed<<Tmuon.innermostMeasurementState();
	    candidate->Nrechits=Tmuon.recHitsSize();
	    if (match){
	    candidate->hasSim=true;
	    LogDebug(_category)<<O(_analyzer.convert(match));
	    candidate->sim<<_analyzer.convert(match);}
	    candidate->hasCpca=true;
	    LogDebug(_category)<<O(Tmuon.initialFreeState());
	    candidate->Cpca<<Tmuon.initialFreeState();
	    candidate->hasOutIn=false;
	    candidate->hasOutOut=false;

	    //fill the road hit information
	    for (IntrusiveAnalyzer::RecHitWithStateContainer::const_iterator iTThit = roadcontent.begin();iTThit != roadcontent.end();iTThit++){
	    if (!iTThit->first->isValid()){edm::LogVerbatim(_category)<<"TTRH not valid !"; continue;}
	    DressedHit * hit = candidate->addRecHit();       
	    std::pair<IntrusiveAnalyzer::RecHit_match_PSimHit, TrajectoryStateOnSurface> matchHit;
	    if (match)
	    matchHit = _analyzer.match(*(iTThit->first),false,match->trackId());//specific trakc ID
	    else{
	    matchHit = _analyzer.match(*(iTThit->first));//no trackId. match code will be funny
	    matchHit.first=IntrusiveAnalyzer::OTHERTRACK;}
	    hit->matchcode = matchHit.first;
	    if (!matchHit.second.isValid()){
	    if( !iTThit->first->isValid()){edm::LogError(_category)<<"rechit is really bad"; continue;}
	    else {hit->simhit<<(*iTThit->first);}
	    }
	    else {
	    LogDebug(_category)<<O(matchHit.second);
	    hit->simhit<<matchHit.second;}
	    if (iTThit->first->isValid()){
	    LogDebug(_category)<<"(*iTThit->first)";
	    hit->rechit<<(*iTThit->first);}
	    LogDebug(_category)<<O(iTThit->second);
	    hit->road<<iTThit->second; //the propagated state
	    }


	    if (match){
	    //fill the simtrack hit information
	    std::map< double, const PSimHit > mccontent = _analyzer.MapCount(match->trackId(),true,false,true);//only ionisation shown
	    for (std::map< double, const PSimHit >::iterator mc_hit=mccontent.begin(); mc_hit!=mccontent.end();++mc_hit){
	    DressedHit * simhit = candidate->addSimHit();
	    LogDebug(_category)<<O(_analyzer.convert(mc_hit->second));
	    simhit->simhit<<_analyzer.convert(mc_hit->second);
	    TransientTrackingRecHit::ConstRecHitPointer rechit= _analyzer.match(mc_hit->second);
	    if (rechit){
	    LogDebug(_category)<<"(*rechit)";
	    simhit->rechit<<(*rechit);}
	    }
	    }

	    //select a best track among those associated
	    double minChi2_ndof = 100000;
	    reco::TrackRef best;
	    for (edm::RefVector<reco::TrackCollection>::const_iterator Tkit = STA_to_Tk_it->val.begin(); Tkit!=STA_to_Tk_it->val.end(); ++Tkit){
	    const reco::Track & tk = **Tkit;
	    //select the best 
	    if (tk.ndof()!=0){if (tk.chi2()/tk.ndof() < minChi2_ndof){minChi2_ndof = tk.chi2()/tk.ndof();best= *Tkit;}}
	    }
	    if (best.isNonnull()) { edm::LogVerbatim(_category)<<"there is a best match with chi2/ndof="<<minChi2_ndof<<" at index: "<<best.index();}
	    
	    //fill the information for all the tracks assocatied to that guy
	    for (edm::RefVector<reco::TrackCollection>::const_iterator Tkit = STA_to_Tk_it->val.begin(); Tkit!=STA_to_Tk_it->val.end(); ++Tkit){
	    const reco::Track & tk = **Tkit;
	    if (*Tkit==best) continue; //skip the best in the other candidate collection
	    reco::TransientTrack Ttk = _analyzer.TTbuilder().build(&tk);
	    Candidate * other_candidate = candidate->addCandidate();
	    LogDebug(_category)<<O(Ttk.initialFreeState());
	    other_candidate->updated<<Ttk.initialFreeState();
	    trackingRecHit_iterator rh = Ttk.recHitsBegin();
	    for (; rh !=Ttk.recHitsEnd();++rh){
	    if (!(*rh)->isValid()) {edm::LogVerbatim(_category)<<"RH not valid ("<<(*rh)->geographicalId().rawId()<<")"; continue;}
	    std::pair<IntrusiveAnalyzer::RecHit_match_PSimHit, TrajectoryStateOnSurface> matchHit;
	    if (match) 
	    matchHit = _analyzer.match(**rh,false,match->trackId());//specific trakc ID 
	    else{
	    matchHit = _analyzer.match(**rh);//no trackId. match code will be funny  
	    matchHit.first=IntrusiveAnalyzer::OTHERTRACK;}
	    DressedHit * hit = other_candidate->add();
	    hit->matchcode = matchHit.first;
	    TransientTrackingRecHit::RecHitPointer Trh = _analyzer.TTRHbuilder().build(&(**rh));
	    if (!matchHit.second.isValid()){
	    if( !Trh->isValid()){edm::LogError(_category)<<"rechit is really bad"; continue;}
	    else {hit->simhit<<(*Trh);}
	    }
	    else    {
	    hit->simhit<<matchHit.second; }
	    
	    if (Trh->isValid()){
	    LogDebug(_category)<<"(*Trh)";
	    hit->rechit<<(*Trh);}
	    }
	    }
	    
	    //fill the best track information
	    if (best.isNonnull()) {
	    reco::TransientTrack Ttk = _analyzer.TTbuilder().build(&(*best));
	    candidate->hasBest=true;
	    LogDebug(_category)<<O(Ttk.initialFreeState());
	    candidate->best.updated<<Ttk.initialFreeState();
	    trackingRecHit_iterator rh = Ttk.recHitsBegin(); 
	    for (; rh !=Ttk.recHitsEnd();++rh){ 
	    if (!(*rh)->isValid()) {edm::LogVerbatim(_category)<<"RH not valid ("<<(*rh)->geographicalId().rawId()<<")"; continue;}
	    std::pair<IntrusiveAnalyzer::RecHit_match_PSimHit, TrajectoryStateOnSurface> matchHit;
	    if (match) 
	    matchHit = _analyzer.match(**rh,false,match->trackId());//specific trakc ID
	    else{
	    matchHit = _analyzer.match(**rh);//no trackId. match code will be funny 
	    matchHit.first=IntrusiveAnalyzer::OTHERTRACK;}
	    
	    DressedHit * hit = candidate->best.add();
	    hit->matchcode = matchHit.first;
	    TransientTrackingRecHit::RecHitPointer Trh = _analyzer.TTRHbuilder().build(&(**rh));
	    if (!matchHit.second.isValid()){
	    if( !Trh->isValid()){edm::LogError(_category)<<"rechit is really bad"; continue;}
	    else {hit->simhit<<(*Trh);}
	    }
	    else	 {
	    hit->simhit<<matchHit.second; }
	    
	    if (Trh->isValid()){
	    LogDebug(_category)<<"(*Trh)";
	 hit->rechit<<(*Trh);}
	 }
	 }
     */

   }//loop over STA to Tkt

   edm::LogVerbatim(_category)<<" ============ Tkt to STA ================= ";
   //loop over the track in the event


   edm::LogVerbatim(_category)<<" ============ STA to seed ================= ";
   //loop over the STA track in the event
   
   for (reco::TrackToManyTrajectorySeedAssociationCollection::const_iterator STA_to_seed=STAtrack_to_Tkseed->begin(); STA_to_seed!=STAtrack_to_Tkseed->end();++STA_to_seed){
     const reco::Track & muon = *STA_to_seed->key;     

     edm::RefVector<TrajectorySeedCollection> seeds = STA_to_seed->val;
     FreeTrajectoryState cIPFTS = _analyzer.transformer().initialFreeState(muon,_analyzer.field().product());
     edm::LogVerbatim(_category)<<"("<<seeds.size()<<") seeds associated with this STA track\n"
				       <<"STA track: \n"<<cIPFTS<<"------";
     uint is=0;
     for (edm::RefVector<TrajectorySeedCollection>::const_iterator seedit = seeds.begin(); seedit!=seeds.end();++seedit){
       edm::LogVerbatim(_category)<<++is<<"] seed with "<<(**seedit).nHits()<<" rechits\n"
					 <<"starting at:\n"<<_analyzer.transformer().transientState((**seedit).startingState(),
												    &_analyzer.tracker()->idToDet(DetId((**seedit).startingState().detId()))->surface(),
												    _analyzer.field().product())
					 <<"on detId: "<<(**seedit).startingState().detId()<<"\n"
					 <<"going: "<<(**seedit).direction();
       reco::TrajectorySeedToManyTrackCandidateAssociationCollection::const_iterator seedit_in_assoc=Tkseed_to_Tkcand->find(*seedit);
       if (seedit_in_assoc!=Tkseed_to_Tkcand->end())
	 {edm::LogVerbatim(_category)<<"this seed has been pursued into ("<< seedit_in_assoc->val.size()<<") track candidates";}
       else
	 {edm::LogVerbatim(_category)<<"this seed was NOT persued.";}

       
     }
   }

   edm::LogVerbatim(_category)<<" ============ SIM to track ================= ";
   //loop over the muon simtracks and find out if they have been reconstructed
   //   std::vector<const SimTrack*> mc_muons=_analyzer.particles(13,0);//muon, whatever the charge
   //   for(std::vector<const SimTrack*> ::iterator mc_it=mc_muons.begin(); mc_it!=mc_muons.end();++mc_it)
   std::vector<SimTrackRef> mc_muons =_analyzer.particles(13,0);
   for (std::vector<SimTrackRef>::iterator mcref_it = mc_muons.begin(); mcref_it!=mc_muons.end();++mcref_it)
     {
       const SimTrack * mc_it=&(**mcref_it);
       edm::LogVerbatim(_category)<<"a simulated muon with state:\n"<<_analyzer.convert(mc_it);

       //find out if any STA track is matched to this simtrack
       reco::TrackRef muon;
       bool matched=false;
       for (reco::TrackToManyTrackAssociationCollection::const_iterator STA_to_Tk_it=  STAtrack_to_Tktrack->begin(); STA_to_Tk_it!=STAtrack_to_Tktrack->end() ;++STA_to_Tk_it){
	 muon = STA_to_Tk_it->key;
	 //build a transient track from this guy
	 reco::TransientTrack Tmuon = _analyzer.TTbuilder().build(&(*muon));
	 const SimTrack * Mmatch = &(*_analyzer.match(Tmuon,0.2));//position dR=0.2
	 GlobalPoint mom(muon->momentum().X(),muon->momentum().Y(),muon->momentum().Z());
	 const SimTrack * Pmatch = &(*_analyzer.match(mom,NULL,0.2));//momentum dR=0.2
	 if (Mmatch){
	   if (Mmatch->momentum() == (mc_it)->momentum()) {matched=true;break;} //that's a proper one
	 }
	 if (Pmatch){
	   if (Pmatch->momentum() == (mc_it)->momentum()) {matched=true;break;} //that's a proper one
	 }
       }

       if (matched){
	 //there is a STA track associated to this simulated muon
	 edm::LogVerbatim(_category)<<" is matched to a STA track:\n"
					   <<_analyzer.transformer().initialFreeState(*muon,_analyzer.field().product());
       }
       else{
	 edm::LogVerbatim(_category)<<" is not matched to a STA track";
	 //simulated muon never reconstructed as a STA track
       }
       edm::LogVerbatim(_category)<<"\n------\n";

       //find out if any Tk track is mathed to this simtrack
       reco::TrackRef track;
       matched=false;
       for (reco::TrackTrackCandidateAssociationCollection::const_iterator Tk_to_Tkc_it = Tktrack_to_Tkcand->begin(); Tk_to_Tkc_it!=Tktrack_to_Tkcand->end(); ++Tk_to_Tkc_it) {
	 track = Tk_to_Tkc_it->key;
	 //	 reco::TransientTrack Ttrack = _analyzer.TTbuilder().build(&(*track));
	 //	 const SimTrack * match = _analyzer.match(Ttrack,0.2); //position dR=0.2
	 GlobalPoint mom(track->momentum().X(),track->momentum().Y(),track->momentum().Z());
	 SimTrackRef matchRef = _analyzer.match(mom,NULL,0.2);  //momentum dR=0.2
	 if (matchRef.isNonnull()){
	   const SimTrack * match = &(*_analyzer.match(mom,NULL,0.2)); //momentum dR=0.2
	   if (!match) continue;
	   if (match->momentum() == (mc_it)->momentum()) {matched=true;break;} //that's a proper one
	 }
       }
       if (matched){
	 //there is indeed a tracker track that is matched to that muon
	 edm::LogVerbatim(_category)<<" is matched to a Tk track:\n"
					   <<_analyzer.transformer().initialFreeState(*track,_analyzer.field().product());
       }
       else{
	 edm::LogVerbatim(_category)<<" is not matched to a Tk track";
	 //nope this simulated track has not been recsontructed as a track
       }
     }


   edm::LogVerbatim(_category)<<" ============ Fill the tree ================= ";
   
   //a list of pair sim/reco whether there's a match or not
   std::vector<std::pair<SimTrackRef,reco::TrackRef> > matchedList = _analyzer.match(mc_muons, STAtrack);

   for (std::vector<std::pair<SimTrackRef,reco::TrackRef> > ::iterator mit = matchedList.begin();mit!=matchedList.end();++mit){
     //create a candidate
     RoadCandidate * candidate = _dump->add();
     
     //simulated track ?
     const SimTrack * match=NULL;
     if (mit->first.isNonnull()){
       match= &(*mit->first);
       candidate->hasSim=true;
       LogDebug(_category)<<O(_analyzer.convert(match));
       candidate->sim<<_analyzer.convert(match);
       
       //fill the simtrack hit information
       //only ionisation rechit
       std::map< double, const PSimHit > mccontent = _analyzer.MapCount(match->trackId(),true,false,true,"LowTof",_skipPixel);
       for (std::map< double, const PSimHit >::iterator mc_hit=mccontent.begin(); mc_hit!=mccontent.end();++mc_hit){
	 DressedHit * simhit = candidate->addSimHit();
	 LogDebug(_category)<<O(_analyzer.convert(mc_hit->second));
	 simhit->simhit<<_analyzer.convert(mc_hit->second);
	 LogDebug(_category)<<"matching PSimHit to a rechit (use ransientTrackingRecHit)";
	 TransientTrackingRecHit::ConstRecHitPointer rechit= _analyzer.match(mc_hit->second);
	 LogDebug(_category)<<"match done. let's look at validity";
	 if (rechit){
	   LogDebug(_category)<<"(*rechit)";
	   simhit->rechit<<(*rechit);}
	 //add the road information here
	 if (mit->second.isNonnull()){
	   reco::TransientTrack Tmuon = _analyzer.TTbuilder().build(&(*mit->second));
	   FreeTrajectoryState cPCA = Tmuon.initialFreeState();
	   if (cPCA.position().mag()!=0 && cPCA.momentum().mag()!=0) {
	     const BoundPlane & surface = _analyzer.tracker()->idToDet(DetId(mc_hit->second.detUnitId()))->surface();
	     TrajectoryStateOnSurface tsos = _analyzer.prop(alongMomentum)->propagate(cPCA,surface);
	     if (tsos.isValid()){
	       simhit->road<<tsos;}
	   }//free state is valid
	 }//cpca state available
       }//loop over simhits
     }//there is a simulated track

     //STA track ?
     if (mit->second.isNonnull()){
       reco::TransientTrack Tmuon = _analyzer.TTbuilder().build(&(*mit->second));
       
       candidate->hasSeed=true;
       LogDebug(_category)<<O(Tmuon.innermostMeasurementState());
       candidate->seed<<Tmuon.innermostMeasurementState();
       candidate->Nrechits=Tmuon.recHitsSize();
	    
       FreeTrajectoryState cPCA = Tmuon.initialFreeState();
       if (cPCA.position().mag()!=0 || cPCA.momentum().mag()!=0){
 	 //otherwise, it has failed to get an initialFreeState, most likely the track does not propagate back to IP...

	 candidate->hasCpca=true;
	 LogDebug(_category)<<O(Tmuon.initialFreeState());
	 candidate->Cpca<<Tmuon.initialFreeState();
	 candidate->hasOutIn=false;
	 candidate->hasOutOut=false;
	 
	 //       edm::LogInfo(_category)<<O(candidate->hasCpca)<<O(candidate->hasSeed)<<O(_dump->run)<<O(_dump->evt);
	 
	 double chi2 =40;
	 IntrusiveAnalyzer::RecHitWithStateContainer roadcontent = _analyzer.roadContent(Tmuon,chi2,true,true);//with pixel and strips
	 //fill the road hit information
	 for (IntrusiveAnalyzer::RecHitWithStateContainer::const_iterator iTThit = roadcontent.begin();iTThit != roadcontent.end();iTThit++){
	   if (!iTThit->first->isValid()){edm::LogVerbatim(_category)<<"TTRH not valid !"; continue;}
	   DressedHit * hit = candidate->addRecHit();       
	   std::pair<IntrusiveAnalyzer::RecHit_match_PSimHit, TrajectoryStateOnSurface> matchHit;
	   if (match)
	     matchHit = _analyzer.match(*(iTThit->first),false,match->trackId());//specific trakc ID
	   else{
	     matchHit = _analyzer.match(*(iTThit->first));//no trackId. match code will be funny
	     matchHit.first=IntrusiveAnalyzer::OTHERTRACK;}
	   hit->matchcode = matchHit.first;
	   if (!matchHit.second.isValid()){
	     if( !iTThit->first->isValid()){edm::LogError(_category)<<"rechit is really bad"; continue;}
	     else {hit->simhit<<(*iTThit->first);}
	   }
	   else {
	     LogDebug(_category)<<O(matchHit.second);
	     hit->simhit<<matchHit.second;}
	   if (iTThit->first->isValid()){
	     LogDebug(_category)<<"(*iTThit->first)";
	     hit->rechit<<(*iTThit->first);}
	   LogDebug(_category)<<O(iTThit->second);
	   hit->road<<iTThit->second; //the propagated state
	 }
       }//not propagatable to IP

       //is the STA into Tk seed
       reco::TrackToManyTrajectorySeedAssociationCollection::const_iterator STA_to_seed=STAtrack_to_Tkseed->find(mit->second);
       if (STA_to_seed!=STAtrack_to_Tkseed->end()){
	 edm::RefVector<TrajectorySeedCollection> seeds = STA_to_seed->val;
	 
	 for (edm::RefVector<TrajectorySeedCollection>::const_iterator seedit = seeds.begin(); seedit!=seeds.end();++seedit){
	   Candidate * seed = candidate->addSeed();
	   
	   seed->updated<<_analyzer.transformer().transientState((**seedit).startingState(),
								&_analyzer.tracker()->idToDet(DetId((**seedit).startingState().detId()))->surface(),
								_analyzer.field().product());
	 }//loop over the associated seeds
       }// STA is into a vector<seed>



       //is the STA into Tk tracks ?
       reco::TrackToManyTrackAssociationCollection::const_iterator STA_to_Tk_it=STAtrack_to_Tktrack->find(mit->second);
       if (STA_to_Tk_it!=STAtrack_to_Tktrack->end()){
	     
	     //select a best track among those associated
	     double minChi2_ndof = 100000;
	     reco::TrackRef best;
	     for (edm::RefVector<reco::TrackCollection>::const_iterator Tkit = STA_to_Tk_it->val.begin(); Tkit!=STA_to_Tk_it->val.end(); ++Tkit){
	       const reco::Track & tk = **Tkit;
	       //select the best 
	       if (tk.ndof()!=0){if (tk.chi2()/tk.ndof() < minChi2_ndof){minChi2_ndof = tk.chi2()/tk.ndof();best= *Tkit;}}
	     }
	     if (best.isNonnull()) { edm::LogVerbatim(_category)<<"there is a best match with chi2/ndof="<<minChi2_ndof<<" at index: "<<best.index();}
	     
	     //fill the information for all the tracks assocatied to that guy
	     for (edm::RefVector<reco::TrackCollection>::const_iterator Tkit = STA_to_Tk_it->val.begin(); Tkit!=STA_to_Tk_it->val.end(); ++Tkit){
	       const reco::Track & tk = **Tkit;
	       if (*Tkit==best) continue; //skip the best in the other candidate collection
	       reco::TransientTrack Ttk = _analyzer.TTbuilder().build(&tk);
	       Candidate * other_candidate = candidate->addCandidate();
	       LogDebug(_category)<<O(Ttk.initialFreeState());
	       other_candidate->updated<<Ttk.initialFreeState();
	       trackingRecHit_iterator rh = Ttk.recHitsBegin();
	       for (; rh !=Ttk.recHitsEnd();++rh){
		 if (!(*rh)->isValid()) {edm::LogVerbatim(_category)<<"RH not valid ("<<(*rh)->geographicalId().rawId()<<")"; continue;}
		 std::pair<IntrusiveAnalyzer::RecHit_match_PSimHit, TrajectoryStateOnSurface> matchHit;
		 if (match) 
		   matchHit = _analyzer.match(**rh,false,match->trackId());//specific trakc ID 
		 else{
		   matchHit = _analyzer.match(**rh);//no trackId. match code will be funny  
		   matchHit.first=IntrusiveAnalyzer::OTHERTRACK;}
		 DressedHit * hit = other_candidate->add();
		 hit->matchcode = matchHit.first;
		 TransientTrackingRecHit::RecHitPointer Trh = _analyzer.TTRHbuilder().build(&(**rh));
		 if (!matchHit.second.isValid()){
		   if( !Trh->isValid()){edm::LogError(_category)<<"rechit is really bad"; continue;}
		   else {hit->simhit<<(*Trh);}
		 }
		 else    {
		   hit->simhit<<matchHit.second; }
		 
		 if (Trh->isValid()){
		   LogDebug(_category)<<"(*Trh)";
		   hit->rechit<<(*Trh);}
	       }
	     }

	     //fill the best track information
	     if (best.isNonnull()) {
	       reco::TransientTrack Ttk = _analyzer.TTbuilder().build(&(*best));
	       candidate->hasBest=true;
	       LogDebug(_category)<<O(Ttk.initialFreeState());
	       candidate->best.updated<<Ttk.initialFreeState();
	       trackingRecHit_iterator rh = Ttk.recHitsBegin(); 
	       for (; rh !=Ttk.recHitsEnd();++rh){ 
		 if (!(*rh)->isValid()) {edm::LogVerbatim(_category)<<"RH not valid ("<<(*rh)->geographicalId().rawId()<<")"; continue;}
		 std::pair<IntrusiveAnalyzer::RecHit_match_PSimHit, TrajectoryStateOnSurface> matchHit;
		 if (match) 
		   matchHit = _analyzer.match(**rh,false,match->trackId());//specific trakc ID
		 else{
		   matchHit = _analyzer.match(**rh);//no trackId. match code will be funny 
		   matchHit.first=IntrusiveAnalyzer::OTHERTRACK;}
		 
		 DressedHit * hit = candidate->best.add();
		 hit->matchcode = matchHit.first;
		 TransientTrackingRecHit::RecHitPointer Trh = _analyzer.TTRHbuilder().build(&(**rh));
		 if (!matchHit.second.isValid()){
		   if( !Trh->isValid()){edm::LogError(_category)<<"rechit is really bad"; continue;}
		   else {hit->simhit<<(*Trh);}
		 }
		 else	 {
		   hit->simhit<<matchHit.second; }
		 
		 if (Trh->isValid()){
		   LogDebug(_category)<<"(*Trh)";
		   hit->rechit<<(*Trh);}
	       }
	     }
	   }//STA is into Tk track
     }//is STA recoed
     fillMe=true;
   }//loop over the sim/reco matches

   edm::LogVerbatim(_category)<<" ============ done ================= ";



   if (fillMe){
     edm::LogInfo(_category)<<"filling an entry";
     _tree->Fill();}
}





// ------------ method called once each job just before starting event loop  ------------
void 
L3MuonAnalyzer::beginJob(const edm::EventSetup&)
{
  edm::LogInfo(_category)<<"opening a dump root file: "<<_fileName;
  _file = TFile::Open(_fileName.c_str(),"recreate");
  _file->cd();
  _dump = new EventDump;
  _tree = new TTree("dump","event dump to monitor muon road search process");
  _tree->Branch("aDump","EventDump",&_dump);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
L3MuonAnalyzer::endJob() {

  edm::LogInfo(_category)<<"writing a dump root file with: "<<_tree->GetEntriesFast()<<" entries(fast)";
  _file->cd();
  _tree->Write();
  _file->Close();
  
  //  delete _dump;
  //  delete _tree;
}

//define this as a plug-in
DEFINE_FWK_MODULE(L3MuonAnalyzer);
