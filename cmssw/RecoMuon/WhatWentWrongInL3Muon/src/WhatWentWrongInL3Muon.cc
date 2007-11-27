#include "RecoMuon/WhatWentWrongInL3Muon/interface/WhatWentWrongInL3Muon.h"

WhatWentWrongInL3Muon::WhatWentWrongInL3Muon(const edm::ParameterSet& iConfig)

{
  category ="WhatWentWrongInL3Muon";

   //now do what ever initialization is needed
  L2Label = iConfig.getParameter<edm::InputTag> ("L2Label");
  L3Labels = iConfig.getParameter<std::vector<edm::InputTag> >("L3Labels");
  trackingParticleLabel = iConfig.getParameter<edm::InputTag> ("trackingParticleLabel");

  theDQM = edm::Service<DaqMonitorBEInterface>().operator->();
  theRootFileName = iConfig.getParameter<std::string>("rootFileName");
  //  theRootFileName = iConfig.getParameter<std::string>("@module_label")+".root";

  theAssocName = iConfig.getParameter<std::string>("associator");
  thePropagatorName = iConfig.getParameter<std::string>("propagator");

  theTellMe = iConfig.getParameter<bool>("tellMe");
}


WhatWentWrongInL3Muon::~WhatWentWrongInL3Muon()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


void WhatWentWrongInL3Muon::one(nodePlotter & node1, reco::TrackRef & refL2, FreeTrajectoryState & l2State){
  node1.element("nRecHit")->Fill(refL2->recHitsSize());
  node1.element("nRecHit_2D")->Fill(l2State.momentum().perp(), l2State.momentum().eta(), refL2->recHitsSize()); 
  node1.element("px")->Fill(l2State.momentum().x());
  node1.element("py")->Fill(l2State.momentum().y());
  node1.element("pz")->Fill(l2State.momentum().z());
  node1.element("pt")->Fill(l2State.momentum().perp());
  node1.element("p")->Fill(l2State.momentum().mag());
  node1.element("m_eta")->Fill(l2State.momentum().eta());
  node1.element("m_phi")->Fill(l2State.momentum().phi());
  node1.element("x")->Fill(l2State.position().x());
  node1.element("y")->Fill(l2State.position().y());
  node1.element("z")->Fill(l2State.position().z());
  node1.element("v_eta")->Fill(l2State.position().eta());
  node1.element("v_phi")->Fill(l2State.position().phi());
  node1.element("chi2")->Fill(refL2->chi2());
  node1.element("Prob")->Fill(TMath::Prob(refL2->chi2(),(int) refL2->ndof()));
}


void WhatWentWrongInL3Muon::oneLocal(nodePlotter & node1, reco::RecoToSimCollection & recSimColl, reco::TrackRef & refL2, FreeTrajectoryState & l2State){
  //find out about the simulated hits
  //is there a sim match
  reco::RecoToSimCollection::const_iterator ass= recSimColl.find(refL2);
  if (ass!=recSimColl.end()){
    //    const reco::TrackRef & track = RtSit->key;
    const std::vector<std::pair<TrackingParticleRef,double> > & tp = ass->val;
    if (tp.size()==0){ edm::LogError(category)<<"track showing up in association map, even though there are no associations.";}
    else{
      std::vector<std::pair<TrackingParticleRef,double> >::const_iterator vector_iterator = tp.begin();
      for (;vector_iterator!=tp.end();++vector_iterator){
	const std::pair<TrackingParticleRef,double> & matching_pair = *vector_iterator;
	//loop simhits
	std::vector<PSimHit>::const_iterator pSimIt= matching_pair.first->pSimHit_begin();
	node1.element("matchId")->Fill(matching_pair.first->pdgId());
	//consider the residual/pull only for muon
	if (fabs(matching_pair.first->pdgId())!= 13) continue;
	
	for (;pSimIt!=matching_pair.first->pSimHit_end();++pSimIt){
	  //get the surface 
	  const GeomDet * gdet = theTrackerGeometry->idToDet(pSimIt->detUnitId());
	  if (!gdet){edm::LogError(category)<<pSimIt->detUnitId()<<" not valid geomdet."; continue;}
	  const BoundPlane & plane = gdet->surface();
	  
	  //propagate to the surface
	  TrajectoryStateOnSurface l2pState = thePropagator->propagate(l2State, plane);
	  if (!l2pState.isValid()){edm::LogError(category)<<"L2 tracker does not propagate to surface of: "<<pSimIt->detUnitId(); continue;}
	  
	  //you can estimate the distance if you want
	  LocalPoint simpoint = pSimIt->localPosition();
	  LocalPoint l2point = l2pState.localPosition();
	  LocalVector D= l2point - simpoint;
	  LocalError error = l2pState.localError().positionError();

	  GlobalPoint l2glbPoint = plane.toGlobal(l2point);
	  GlobalPoint simglbPoint = plane.toGlobal(simpoint);
	  
	  node1.element("lDx")->Fill(D.x());
	  node1.element("lDx_rho")->Fill(D.x(),l2glbPoint.mag());
	  node1.element("lDx_r")->Fill(D.x(),l2glbPoint.perp());
	  node1.element("lDy")->Fill(D.y());
	  node1.element("lDy_rho")->Fill(D.y(),l2glbPoint.mag());
	  node1.element("lDy_r")->Fill(D.y(),l2glbPoint.perp());
	  node1.element("lDxDy")->Fill(D.x(),D.y());
	  node1.element("lD")->Fill(D.perp());
	  node1.element("lDx_pull")->Fill(D.x()/sqrt(error.xx()));
	  node1.element("lDx_pull_rho")->Fill(D.x()/sqrt(error.xx()),l2glbPoint.mag());
	  node1.element("lDx_pull_r")->Fill(D.x()/sqrt(error.xx()),l2glbPoint.perp());
	  node1.element("lDy_pull")->Fill(D.y()/sqrt(error.yy()));
	  node1.element("lDy_pull_rho")->Fill(D.y()/sqrt(error.yy()),l2glbPoint.mag());
	  node1.element("lDy_pull_r")->Fill(D.y()/sqrt(error.yy()),l2glbPoint.perp());
	  
	}//loop over association to tracking particle
      }//loop over simulated hit
    }//no association to tracking particle
  }//no entry for this track in the associatin map
}

//
// member functions
//

// ------------ method called to for each event  ------------
void
WhatWentWrongInL3Muon::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   //get the mag field
   edm::ESHandle<MagneticField> field;
   iSetup.get<IdealMagneticFieldRecord>().get(field);
   
   //get the associator and tracking particle
   iSetup.get<TrackAssociatorRecord>().get(theAssocName,theAssociator);
   edm::Handle<TrackingParticleCollection> TPtracks;
   iEvent.getByLabel(trackingParticleLabel,TPtracks);

   //get the propagator
   iSetup.get<TrackingComponentsRecord>().get(thePropagatorName, thePropagator);

   //get the tracker geometry
   iSetup.get<TrackerDigiGeometryRecord>().get(theTrackerGeometry);

   TrajectoryStateTransform transformer;
   
   //open the L2 collection
   edm::Handle<reco::TrackCollection> L2tracks;
   iEvent.getByLabel(L2Label, L2tracks);

   //associate
   reco::RecoToSimCollection L2recSimColl = theAssociator->associateRecoToSim(L2tracks, TPtracks, &iEvent);

   //for different names of L3 aglorithms
   for (std::vector<edm::InputTag>::iterator l3it=L3Labels.begin(); l3it!=L3Labels.end(); ++l3it){

     // open the muon links for that label
     Handle<reco::MuonTrackLinksCollection> pMuounlinks;
     iEvent.getByLabel(*l3it, pMuounlinks);

     if (theTellMe){
       edm::LogInfo(category)<<iEvent.id().run()<<" "<<iEvent.id().event()
			     <<" for: "<<(*l3it)<<" L2: "<<L2tracks->size()
			     <<" L3: "<<pMuounlinks->size();
     }

     //loop L2 collection. 
     for (uint iL2=0;iL2!=L2tracks->size();iL2++){
       //   make reference to l2 track
       reco::TrackRef refL2(L2tracks, iL2);

       //get the L2 IP state
       FreeTrajectoryState l2State = transformer.initialFreeState(*refL2, field.product());       
       if (l2State.position().mag()==0)
	 {edm::LogError(category)<<"invalid state from l2 track initial state. skipping.\n"<<l2State;continue;}
       
       bool hasBeenReconstructed=false;
       uint matchingLink=0;
       //   search ref in tracklinks
       for (uint iLink=0;iLink!=pMuounlinks->size();iLink++){
	 if ((*pMuounlinks)[iLink].standAloneTrack() == refL2){
	   //this L3muon is coming from this L2 track
	   hasBeenReconstructed =true;
	   matchingLink=iLink;
	 }//same track ref
       }//loop over track links
       
       if (hasBeenReconstructed){
	 edm::LogVerbatim(category)<<"the L2 track: "<<iL2<<" in: "<<L2Label.encode()<<" has been reconstructed in: "<< l3it->encode()<<" at index: "<<matchingLink;
	 reco::TrackRef refL3 = (*pMuounlinks)[matchingLink].globalTrack();
	 reco::TrackRef refTkL3 = (*pMuounlinks)[matchingLink].trackerTrack();

	 FreeTrajectoryState l3State = transformer.initialFreeState(*refL3, field.product());
	 if (l3State.position().mag()==0)
	   {edm::LogError(category)<<"invalid state from l3 track initial state. skipping.\n"<<l3State;continue;}
	   
	 FreeTrajectoryState l3TkState = transformer.initialFreeState(*refTkL3, field.product());
	 if (l3TkState.position().mag()==0 || l3TkState.momentum().mag()==0)
	   {edm::LogError(category)<<"invalid state from l3 tracker track initial state. skipping.\n"<<l3TkState;continue;}
	   
	 edm::LogVerbatim(category)<<"L2 state:\n"
				   <<l2State
				   <<"\n l3 state: \n"
				   <<l3State
				   <<"\n l3 tracker state: \n"
				   <<l3TkState;


	 nodePlotter & node0 = plotter.dir(*l3it).dir("good");
	 node0.element("L3L2pt")->Fill(l3State.momentum().perp(),l2State.momentum().perp());
	 node0.element("TkL3L2pt")->Fill(l3TkState.momentum().perp(),l2State.momentum().perp());
	 node0.element("TkL3L3pt")->Fill(l3TkState.momentum().perp(),l3State.momentum().perp());

	 nodePlotter  & node1 = plotter.dir(*l3it).dir("good").dir("L2Track");
	 one(node1,refL2,l2State);

	 oneLocal(node1, L2recSimColl, refL2, l2State);

	 nodePlotter  & node2 = plotter.dir(*l3it).dir("good").dir("L3Track");
	 one(node2,refL3,l3State);

	 nodePlotter  & node3 = plotter.dir(*l3it).dir("good").dir("L3TkTrack");
	 one(node3,refTkL3,l3TkState);

       }
       else{
	 //it has not been recoed
	 edm::LogVerbatim(category)<<"the L2 track: "<<iL2<<" in: "<< L2Label.encode()<<" has NOT  been reconstructed\n"
				   <<"L2 state:\n"
				   <<l2State;

	 nodePlotter  & node1 = plotter.dir(*l3it).dir("failed").dir("L2Track");
	 one(node1,refL2,l2State);
	 
	 oneLocal(node1, L2recSimColl, refL2, l2State);

       }//l2 not recoed


     }//loop over l2s
   }//loop over L3 tags
}


// ------------ method called once each job just before starting event loop  ------------
void 
WhatWentWrongInL3Muon::beginJob(const edm::EventSetup&)
{
  //book plots
  plotter.setName("L3Muon");

  //for each label
  for (std::vector<edm::InputTag>::iterator l3it=L3Labels.begin(); l3it!=L3Labels.end(); ++l3it){
    nodePlotter & node1 = plotter.registerNode(*l3it);
      
    const uint Nstatus=2;
    std::string state[Nstatus] = { "good", "failed" };
    for (uint iS=0;iS!=Nstatus;iS++){
      nodePlotter & node2 = node1.registerNode(state[iS]);

      const uint Nwhich=3;
      std::string which[Nwhich] = { "L2Track" , "L3Track", "L3TkTrack" };
      if (state[iS] == "good"){
	theDQM->setCurrentFolder(node2.fullName());
	node2.registerElement(theDQM->book2D("L3L2pt","p_{T} in GeV, L2 versus L3",500,0,500, 500,0,500));
	node2.registerElement(theDQM->book2D("TkL3L2pt","p_{T} in GeV, L2 versus L3 tracker leg",500,0,500, 500,0,500));
	node2.registerElement(theDQM->book2D("TkL3L3pt","p_{T} in GeV, L3 versus L3 tracker leg",500,0,500, 500,0,500));
      }
      for (uint iW=0;iW!=Nwhich;iW++){
	if (iW>=1 && state[iS] == "failed") continue;

	nodePlotter & node3= node2.registerNode(which[iW]);

	theDQM->setCurrentFolder(node3.fullName());
	edm::LogVerbatim(category)<<node3.fullName();
	node3.registerElement(theDQM->book1D("nRecHit","number of muon rechits", 100, 0, 100));
	node3.registerElement(theDQM->bookProfile2D("nRecHit_2D","number of muon rechits versus eta and pt", 500,0,500, 50,-5,5, 100, 0, 100));
	node3.registerElement(theDQM->book1D("px","p_{x} in GeV",300,-300,300));
	node3.registerElement(theDQM->book1D("py","p_{y} in GeV",300,-300,300));
	node3.registerElement(theDQM->book1D("pz","p_{z} in GeV",300,-300,300));
	node3.registerElement(theDQM->book1D("pt","p_{T} in GeV",500,0,500));
	node3.registerElement(theDQM->book1D("p","p in GeV",300,-300,300));
	node3.registerElement(theDQM->book1D("m_eta","momentum #eta ",50, -5, 5));
	node3.registerElement(theDQM->book1D("m_phi","momentum #varphi ",50, -TMath::Pi(), TMath::Pi()));
	node3.registerElement(theDQM->book1D("x","x in cm",500,-100,100));
	node3.registerElement(theDQM->book1D("y","y in cm",500,-100,100));
	node3.registerElement(theDQM->book1D("z","z in cm",300,-300,300));
	node3.registerElement(theDQM->book1D("v_eta","position #eta ",50, -5, 5));
	node3.registerElement(theDQM->book1D("v_phi","position #varphi ",50, -TMath::Pi(), TMath::Pi()));
	node3.registerElement(theDQM->book1D("chi2","track #chi^{2}",400,0,400));
	node3.registerElement(theDQM->book1D("Prob","track #chi^{2} probability",100,0,1));
	
	//	if (state[iS] == "failed"){
	  node3.registerElement(theDQM->book1D("lDx","local #Delta_{x}(track-sim) in cm",200,-20,20));
	  node3.registerElement(theDQM->book1D("lDy","local #Delta_{y}(track-sim) in cm",200,-20,20));
	  node3.registerElement(theDQM->book2D("lDxDy","local #Delta_{y} versus #Delta_{x} (track-sim) in cm",200,-20,20, 200,-20,20));
	  node3.registerElement(theDQM->book1D("lD","local distance (track-sim) in cm",200,-20,20));
	  node3.registerElement(theDQM->book1D("lDx_pull","local pull of #Delta_{x}(track-sim)",100,-20,20));
	  node3.registerElement(theDQM->book1D("lDy_pull","local pull of #Delta_{y}(track-sim)",100,-20,20));
	  node3.registerElement(theDQM->book1D("matchId","pdgId of the match track",2000,-1000,1000));

	  node3.registerElement(theDQM->book2D("lDx_rho","local #Delta_{x}(track-sim) in cm function of #rho",200,-20,20, 300,0,300));
	  node3.registerElement(theDQM->book2D("lDy_rho","local #Delta_{y}(track-sim) in cm function of #rho",200,-20,20, 300,0,300));
	  node3.registerElement(theDQM->book2D("lDx_pull_rho","local pull of #Delta_{x}(track-sim) versus #rho",100,-20,20, 300,0,300));
	  node3.registerElement(theDQM->book2D("lDy_pull_rho","local pull of #Delta_{y}(track-sim) versus #rho",100,-20,20, 300,0,300));

	  node3.registerElement(theDQM->book2D("lDx_r","local #Delta_{x}(track-sim) in cm function of r",200,-20,20, 200,0,200));
	  node3.registerElement(theDQM->book2D("lDy_r","local #Delta_{y}(track-sim) in cm function of r",200,-20,20, 200,0,200));
	  node3.registerElement(theDQM->book2D("lDx_pull_r","local pull of #Delta_{x}(track-sim) versus r",100,-20,20, 200,0,200));
	  node3.registerElement(theDQM->book2D("lDy_pull_r","local pull of #Delta_{y}(track-sim) versus r",100,-20,20, 200,0,200));

	  //	}

      }//node3. l2/l3/l3tk
    }//node2. failed/good
  }//node1. label
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
WhatWentWrongInL3Muon::endJob() {
  if (!theRootFileName.empty()){
    theDQM->save(theRootFileName);
  }
}
