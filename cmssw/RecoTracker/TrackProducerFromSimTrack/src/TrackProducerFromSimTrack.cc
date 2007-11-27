// -*- C++ -*-
//
// Package:    TrackProducerFromSimTrack
// Class:      TrackProducerFromSimTrack
// 
/**\class TrackProducerFromSimTrack TrackProducerFromSimTrack.cc RecoTracker/TrackProducerFromSimTrack/src/TrackProducerFromSimTrack.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jean-Roch Vlimant
//         Created:  Fri Apr 13 17:42:20 CDT 2007
// $Id$
//
//

#include "RecoTracker/TrackProducerFromSimTrack/interface/TrackProducerFromSimTrack.h"

//
// constructors and destructor
//
TrackProducerFromSimTrack::TrackProducerFromSimTrack(const edm::ParameterSet& iConfig)
{
  category="TrackProducerFromSimTrack";

  instanceName = iConfig.getParameter<std::string>("instanceName");
   //register your products
  produces<reco::TrackCollection>(instanceName);
  produces<TrackingRecHitCollection>(instanceName); //will be empty for the time being
  produces<reco::TrackExtraCollection>(instanceName);

  //now do what ever other initialization is needed  
  //  errorMatrixRootFile = iConfig.getParameter<std::string>("errorMatrixRootFile");
  matrixProvider_pset = iConfig.getParameter<edm::ParameterSet>("errorMatrix_pset");
  errorMatrixOverEstimate = iConfig.getParameter<double>("errorMatrixOverEstimate");

  plotFileName = iConfig.getParameter<std::string>("plotFileName");
  

  
}


TrackProducerFromSimTrack::~TrackProducerFromSimTrack()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
TrackProducerFromSimTrack::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  ievent=&iEvent;
  isetup=&iSetup;

  //Retreive the collection of simtracks
  iEvent.getByType<edm::SimTrackContainer>(simTracks);
  //retreive the collection of simvertex
  iEvent.getByType<edm::SimVertexContainer>(simVtx);

  //get the magnetic field
  iSetup.get<IdealMagneticFieldRecord>().get(field);

  //get the propagator along momentum
  std::string propagatorAlongName = "SteppingHelixPropagatorAlong";
  iSetup.get<TrackingComponentsRecord>().get(propagatorAlongName,propagatorAlong);

  //prepare the output collection
  std::auto_ptr<reco::TrackCollection> Toutput(new reco::TrackCollection());
  std::auto_ptr<TrackingRecHitCollection> TRHoutput(new TrackingRecHitCollection());
  std::auto_ptr<reco::TrackExtraCollection> TEoutput(new reco::TrackExtraCollection());
  refprodTE = iEvent.getRefBeforePut<reco::TrackExtraCollection>(instanceName);
  TEi=0;
  refprodRH =iEvent.getRefBeforePut<TrackingRecHitCollection>(instanceName);
  RHi=0;

  //loop over the sim track
  edm::LogInfo("TrackProducerFromSimTrack")<<"looking at ("<<simTracks->size()<<") simtracks";
  for (uint isim=0;isim!=simTracks->size();++isim){
    //get a reference to the simtrack
    const SimTrack & simtrack=(*simTracks)[isim];

    //make a selection on it
    if (!selectSimTrack(simtrack)) continue;

    //create a reco track from the simtrack
    reco::Track recotrack = makeTrack(simtrack);

    //make a selection on the create reco::Track
    if (!selectTrack(recotrack)) continue;
    
    //make it in the collection of track
    Toutput->push_back(recotrack);
    reco::Track & recotrackref = Toutput->back();

    //create the track extra
    reco::TrackExtra * extra = makeTrackExtra(simtrack,recotrackref,*TEoutput);
    if (!extra){
      edm::LogError(category)<<"cannot create the track extra for this track.";
      //pop the inserted track
      Toutput->pop_back();
      continue;}

    //attach the collection of rechits
    if (!attachRecHits(simtrack,recotrackref,*extra,*TRHoutput)){
      edm::LogError(category)<<"cannot attach any rechits on this track";
      //pop the inserted track
      Toutput->pop_back();
      //pop the track extra
      TEoutput->pop_back();
      continue;}
    
  }//loop over simtracks
  
  //write things to the event.
  edm::LogInfo(category)<<"writing: \n a collection of ("<< Toutput->size() <<") reco::Track to the event."
			<<"\n a collection of ("<< TEoutput->size() <<") reco::TrackExtra to the event."
			<<"\n a collection of ("<< TRHoutput->size() <<") associated rechits to the event.";
			    
  iEvent.put(Toutput,instanceName);
  iEvent.put(TEoutput,instanceName);
  iEvent.put(TRHoutput,instanceName);
}

// ------------ method called once each job just before starting event loop  ------------



TH1 * TrackProducerFromSimTrack::get(const char * hname){
  TH1* h = (TH1*) plotFile->Get(hname);
  if (!h){edm::LogError(category)<<"could not get: "<<hname<<" from the histogram file. expect troubles. like a seg fault.";}
  return h;}


void 
TrackProducerFromSimTrack::beginJob(const edm::EventSetup&)
{
  //  matrixProvider = new ErrorMatrix(errorMatrixRootFile.c_str());
  matrixProvider = new ErrorMatrix(matrixProvider_pset);

  if (plotFileName =="")
    plotFile=0;
  else{
    
    plotFile = TFile::Open(plotFileName.c_str(),"recreate");
    plotFile->cd();
    //book histograms in this section
    //you will get them back from their name

    const   double pull_max = 3;
    TH1 *hist;

    const TString local_vars[5] = {"0","1","2","3","4"};
    const TString local_units[5] = {"[?]","[?]","[?]","[?]","[?]"};
    const TString local_tags[5] = {"0","1","2","3","4"};
    const int local_n[5] = { 100, 100, 100, 100, 100 };
    const double local_min[5] = { -0.3, -0.3, -0.3, -0.3, -0.3 } ;
    const double local_max[5] = { 0.3, 0.3, 0.3, 0.3, 0.3 } ;

    for (uint ihist=0;ihist!=5;++ihist) {
      hist= new TH1D("hist_local_sim_"+local_tags[ihist],local_vars[ihist]+" distribution for SimTrack",local_n[ihist],local_min[ihist],local_max[ihist]);
      hist->SetXTitle(local_vars[ihist]+" "+local_units[ihist]);
      hist= new TH1D("hist_local_track_"+local_tags[ihist],local_vars[ihist]+" distribution for SimTrack",local_n[ihist],local_min[ihist],local_max[ihist]);
      hist->SetXTitle(local_vars[ihist]+" "+local_units[ihist]);
      //difference
      hist= new TH1D("hist_local_residual_"+local_tags[ihist],local_vars[ihist]+" residual",local_n[ihist],local_min[ihist],local_max[ihist]);
      hist->SetXTitle("#Delta_{gen-reco}("+local_vars[ihist]+") "+local_units[ihist]);
      //difference/error
      hist= new TH1D("hist_local_pull_"+local_tags[ihist],local_vars[ihist]+" pull",local_n[ihist],-pull_max,pull_max);
      hist->SetXTitle("#frac{#Delta_{gen-reco}}{#sigma}("+local_vars[ihist]+") [100%]");
      //error
      hist= new TH1D("hist_local_sigma_"+local_tags[ihist],local_vars[ihist]+" sigma",local_n[ihist],local_min[ihist],local_max[ihist]);
      hist->SetXTitle("#sigma("+local_vars[ihist]+") "+local_units[ihist]);


    }

    const TString rotated_vars[5] = {"0","1","2","3","4"};
    const TString rotated_units[5] = {"[?]","[?]","[?]","[?]","[?]"};
    const TString rotated_tags[5] = {"0","1","2","3","4"};
    const int rotated_n[5] = { 100, 100, 100, 100, 100 };
    const double rotated_min[5] = { -0.3, -0.3, -0.3, -0.3, -0.3 } ;
    const double rotated_max[5] = { 0.3, 0.3, 0.3, 0.3, 0.3 } ;
    for (uint ihist=0;ihist!=5;++ihist) {
      hist= new TH1D("hist_rotated_sim_"+rotated_tags[ihist],rotated_vars[ihist]+" distribution for SimTrack",rotated_n[ihist],rotated_min[ihist],rotated_max[ihist]);
      hist->SetXTitle(rotated_vars[ihist]+" "+rotated_units[ihist]);
      hist= new TH1D("hist_rotated_track_"+rotated_tags[ihist],rotated_vars[ihist]+" distribution for SimTrack",rotated_n[ihist],rotated_min[ihist],rotated_max[ihist]);
      hist->SetXTitle(rotated_vars[ihist]+" "+rotated_units[ihist]);
      //difference
      hist= new TH1D("hist_rotated_residual_"+rotated_tags[ihist],rotated_vars[ihist]+" residual",rotated_n[ihist],rotated_min[ihist],rotated_max[ihist]);
      hist->SetXTitle("#Delta_{gen-reco}("+rotated_vars[ihist]+") "+rotated_units[ihist]);
      //difference/error
      hist= new TH1D("hist_rotated_pull_"+rotated_tags[ihist],rotated_vars[ihist]+" pull",rotated_n[ihist],-pull_max,pull_max);
      hist->SetXTitle("#frac{#Delta_{gen-reco}}{#sigma}("+rotated_vars[ihist]+") [100%]");
      //error
      hist= new TH1D("hist_rotated_sigma_"+rotated_tags[ihist],rotated_vars[ihist]+" sigma",rotated_n[ihist],rotated_min[ihist],rotated_max[ihist]);
      hist->SetXTitle("#sigma("+rotated_vars[ihist]+") "+rotated_units[ihist]);

    }

  }

}

// ------------ method called once each job just after ending the event loop  ------------
void 
TrackProducerFromSimTrack::endJob() {
  if (plotFile){
    plotFile->Write();
    plotFile->Close();}
}

bool TrackProducerFromSimTrack::selectSimTrack(const SimTrack & simtrack){
  //make a selection on the simulated track.
  //pT, eta, type, whatever you want. 
  LogDebug("TrackProducerFromSimTrack")<<"dealing with a particle of type: "<<simtrack.type();
  //only muons
  if (abs(simtrack.type())!=13) return false;
  return true;  
}

reco::TrackExtra * TrackProducerFromSimTrack::makeTrackExtra(const SimTrack & simtrack,
							     reco::Track & recotrack,
							     reco::TrackExtraCollection& TEcol){
  //get the initial free state of the track
  TrajectoryStateTransform transformer;
  TrajectoryStateOnSurface current;
  FreeTrajectoryState initialState = transformer.initialFreeState(recotrack,field.product());
  GlobalPoint point;
  GlobalVector vector;

  /* this all thing could have worked out if surface was stored directly and not via detId...
  //rought geometry of the muon system, to get outer and inner measurement poistion/momentum
  const double r_min=385;
  const double r_max=728;
  const double z_min=721;
  const double z_max=1060;
  Surface::RotationType nullRot;
  Surface::PositionType center(0,0,0);
  Surface::PositionType pos_max(0,0,z_max);
  Surface::PositionType pos_min(0,0,z_min);
  Surface::PositionType neg_max(0,0,-z_max);
  Surface::PositionType neg_min(0,0,-z_min);
  Plane::PlanePointer min_pos_plane = Plane::build(pos_min,nullRot);
  Plane::PlanePointer max_pos_plane = Plane::build(pos_max,nullRot);
  Plane::PlanePointer min_neg_plane = Plane::build(neg_min,nullRot);
  Plane::PlanePointer max_neg_plane = Plane::build(neg_max,nullRot);

  Cylinder::CylinderPointer max_cyl = Cylinder::build(center,nullRot,r_max);
  Cylinder::CylinderPointer min_cyl = Cylinder::build(center,nullRot,r_min);
  
  //outer measurement
  TrajectoryStateOnSurface outertsos;
  current = propagatorAlong->propagate(initialState,*max_cyl);
  if (current.isValid()){
    if (fabs(current.globalPosition().z())<z_min){
      outertsos = current;}}
  if (!outertsos.isValid()){
    current = propagatorAlong->propagate(initialState,
					 ((initialState.momentum().z()>0)?*max_pos_plane:*max_neg_plane));
    if (!current.isValid()){
      edm::LogError(category)<<"cannot get an outer measurement for this track:\n"<<initialState;
      return 0;}
    else{
      outertsos = current;}
  }
  point = outertsos.globalPosition();
  vector = outertsos.globalMomentum();
  unsigned int outerId=0;
  math::XYZVector outmom(vector.x(),vector.y(),vector.z());
  math::XYZPoint  outpos(point.x(),point.y(),point.z());
  
  //inner measurement
  TrajectoryStateOnSurface innertsos;
  current = propagatorAlong->propagate(initialState, *min_cyl);
  if (current.isValid()){
    if (fabs(current.globalPosition().z())<z_min){
      innertsos = current;}}
  if (!innertsos.isValid()){
    current = propagatorAlong->propagate(initialState,
					 ((initialState.momentum().z()>0)?*min_pos_plane:*min_neg_plane));
    if (!current.isValid()){
      edm::LogError(category)<<"cannot get an inner measurement for this track:\n"<<initialState;
      return 0;}
    else{
      innertsos = current;}
  }
  point = innertsos.globalPosition();
  vector = innertsos.globalMomentum();
  unsigned int innerId=0;
  math::XYZVector inmom(vector.x(),vector.y(),vector.z());
  math::XYZPoint  inpos(point.x(),point.y(),point.z()); 

  edmLogInfo(category)<<"outer state for outer measurement:\n"<<outertsos
			<<"\ninner state for inner measurement:\n"<<innertsos;
  */

  //----------------
  //  the story goes:
  //  get the inner and outer most surface of the psimhits associated with this muon
  //----------------

  //get the list of PSimHit for this muon.
  unsigned int id_to_match = simtrack.trackId();
  const unsigned int index_fix=0;
  double dmin=1000000;
  double dmax=0;
  const Surface * outerSurface =0;
  unsigned int outerId=0;
  const Surface * innerSurface =0;
  unsigned int innerId=0;

  const unsigned int Ndetector=3;
  const std::string geometries_name[Ndetector] = {"MuonDTHits", "MuonCSCHits", "MuonRPCHits" };
  const std::string pSimHitSource = "g4SimHits";

  edm::ESHandle<DTGeometry> dtGeomESH;
  edm::ESHandle<CSCGeometry> cscGeomESH;
  edm::ESHandle<RPCGeometry> rpcGeomESH;
  //get the geometry of the muon system
  isetup->get<MuonGeometryRecord>().get(dtGeomESH);
  isetup->get<MuonGeometryRecord>().get(cscGeomESH);
  isetup->get<MuonGeometryRecord>().get(rpcGeomESH);
  
  const TrackingGeometry * geometry=0;
  edm::Handle<edm::PSimHitContainer> simHitsAllMuons;

  for (uint i=0;i!=Ndetector;++i){
    //select the geometry
    if (i==0) geometry = dtGeomESH.product();
    else if (i==1)geometry = cscGeomESH.product();
    else if (i==2)geometry = rpcGeomESH.product();
    else{edm::LogError(category)<<"geometry index going too far."; break;}
    //get the psimhits
    ievent->getByLabel(pSimHitSource,geometries_name[i],simHitsAllMuons);
    //loop over them
    for (edm::PSimHitContainer::const_iterator ips=simHitsAllMuons->begin();ips!=simHitsAllMuons->end();++ips){
      //reject hits from other particles
      if (ips->trackId() != id_to_match + index_fix) continue; 
      //get the det id on which the hit is
      DetId id(ips->detUnitId());
      const GeomDetUnit * gdu = geometry->idToDetUnit(id);
      if (!gdu){edm::LogError(category)<<"cannot get GeomDetUnit: "<<ips->detUnitId(); continue;}
      //get its global position
      GlobalPoint position = gdu->surface().toGlobal(ips->entryPoint());
      double mag= position.mag();
      //select inner most point
      if (mag<dmin){
	dmin = mag; innerId = ips->detUnitId();
	innerSurface = &gdu->surface();}
      //select outer most point
      if (mag>dmax){
	dmax = mag; outerId = ips->detUnitId();
	outerSurface = &gdu->surface();}
      }
  }

  //now make the states outer and inner
  if (!outerSurface) {edm::LogError(category)<<"cannot get outer surface. no extra."; return 0;}
  TrajectoryStateOnSurface outertsos = propagatorAlong->propagate(initialState,*outerSurface);
  if (!outertsos.isValid()){edm::LogError(category)<<"cannot get outer state. no extra."; return 0;}
  point = outertsos.globalPosition();
  vector = outertsos.globalMomentum();
  math::XYZVector outmom(vector.x(),vector.y(),vector.z());
  math::XYZPoint  outpos(point.x(),point.y(),point.z());

  if (!innerSurface) {edm::LogError(category)<<"cannot get inner surface. no extra."; return 0;}
  TrajectoryStateOnSurface innertsos = propagatorAlong->propagate(initialState,*innerSurface);
  if (!innertsos.isValid()){edm::LogError(category)<<"cannot get inner state. no extra."; return 0;}
  point = innertsos.globalPosition();
  vector = innertsos.globalMomentum();
  math::XYZVector inmom(vector.x(),vector.y(),vector.z());
  math::XYZPoint  inpos(point.x(),point.y(),point.z()); 

  LogDebug(category)<<"outer state for outer measurement:\n"<<outertsos
		    <<"\ninner state for inner measurement:\n"<<innertsos;  
  

  TEcol.push_back(reco::TrackExtra (outpos, outmom, true, inpos, inmom, true,
				    outertsos.curvilinearError().matrix(), outerId,
				    innertsos.curvilinearError().matrix(), innerId,
				    alongMomentum));
  //add a reference to the trackextra on the track
  recotrack.setExtra(edm::Ref<reco::TrackExtraCollection>(refprodTE, TEi++));

  //return the reference to the last inserted then
  return &(TEcol.back());}

bool TrackProducerFromSimTrack::attachRecHits(const SimTrack & simtrack,
					      reco::Track & recotrack,
					      reco::TrackExtra & trackextra,
					      TrackingRecHitCollection& RHcol){
  //this function is dummy so far. do not care about this
  //maybe in the future, we would want to fack the rechits, or attached the rechits trully associated with this simtrack
  return true;
}


reco::Track TrackProducerFromSimTrack::makeTrack(const SimTrack & simtrack){
  //concrete implementation of the conversion from SimTrack to reco::Track
  double chi2=0;//do not worry about that yet
  double ndof=0;//do not worry about that yet

  if (simtrack.vertIndex()==-1){
    std::cout <<" this is an error: no vertex."<<std::endl;
    return reco::Track();}
  const SimVertex & vertex=(*simVtx)[simtrack.vertIndex()];
 
  GlobalPoint position(vertex.position().x(),
		       vertex.position().y(),
		       vertex.position().z());
  GlobalVector momentum(simtrack.momentum().px(),
			simtrack.momentum().py(),
			simtrack.momentum().pz());
  int charge = (int) simtrack.charge();
 
  GlobalTrajectoryParameters unsmeared_glb_parameters(position,
						      momentum,
						      charge,
						      field.product());
 

  //Make a FTS from the vertex point and momentum direction obtained from simtrack and simvetex
  
  FreeTrajectoryState origin_state(unsmeared_glb_parameters);
  LogDebug(category)<<"original state:\n"<<origin_state;
 
  //go to pca with analytical extrapolator
  GlobalPoint vtx(0,0,0);
  AnalyticalImpactPointExtrapolator extrapolator(field.product());
  TrajectoryStateOnSurface tsos=extrapolator.extrapolate(origin_state,vtx);
 
  //this FTS is the PCA state. implementation later on
  //  FreeTrajectoryState PCAstate=origin_state;
  FreeTrajectoryState PCAstate=*tsos.freeState();
  LogDebug(category)<<"PCA state:\n"<<PCAstate;
 
  //compute the implicit surface of the FTS
  Surface::PositionType center = PCAstate.position();
  //axis from the convention cartesian->local
  GlobalVector aZ = PCAstate.momentum().unit();
  GlobalVector aX = GlobalVector(aZ.y(), -aZ.x() ,0).unit();
  GlobalVector aY = aZ.cross(aX);
 
  Surface::RotationType rot(aX,aY,aZ);
  Plane::PlanePointer plane = Plane::build(center,rot);
 
  //convert global things into local things
  LogDebug("TrackProducerFromSimTrack")<<"global position: "<<PCAstate.position()<<" global momentum: "<<PCAstate.momentum();
  LocalVector lmomentum = plane->toLocal(PCAstate.momentum());
  LocalPoint lposition = plane->toLocal(PCAstate.position()); 
  LogDebug(category)<<"local position: "<<lposition<<" local momentum: "<<lmomentum;
 
  //get the 5 local parameters vector
  LocalTrajectoryParameters local_parameters(lposition,lmomentum,PCAstate.charge());
  AlgebraicVector unsmeared_parameters = local_parameters.vector_old();
  LogDebug(category)<<"unsmeared parameters: "<<unsmeared_parameters;
 
  //those are the parameters of the smearing
  CurvilinearTrajectoryError curv_matrix = smearingParameters(origin_state);
  
  AlgebraicSymMatrix c_smearing_parameters = curv_matrix.matrix_old();

  /*  AlgebraicSymMatrix c_smearing_parameters(5); 
      c_smearing_parameters[0][0] = 0.01;
      c_smearing_parameters[1][1] = 0.01;
      c_smearing_parameters[2][2] = 0.01;
      c_smearing_parameters[3][3] = 0.01;
      c_smearing_parameters[4][4] = 0.01;
      CurvilinearTrajectoryError curv_matrix(c_smearing_parameters);
  */
  
  AlgebraicSymMatrix smearing_parameters = smearingParameters(curv_matrix,origin_state,*plane);
 
  LogDebug(category)<<"smearing parameters: "<<smearing_parameters;
 
  //do the smearing
  AlgebraicVector smeared_parameters = smear(unsmeared_parameters, smearing_parameters);
 
  LogDebug(category)<<"smeared parameters: "<<smeared_parameters;
 
  //get the smeared local parameters
  LocalTrajectoryParameters smeared_local_parameters(smeared_parameters,
						     local_parameters.pzSign());
 

  if (plotFile){
  for (int i=0;i!=5;++i){
    get(Form("hist_local_sim_%d",i))->Fill(unsmeared_parameters[i]);
    get(Form("hist_local_track_%d",i))->Fill(smeared_parameters[i]);
    get(Form("hist_local_sigma_%d",i))->Fill(sqrt(c_smearing_parameters[i][i]));
    get(Form("hist_local_sigma_%d",i))->Fill(sqrt(c_smearing_parameters[i][i]));
    get(Form("hist_local_residual_%d",i))->Fill(unsmeared_parameters[i]-smeared_parameters[i]);
    get(Form("hist_local_pull_%d",i))->Fill((unsmeared_parameters[i]-smeared_parameters[i])/sqrt(c_smearing_parameters[i][i]));
    get(Form("hist_local_pull_%d",i))->Fill((unsmeared_parameters[i]-smeared_parameters[i])/sqrt(c_smearing_parameters[i][i]));
  }}
  
 
  //convert them back to global frame
  LogDebug(category)<<"smeared local position: "<<smeared_local_parameters.position()<<" smeared local momentum: "<<smeared_local_parameters.momentum();
  LogDebug(category)<<"smeared global position: "<<plane->toGlobal(smeared_local_parameters.position())<<" smeared global momentum: "<<plane->toGlobal(smeared_local_parameters.momentum());
 
  //this is the charge of the output
  int refcharge = smear_charge(origin_state);  
 
  //make a free trajectory state back
  GlobalTrajectoryParameters smeared_glb_parameters(plane->toGlobal(smeared_local_parameters.position()),
						   plane->toGlobal(smeared_local_parameters.momentum()),
						   refcharge,
						   field.product());
  FreeTrajectoryState smeared_origin_state(smeared_glb_parameters,curv_matrix);
  LogDebug(category)<<"smeared state at its origin:\n"<<smeared_origin_state;
 

  /*
     //get it to it pca again
     TrajectoryStateOnSurface smeared_pca_state = extrapolator.extrapolate(smeared_origin_state,vtx);
     LogDebug(category)<<"smeared state at its pca:\n"<<smeared_pca_state;
     GlobalPoint refpos = smeared_pca_state.globalPosition();
     GlobalVector refmom = smeared_pca_state.globalMomentum();
     
     //construct the track parameters
     reco::TrackBase::Point reference_point(refpos.x(),
     refpos.y(),
     refpos.z());
     reco::TrackBase::Vector reference_momentum(refmom.x(),
     refmom.y(),
     refmom.z());

     //should be return by a function. 
     reco::TrackBase::CovarianceMatrix covariance_matrix;
     
     //  assign propagated smearing parameters for the time being
     AlgebraicSymMatrix refError = smeared_pca_state.curvilinearError().matrix_old();
     
     for (int i=0;i!=5;++i){for(int j=0;j!=5;++j){covariance_matrix(i,j)=errorMatrixOverEstimate*refError[i][j];}}
     
     
     edm::LogInfo(category)<<"From SimTrack (PCA), x: "<<PCAstate.position()<<" p: "<<PCAstate.momentum()<<"\n"
     <<" making reco::Track, x: "<<smeared_pca_state.globalPosition()<<" p: "<<smeared_pca_state.globalMomentum()<<"\n"
     <<covariance_matrix;
     
     return reco::Track(chi2,ndof,reference_point,reference_momentum,refcharge,covariance_matrix);
*/
  

  TSCPBuilderNoMaterial tscpBuilder;
  TrajectoryStateClosestToPoint tscp = tscpBuilder(smeared_origin_state,
						   GlobalPoint(0,0,0));

  GlobalPoint refpos = tscp.theState().position();
  GlobalVector refmom = tscp.theState().momentum();
  //construct the track parameters
  reco::TrackBase::Point reference_point(refpos.x(),
					 refpos.y(),
					 refpos.z());
  reco::TrackBase::Vector reference_momentum(refmom.x(),
					     refmom.y(),
					     refmom.z());

  //should be return by a function.
  reco::TrackBase::CovarianceMatrix covariance_matrix;
  //  assign propagated smearing parameters for the time being
  AlgebraicSymMatrix refError = tscp.theState().curvilinearError().matrix_old();

  for (int i=0;i!=5;++i){for(int j=0;j!=5;++j){covariance_matrix(i,j)=errorMatrixOverEstimate*refError[i][j];}}

  //get it to it pca again to check it.
  TrajectoryStateOnSurface smeared_pca_state = extrapolator.extrapolate(smeared_origin_state,vtx);

  edm::LogInfo(category)<<"From SimTrack (PCA), x: "<<PCAstate.position()<<" p: "<<PCAstate.momentum()<<"\n"
			<<" (GEN)               x: "<<origin_state.position()<<" p: "<<origin_state.momentum()<<"\n"
			<<"make reco::Track(PCA), x: "<<smeared_pca_state.globalPosition()<<" p: "<<smeared_pca_state.globalMomentum()<<"\n"
			<<" (TSPC)                x: "<<tscp.theState().position()<<" p: "<<tscp.theState().momentum()<<"\n"
			<<covariance_matrix<<"\n";

  return reco::Track(chi2,ndof,reference_point,reference_momentum,refcharge,covariance_matrix);
}




bool TrackProducerFromSimTrack::selectTrack(const reco::Track & recotrack){
  //make a selection on the simulated track.
  //pT, eta, whatever you want. 
  return true;
}





int TrackProducerFromSimTrack::smear_charge(const FreeTrajectoryState & state)
{
  //will implement later on
  return state.charge();
}

CurvilinearTrajectoryError TrackProducerFromSimTrack::smearingParameters(const FreeTrajectoryState & state){
  return matrixProvider->get(state.momentum());}

AlgebraicSymMatrix TrackProducerFromSimTrack::smearingParameters(const CurvilinearTrajectoryError & curv_matrix, const FreeTrajectoryState & state,const Plane & plane){
  LogDebug(category)<<"errror matrix from provider: "<<curv_matrix.matrix_old();
  
  TrajectoryStateOnSurface tsos(state,plane);
  JacobianCurvilinearToLocal j(plane,tsos.localParameters(),*field.product());
  
  return curv_matrix.matrix_old().similarity(j.jacobian_old());}



AlgebraicVector TrackProducerFromSimTrack::smear(const AlgebraicVector & vector, const AlgebraicSymMatrix & matrix){
 
  int ifail=0;
 
  LogDebug(category)<<"original matrix: "<<matrix;
 
  AlgebraicSymMatrix inverse_cov = matrix.inverse(ifail);
  //check on ifail maybe
  LogDebug(category)<<"inverted matrix: "<<inverse_cov<<"\n ifail status: "<<ifail;
 
  //get U and delta
  HepSymMatrix delta = inverse_cov;
  HepMatrix U_matrix = diagonalize(&delta);
  LogDebug(category)<<"diagonalized matrix: "<<delta<<"\nU matrix: "<<U_matrix;
 
  //get x_0'
  LogDebug(category)<<"vector of parameter in the original base: "<<vector;
  AlgebraicVector vector_prime =U_matrix.T()*vector;
  LogDebug(category)<<"vector of parameter in the rotated base: "<<vector_prime;
 
  //smear the parameters one by one according to corresponding sigma in the 
  static TRandom2 rndGenerator;
  AlgebraicVector vector_prime_smeared(5);
  for (int i=1; i<=5; ++i){
    double unsmeared_parameter = vector_prime(i);
    double delta_i = delta(i,i);
    if (delta_i==0){ edm::LogError(category)<<"eigen value of diagonalized weight matrix is zero. skipping the parameter"; continue;}
    if (delta_i<0) { edm::LogError(category)<<"eigen value of diagonalized weight matrix is negative. skipping the parameter"; continue;}
    double sigma_parameter = 1./sqrt(delta_i);
    LogDebug(category)<<"smearing around: "<<unsmeared_parameter<<" with width: "<<sigma_parameter;
    double smeared_parameter = rndGenerator.Gaus(unsmeared_parameter, sigma_parameter);
    LogDebug(category)<<"smeared result: "<<smeared_parameter;
    vector_prime_smeared(i)=smeared_parameter;
 
  }
 

  if (plotFile){
  for (int i=0;i!=5;++i){
    get(Form("hist_rotated_sim_%d",i))->Fill(vector_prime[i]);
    get(Form("hist_rotated_track_%d",i))->Fill(vector_prime_smeared[i]);
    get(Form("hist_rotated_sigma_%d",i))->Fill(1./sqrt(delta[i][i]));
    get(Form("hist_rotated_residual_%d",i))->Fill(vector_prime[i]-vector_prime_smeared[i]);
    get(Form("hist_rotated_pull_%d",i))->Fill((vector_prime[i]-vector_prime_smeared[i])*sqrt(delta[i][i]));
  }}
 



    
  //get back to original parameter reference
  LogDebug(category)<<"vector of smeared parameter in the rotated base: "<<vector_prime_smeared;
  AlgebraicVector vector_smeared = U_matrix*vector_prime_smeared;
  LogDebug(category)<<"vector of smeared parameter in the original base: "<<vector_smeared;      

  return vector_smeared;
}

