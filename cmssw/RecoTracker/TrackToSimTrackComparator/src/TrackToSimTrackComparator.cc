// -*- C++ -*-
//
// Package:    TrackToSimTrackComparator
// Class:      TrackToSimTrackComparator
// 
/**\class TrackToSimTrackComparator TrackToSimTrackComparator.cc RecoTracker/TrackToSimTrackComparator/src/TrackToSimTrackComparator.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jean-Roch Vlimant
//         Created:  Fri May  4 19:07:59 CDT 2007
// $Id: TrackToSimTrackComparator.cc,v 1.2 2007/11/27 20:47:35 vlimant Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include <FWCore/Framework/interface/EventSetup.h>
#include <FWCore/Framework/interface/ESHandle.h>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h>
#include <TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h>
#include <TrackingTools/PatternTools/interface/TSCPBuilderNoMaterial.h>


#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

#include <MagneticField/Records/interface/IdealMagneticFieldRecord.h>
#include <MagneticField/Engine/interface/MagneticField.h>

#include <TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h>

#include <DataFormats/GeometrySurface/interface/Plane.h>

#include "DQMServices/Core/interface/DaqMonitorBEInterface.h"
#include "DQMServices/Daemon/interface/MonitorDaemon.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TDirectory.h"


//
// class decleration
//

class TrackToSimTrackComparator : public edm::EDAnalyzer {
   public:
      explicit TrackToSimTrackComparator(const edm::ParameterSet&);
      ~TrackToSimTrackComparator();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  TH1* get(const char * hname);
  TH1F * book1D(TString name, 
		TString title, 
		int nchX, double lowX, double highX);

  // ----------member data ---------------------------
  edm::InputTag theTrackLabel;
  edm::InputTag theTrackingParticleLabel;

  //name of the associator
  std::string theAssocLabel;
  //add this for track association
  edm::ESHandle<TrackAssociatorBase> theAssociator;

  edm::ESHandle<MagneticField> theField;

  std::string thePlotDirectoryName;
  std::string theCategory;

  bool doWithExtrapolator;
  bool doWithPerigee;
  bool doWithLocal;
  bool doWithPlane;

  std::map<std::string, TH1*> whoIsWhat;
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
TrackToSimTrackComparator::TrackToSimTrackComparator(const edm::ParameterSet& iConfig)
{
  theCategory = "TrackToSimTrackComparator";

   //now do what ever initialization is needed
  thePlotDirectoryName = iConfig.getParameter<std::string>("plotDirectoryName");

  theAssocLabel = iConfig.getParameter<std::string>("associatorName");

  theTrackLabel = iConfig.getParameter<edm::InputTag>("trackLabel");
  theTrackingParticleLabel = iConfig.getParameter<edm::InputTag>("trackingParticleLabel");

  doWithPerigee=true;
  doWithExtrapolator=false;
  doWithLocal=false;
  doWithPlane=false;

}


TrackToSimTrackComparator::~TrackToSimTrackComparator()
{
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
TrackToSimTrackComparator::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   //retreive track collection
   edm::Handle<reco::TrackCollection> tracks;
   iEvent.getByLabel(theTrackLabel,tracks);

   //open a collection of Tracking particle
   edm::Handle<TrackingParticleCollection> TPtracks;
   iEvent.getByLabel(theTrackingParticleLabel,TPtracks);

   //get the associator
   iSetup.get<TrackAssociatorRecord>().get(theAssocLabel,theAssociator);

   LogDebug(theCategory)<<"going to make the association.";
   //associate
   reco::RecoToSimCollection recSimColl = theAssociator->associateRecoToSim(tracks,TPtracks, &iEvent);
   
   LogDebug(theCategory)<<"I have found: "<<recSimColl.size()<<" associations in total.\n"
			<<"from: "<<tracks->size()<<" reco Tracks.";

   //get the mag field
   iSetup.get<IdealMagneticFieldRecord>().get(theField);
   
   //now loop over the association
   for (reco::RecoToSimCollection::const_iterator RtSit=recSimColl.begin();RtSit!=recSimColl.end();++RtSit) {
     const std::vector<std::pair<TrackingParticleRef,double> > & tp = RtSit->val;
       
     if(tp.size()==0) {continue;}

     //get references for convenience
     //     const SimTrack & sim = tp.begin()->first->g4Tracks().front();
     const reco::Track & track = *RtSit->key;
     const TrackingParticleRef & trp = tp.begin()->first;
     
     //-----------------
     //       the following is just supposed to put the two states in the same plane
     //       that might be silly though
     //-----------------
     //first take the reco track state at pca
     
     TrajectoryStateTransform transformer;
     FreeTrajectoryState trackTSCPstate = transformer.initialFreeState(track,theField.product());
     
     if (trackTSCPstate.position().mag()==0)
       {edm::LogError(theCategory)<<"invalid state from track intial state.skipping."; continue;}

     AlgebraicVector6 trackTSCPstate_vector =  trackTSCPstate.parameters().vector();
     AlgebraicSymMatrix66 global_FTS_comp = trackTSCPstate.cartesianError().matrix();
     
     //now the simtrack is at its origin
     GlobalPoint position_sim_fts(trp->vx(),
				  trp->vy(),
				  trp->vz());
     GlobalVector momentum_sim_fts(trp->px(),
				   trp->py(),
				   trp->pz());
     int charge_sim = (int)trp->charge();
     
     GlobalTrajectoryParameters par_sim(position_sim_fts,momentum_sim_fts,charge_sim,theField.product()); 
     FreeTrajectoryState simOriginalstate(par_sim);

     LogDebug(theCategory)<<"track initial state:\n"<<trackTSCPstate
			  <<"\n sim initial state:\n"<<simOriginalstate;
     
     get("h_fts_origin_state_x")->Fill(simOriginalstate.position().x());
     get("h_fts_origin_state_y")->Fill(simOriginalstate.position().y());
     get("h_fts_origin_state_z")->Fill(simOriginalstate.position().z());
     get("h_track_pT")->Fill(simOriginalstate.momentum().perp());

     GlobalPoint vtx(0,0,0);
     
     if (doWithExtrapolator){
       
       //take it at its pca with extrapolator
       AnalyticalImpactPointExtrapolator extrapolator(theField.product());
       FreeTrajectoryState simPCAstate = *extrapolator.extrapolate(simOriginalstate,vtx).freeState();
       
       AlgebraicVector6 simPCAstate_vector =  simPCAstate.parameters().vector();
       
       //compare FreeTrajectory states in global frame
       get("hist_global_FTS_sim_x")->Fill(simPCAstate_vector[0]);
       get("hist_global_FTS_track_x")->Fill(trackTSCPstate_vector[0]);
       get("hist_global_FTS_residual_x")->Fill(simPCAstate_vector[0]-trackTSCPstate_vector[0]);	
       if(global_FTS_comp(0,0)>0){
	 get("hist_global_FTS_sigma_x")->Fill(sqrt(global_FTS_comp(0,0)));
	 get("hist_global_FTS_pull_x")->Fill((simPCAstate_vector[0]-trackTSCPstate_vector[0])/sqrt(global_FTS_comp(0,0)));
       }
       
       
       get("hist_global_FTS_sim_y")->Fill(simPCAstate_vector[1]);
       get("hist_global_FTS_track_y")->Fill(trackTSCPstate_vector[1]);
       get("hist_global_FTS_residual_y")->Fill(simPCAstate_vector[1]-trackTSCPstate_vector[1]);
       if(global_FTS_comp(1,1)>0){
	 get("hist_global_FTS_sigma_y")->Fill(sqrt(global_FTS_comp(1,1)));
	 get("hist_global_FTS_pull_y")->Fill((simPCAstate_vector[1]-trackTSCPstate_vector[1])/sqrt(global_FTS_comp(1,1)));
       } 
       
       get("hist_global_FTS_sim_z")->Fill(simPCAstate_vector[2]);
       get("hist_global_FTS_track_z")->Fill(trackTSCPstate_vector[2]);
       get("hist_global_FTS_residual_z")->Fill(simPCAstate_vector[2]-trackTSCPstate_vector[2]);
       if(global_FTS_comp(2,2)>0){
	 get("hist_global_FTS_sigma_z")->Fill(sqrt(global_FTS_comp(2,2)));
	 get("hist_global_FTS_pull_z")->Fill((simPCAstate_vector[2]-trackTSCPstate_vector[2])/sqrt(global_FTS_comp(2,2)));
       }
	 
       get("hist_global_FTS_sim_px")->Fill(simPCAstate_vector[3]);
       get("hist_global_FTS_track_px")->Fill(trackTSCPstate_vector[3]);
       get("hist_global_FTS_residual_px")->Fill(simPCAstate_vector[3]-trackTSCPstate_vector[3]);
       if(global_FTS_comp(3,3)>0){
	 get("hist_global_FTS_sigma_px")->Fill(sqrt(global_FTS_comp(3,3)));
	 get("hist_global_FTS_pull_px")->Fill((simPCAstate_vector[3]-trackTSCPstate_vector[3])/sqrt(global_FTS_comp(3,3)));
       }
       
       get("hist_global_FTS_sim_py")->Fill(simPCAstate_vector[4]);
       get("hist_global_FTS_track_py")->Fill(trackTSCPstate_vector[4]);
       get("hist_global_FTS_residual_py")->Fill(simPCAstate_vector[4]-trackTSCPstate_vector[4]);
       if(global_FTS_comp(4,4)>0){
	 get("hist_global_FTS_sigma_py")->Fill(sqrt(global_FTS_comp(4,4)));
	 get("hist_global_FTS_pull_py")->Fill((simPCAstate_vector[4]-trackTSCPstate_vector[4])/sqrt(global_FTS_comp(4,4)));
       }
       
       get("hist_global_FTS_sim_pz")->Fill(simPCAstate_vector[5]);
       get("hist_global_FTS_track_pz")->Fill(trackTSCPstate_vector[5]);
       get("hist_global_FTS_residual_pz")->Fill(simPCAstate_vector[5]-trackTSCPstate_vector[5]);
       if(global_FTS_comp(5,5)>0){
	 get("hist_global_FTS_sigma_pz")->Fill(sqrt(global_FTS_comp(5,5)));
	 get("hist_global_FTS_pull_pz")->Fill((simPCAstate_vector[5]-trackTSCPstate_vector[5])/sqrt(global_FTS_comp(5,5)));
       }
       


     if (doWithPlane){
       //generate the plane according to unique convention
       //compute the implicit surface of the FTS
       Surface::PositionType center = simPCAstate.position();
       //axis from the convention cartesian->local
       GlobalVector aZ = simPCAstate.momentum().unit();
       GlobalVector aX = GlobalVector(aZ.y(), -aZ.x() ,0).unit();
       GlobalVector aY = aZ.cross(aX);

       Surface::RotationType rot(aX,aY,aZ);
       Plane::PlanePointer plane = Plane::build(center,rot);
       
       //now put both trajectory in the same frame
       TrajectoryStateOnSurface simPCAtsos(simPCAstate,*plane);
       //       TrajectoryStateOnSurface trackPCAtsos(trackTSCPstate,*plane); //this is completely silly local z will be non-zero
     
       AlgebraicVector5 simPCAtsos_vector = simPCAtsos.localParameters().vector();

	 
       AnalyticalPropagator propagator(theField.product(),anyDirection);
       TrajectoryStateOnSurface trackPCAtsos = propagator.propagate(trackTSCPstate,*plane);

       AlgebraicVector5 trackPCAtsos_vector = trackPCAtsos.localParameters().vector();
       AlgebraicSymMatrix55 local_covariance = trackPCAtsos.localError().matrix();

       //compare TrajectoryStateOnSurface states in local frame

       get("hist_local_TSOS_sim_q_p")->Fill(simPCAtsos_vector[0]);
       get("hist_local_TSOS_track_q_p")->Fill(trackPCAtsos_vector[0]);
       get("hist_local_TSOS_residual_q_p")->Fill(simPCAtsos_vector[0]-trackPCAtsos_vector[0]);	
       if(local_covariance(0,0)>0){
	 get("hist_local_TSOS_sigma_q_p")->Fill(sqrt(local_covariance(0,0)));
	 get("hist_local_TSOS_pull_q_p")->Fill((simPCAtsos_vector[0]-trackPCAtsos_vector[0])/sqrt(local_covariance(0,0)));
       }
	 

       get("hist_local_TSOS_sim_dx_dz")->Fill(simPCAtsos_vector[1]);
       get("hist_local_TSOS_track_dx_dz")->Fill(trackPCAtsos_vector[1]);
       get("hist_local_TSOS_residual_dx_dz")->Fill(simPCAtsos_vector[1]-trackPCAtsos_vector[1]);
       if(local_covariance(1,1)>0){
	 get("hist_local_TSOS_sigma_dx_dz")->Fill(sqrt(local_covariance(1,1)));
	 get("hist_local_TSOS_pull_dx_dz")->Fill((simPCAtsos_vector[1]-trackPCAtsos_vector[1])/sqrt(local_covariance(1,1)));
       } 
	 
       get("hist_local_TSOS_sim_dy_dz")->Fill(simPCAtsos_vector[2]);
       get("hist_local_TSOS_track_dy_dz")->Fill(trackPCAtsos_vector[2]);
       get("hist_local_TSOS_residual_dy_dz")->Fill(simPCAtsos_vector[2]-trackPCAtsos_vector[2]);
       if(local_covariance(2,2)>0){
	 get("hist_local_TSOS_sigma_dy_dz")->Fill(sqrt(local_covariance(2,2)));
	 get("hist_local_TSOS_pull_dy_dz")->Fill((simPCAtsos_vector[2]-trackPCAtsos_vector[2])/sqrt(local_covariance(2,2)));
       }
	 
       get("hist_local_TSOS_sim_x")->Fill(simPCAtsos_vector[3]);
       get("hist_local_TSOS_track_x")->Fill(trackPCAtsos_vector[3]);
       get("hist_local_TSOS_residual_x")->Fill(simPCAtsos_vector[3]-trackPCAtsos_vector[3]);
       if(local_covariance(3,3)>0){
	 get("hist_local_TSOS_sigma_x")->Fill(sqrt(local_covariance(3,3)));
	 get("hist_local_TSOS_pull_x")->Fill((simPCAtsos_vector[3]-trackPCAtsos_vector[3])/sqrt(local_covariance(3,3)));
       }

       get("hist_local_TSOS_sim_y")->Fill(simPCAtsos_vector[4]);
       get("hist_local_TSOS_track_y")->Fill(trackPCAtsos_vector[4]);
       get("hist_local_TSOS_residual_y")->Fill(simPCAtsos_vector[4]-trackPCAtsos_vector[4]);
       if(local_covariance(4,4)>0){
	 get("hist_local_TSOS_sigma_y")->Fill(sqrt(local_covariance(4,4)));
	 get("hist_local_TSOS_pull_y")->Fill((simPCAtsos_vector[4]-trackPCAtsos_vector[4])/sqrt(local_covariance(4,4)));
       }

       //you can compare them now, they are are consistently defined on the same plane
       //does not mean same origin though
       
       LogDebug(theCategory)<<"SimTrack state is:\n"<<simPCAtsos;
       LogDebug(theCategory)<<"reco::Track state is:\n"<<trackPCAtsos;
     
     }
     }

     if(doWithPerigee){
       //-----------------
       //       the following is to get perigee parameters
       //-----------------

       TSCPBuilderNoMaterial tscpBuilder;
       //       TrajectoryStateClosestToPoint simTSCP = tscpBuilder(simPCAstate,vtx);
       TrajectoryStateClosestToPoint simTSCP = tscpBuilder(simOriginalstate,vtx);
       TrajectoryStateClosestToPoint trackTSCP = tscpBuilder(trackTSCPstate,vtx);

       //now you can compare the perigee parameters
       //rho :transverse curvature (signed) 
       //theta, phi
       //transverse impact point paramater
       //longitudinal impact point parameter
       //the error matrixes might be defined in the same frame. so one can make pulls

       get("h_perigee_sim_rho")->Fill(simTSCP.perigeeParameters().transverseCurvature());
       get("h_perigee_track_rho")->Fill(trackTSCP.perigeeParameters().transverseCurvature());
       get("h_perigee_residual_rho")->Fill(simTSCP.perigeeParameters().transverseCurvature()-trackTSCP.perigeeParameters().transverseCurvature());
       double rho_e=trackTSCP.perigeeError().transverseCurvatureError();
       get("h_perigee_sigma_rho")->Fill(rho_e);
       if (rho_e!=0)
	 get("h_perigee_pull_rho")->Fill((simTSCP.perigeeParameters().transverseCurvature()-trackTSCP.perigeeParameters().transverseCurvature())/rho_e);

       get("h_perigee_sim_theta")->Fill(simTSCP.perigeeParameters().theta());
       get("h_perigee_track_theta")->Fill(trackTSCP.perigeeParameters().theta());
       get("h_perigee_residual_theta")->Fill(simTSCP.perigeeParameters().theta()-trackTSCP.perigeeParameters().theta());
       double theta_e=trackTSCP.perigeeError().thetaError();
       get("h_perigee_sigma_theta")->Fill(theta_e);
       if (theta_e!=0)
	 get("h_perigee_pull_theta")->Fill((simTSCP.perigeeParameters().theta()-trackTSCP.perigeeParameters().theta())/theta_e);



       get("h_perigee_sim_phi")->Fill(simTSCP.perigeeParameters().phi());
       get("h_perigee_track_phi")->Fill(trackTSCP.perigeeParameters().phi());
       get("h_perigee_residual_phi")->Fill(simTSCP.perigeeParameters().phi()-trackTSCP.perigeeParameters().phi());
       double phi_e=trackTSCP.perigeeError().phiError();
       get("h_perigee_sigma_phi")->Fill(phi_e);
       if (phi_e!=0)
	 get("h_perigee_pull_phi")->Fill((simTSCP.perigeeParameters().phi()-trackTSCP.perigeeParameters().phi())/phi_e);

       
       get("h_perigee_sim_d0")->Fill(simTSCP.perigeeParameters().transverseImpactParameter());
       get("h_perigee_track_d0")->Fill(trackTSCP.perigeeParameters().transverseImpactParameter());
       get("h_perigee_residual_d0")->Fill(simTSCP.perigeeParameters().transverseImpactParameter()-trackTSCP.perigeeParameters().transverseImpactParameter());
       double d0_e=trackTSCP.perigeeError().transverseImpactParameterError();
       get("h_perigee_sigma_d0")->Fill(d0_e);
       if (d0_e!=0)
	 get("h_perigee_pull_d0")->Fill((simTSCP.perigeeParameters().transverseImpactParameter()-trackTSCP.perigeeParameters().transverseImpactParameter())/d0_e);



       get("h_perigee_sim_z0")->Fill(simTSCP.perigeeParameters().longitudinalImpactParameter());
       get("h_perigee_track_z0")->Fill(trackTSCP.perigeeParameters().longitudinalImpactParameter());
       get("h_perigee_residual_z0")->Fill(simTSCP.perigeeParameters().longitudinalImpactParameter()-trackTSCP.perigeeParameters().longitudinalImpactParameter());
       double z0_e=trackTSCP.perigeeError().longitudinalImpactParameterError();
       get("h_perigee_sigma_z0")->Fill(z0_e);
       if (z0_e!=0)
	 get("h_perigee_pull_z0")->Fill((simTSCP.perigeeParameters().longitudinalImpactParameter()-trackTSCP.perigeeParameters().longitudinalImpactParameter())/z0_e);
     }

   }//both reference are not empty
}

TH1* TrackToSimTrackComparator::get(const char * hname){
  std::map<std::string, TH1*>::iterator f = whoIsWhat.find(hname);
  if (f!=whoIsWhat.end()){
    return f->second;}
  else {
    edm::LogError(theCategory)<<"\n\n\n\n\n################################################################################\n"
			      <<"################################################################################\n"
			      <<"################################################################################\n"
			      <<"could not get: "<<hname<<". expect troubles. like a seg fault."
			      <<"################################################################################\n";
    return 0;
  }
}


TH1F * TrackToSimTrackComparator::book1D(TString n, 
					 TString t, 
					 int nchX, double lowX, double highX)
{
  std::string name(n);
  std::string title(t);
  MonitorElement * me = edm::Service<DaqMonitorBEInterface>()->book1D(name, title, nchX, lowX, highX);
  TH1F * h=dynamic_cast<TH1F*>(&(**((MonitorElementRootH2 *)me)));
  whoIsWhat[name]=h;
  return h;
}



// ------------ method called once each job just before starting event loop  ------------
void 
TrackToSimTrackComparator::beginJob(const edm::EventSetup& setup)
{
  edm::Service<DaqMonitorBEInterface>()->setCurrentFolder(thePlotDirectoryName);

  //book histograms in this section
  //you will get them back from their name
  TH1 * h;
 
  h= book1D("h_track_pT","reco::Track transverse momentum",200,0,200);
  h->SetXTitle("track p_{T} [GeV]");
  // and so on

  h= book1D("h_fts_origin_state_x","x distribution for fts origin_state of SimTrack",100,-.006,.006);
  h->SetXTitle("x [cm]");
  h= book1D("h_fts_origin_state_y","y distribution for fts origin_state of SimTrack",100,-.006,.006);
  h->SetXTitle("y [cm]");
  h= book1D("h_fts_origin_state_z","z distribution for fts origin_state of SimTrack",100,-.006,.006);
  h->SetXTitle("z [cm]");


  if (doWithPerigee){

    const TString perigee_vars[5] = {"#rho","#theta","#varphi","d_{0}","z_{0}"};
    const TString perigee_units[5] = {"","[rad]","[rad]","[cm]","[cm]"};
    const TString perigee_tags[5] = {"rho","theta","phi","d0","z0"};
    const int perigee_n[5] = { 100, 100, 100, 100, 100 };
    const double perigee_min[5] = { -0.003, -TMath::Pi()/2., -TMath::Pi()/2., -0.1, -30. } ;
    const double perigee_max[5] = { 0.003, TMath::Pi()/2., TMath::Pi()/2., 0.1, 30. } ;
    const double perigee_residual_min[5] = {-0.003, -0.1, -0.5, -.2, -5};
    const double perigee_residual_max[5] ={0.003, 0.1, 0.5, .2, 5};

    for (uint ih=0;ih!=5;++ih) {
      h= book1D("h_perigee_sim_"+perigee_tags[ih],perigee_vars[ih]+" distribution for SimTrack",perigee_n[ih],perigee_min[ih],perigee_max[ih]);
      h->SetXTitle(perigee_vars[ih]+" "+perigee_units[ih]);
      h= book1D("h_perigee_track_"+perigee_tags[ih],perigee_vars[ih]+" distribution for Track",perigee_n[ih],perigee_min[ih],perigee_max[ih]);
      h->SetXTitle(perigee_vars[ih]+" "+perigee_units[ih]);
      //difference



      //difference/error
      h= book1D("h_perigee_pull_"+perigee_tags[ih],perigee_vars[ih]+" pull",perigee_n[ih],-6,6);
      h->SetXTitle("#frac{#Delta_{gen-reco}}{#sigma}("+perigee_vars[ih]+") [100%]");
      //error
      h= book1D("h_perigee_sigma_"+perigee_tags[ih],perigee_vars[ih]+" sigma",perigee_n[ih],perigee_min[ih],perigee_max[ih]);
      h->SetXTitle("#sigma("+perigee_vars[ih]+") "+perigee_units[ih]);


      h= book1D("h_perigee_residual_"+perigee_tags[ih],perigee_vars[ih]+" residual",perigee_n[ih],perigee_residual_min[ih],perigee_residual_max[ih]);
      h->SetXTitle("#Delta_{gen-reco}("+perigee_vars[ih]+") "+perigee_units[ih]);

    }
  }

  /*
    if(doWithLocal){
    const TString local_vars[5] = {"q_p","dx_dz","dy_dz","x","y"};
    const TString local_units[5] = {"","[rad]","[rad]","[cm]","[cm]"};
    const TString local_tags[5] = {"q_p","dx_dz","dy_dz","x","y"};
    const int local_n[5] = { 100, 100, 100, 100, 100 };
    const double local_min[5] = { -100, -TMath::Pi()/2., -TMath::Pi()/2., -10, -30 } ;
    const double local_max[5] = { 100, TMath::Pi()/2., TMath::Pi()/2., 10, 30 } ;
    for (uint ihist=0;ihist!=5;++ihist) {
      h= book1D("hist_local_sim_"+local_tags[ihist],local_vars[ihist]+" distribution for SimTrack",local_n[ihist],local_min[ihist],local_max[ihist]);
      h->SetXTitle(local_vars[ihist]+" "+local_units[ihist]);
      h= book1D("hist_local_track_"+local_tags[ihist],local_vars[ihist]+" distribution for Track",local_n[ihist],local_min[ihist],local_max[ihist]);
      h->SetXTitle(local_vars[ihist]+" "+local_units[ihist]);
      //difference
      h= book1D("hist_local_residual_"+local_tags[ihist],local_vars[ihist]+" residual",local_n[ihist],local_min[ihist],local_max[ihist]);
      h->SetXTitle("#Delta_{gen-reco}("+local_vars[ihist]+") "+local_units[ihist]);
      //difference/error
      h= book1D("hist_local_pull_"+local_tags[ihist],local_vars[ihist]+" pull",local_n[ihist],local_min[ihist],local_max[ihist]);
      h->SetXTitle("#frac{#Delta_{gen-reco}}{#sigma}("+local_vars[ihist]+") [100%]");
      //error
      h= book1D("hist_local_sigma_"+local_tags[ihist],local_vars[ihist]+" sigma",local_n[ihist],local_min[ihist],local_max[ihist]);
      h->SetXTitle("#sigma("+local_vars[ihist]+") "+local_units[ihist]);

    }
  }
  */

  if (doWithExtrapolator){
    
    const TString global_FTS_vars[6] = {"x","y","z","p_{x}","p_{y}","p_{z}"};
    const TString global_FTS_units[6] = {"[cm]","[cm]","[cm]","[GeV]","[GeV]","[GeV]"};
    const TString global_FTS_tags[6] = {"x","y","z","px","py","pz"};
    const int global_FTS_n[6] = { 100, 100, 100, 100, 100 ,100};
    const double global_FTS_min[6] = { -6, -6, -6, -200, -200,-200 } ;
    const double global_FTS_max[6] = { 6, 6, 6, 200, 200 , 200 } ;
    const double global_FTS_pull_min[6] = { -6, -6, -6, -6, -6, -6 } ;
    const double global_FTS_pull_max[6] = { 6, 6, 6, 6, 6, 6 } ;
    
    for (uint ihist=0;ihist!=6;++ihist) {
      h= book1D("hist_global_FTS_sim_"+global_FTS_tags[ihist],global_FTS_vars[ihist]+" distribution for SimTrack",global_FTS_n[ihist],global_FTS_min[ihist],global_FTS_max[ihist]);
      h->SetXTitle(global_FTS_vars[ihist]+" "+global_FTS_units[ihist]);
      h= book1D("hist_global_FTS_track_"+global_FTS_tags[ihist],global_FTS_vars[ihist]+" distribution for Track",global_FTS_n[ihist],global_FTS_min[ihist],global_FTS_max[ihist]);
      h->SetXTitle(global_FTS_vars[ihist]+" "+global_FTS_units[ihist]);
      //difference
      h= book1D("hist_global_FTS_residual_"+global_FTS_tags[ihist],global_FTS_vars[ihist]+" residual",global_FTS_n[ihist],global_FTS_min[ihist],global_FTS_max[ihist]);
      h->SetXTitle("#Delta_{gen-reco}("+global_FTS_vars[ihist]+") "+global_FTS_units[ihist]);
      //difference/error
      h= book1D("hist_global_FTS_pull_"+global_FTS_tags[ihist],global_FTS_vars[ihist]+" pull",global_FTS_n[ihist],global_FTS_pull_min[ihist],global_FTS_pull_max[ihist]);
      h->SetXTitle("#frac{#Delta_{gen-reco}}{#sigma}("+global_FTS_vars[ihist]+") [100%]");
      //error
      h= book1D("hist_global_FTS_sigma_"+global_FTS_tags[ihist],global_FTS_vars[ihist]+" sigma",global_FTS_n[ihist],0.,global_FTS_max[ihist]);
      h->SetXTitle("#sigma("+global_FTS_vars[ihist]+") "+global_FTS_units[ihist]);
      
    }
  }

  if(doWithPlane){
    
    const TString local_TSOS_vars[5] = {"q_p","dx_dz","dy_dz","x","y"};
    const TString local_TSOS_units[5] = {"[?]","[rad]","[rad]","[cm]","[cm]"};
    const TString local_TSOS_tags[5] = {"q_p","dx_dz","dy_dz","x","y"};
    const int local_TSOS_n[5] = { 100, 100, 100, 100, 100 };
    const double local_TSOS_min[5] = { -100, -TMath::Pi()/2., -TMath::Pi()/2., -10, -30 } ;
    const double local_TSOS_max[5] = { 100, TMath::Pi()/2., TMath::Pi()/2., 10, 30 } ;
    const double local_TSOS_pull_min[5] = { -6, -6, -6, -6, -6 };
    const double local_TSOS_pull_max[5] = {6, 6, 6, 6, 6 } ;

    for (uint ihist=0;ihist!=5;++ihist) {
      h= book1D("hist_local_TSOS_sim_"+local_TSOS_tags[ihist],local_TSOS_vars[ihist]+" distribution for SimTrack",local_TSOS_n[ihist],local_TSOS_min[ihist],local_TSOS_max[ihist]);
      h->SetXTitle(local_TSOS_vars[ihist]+" "+local_TSOS_units[ihist]);
      h= book1D("hist_local_TSOS_track_"+local_TSOS_tags[ihist],local_TSOS_vars[ihist]+" distribution for Track",local_TSOS_n[ihist],local_TSOS_min[ihist],local_TSOS_max[ihist]);
      h->SetXTitle(local_TSOS_vars[ihist]+" "+local_TSOS_units[ihist]);
      //difference
      h= book1D("hist_local_TSOS_residual_"+local_TSOS_tags[ihist],local_TSOS_vars[ihist]+" residual",local_TSOS_n[ihist],local_TSOS_min[ihist],local_TSOS_max[ihist]);
      h->SetXTitle("#Delta_{gen-reco}("+local_TSOS_vars[ihist]+") "+local_TSOS_units[ihist]);
      //difference/error
      h= book1D("hist_local_TSOS_pull_"+local_TSOS_tags[ihist],local_TSOS_vars[ihist]+" pull",local_TSOS_n[ihist],local_TSOS_pull_min[ihist],local_TSOS_pull_max[ihist]);
      h->SetXTitle("#frac{#Delta_{gen-reco}}{#sigma}("+local_TSOS_vars[ihist]+") [100%]");
      //error
      h= book1D("hist_local_TSOS_sigma_"+local_TSOS_tags[ihist],local_TSOS_vars[ihist]+" sigma",local_TSOS_n[ihist],0.,local_TSOS_max[ihist]);
      h->SetXTitle("#sigma("+local_TSOS_vars[ihist]+") "+local_TSOS_units[ihist]);

    }
  }
  //histogram booked


  edm::LogInfo(theCategory)<<" histograms book.";

}

// ------------ method called once each job just after ending the event loop  ------------
void 
TrackToSimTrackComparator::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackToSimTrackComparator);
