
// -*- C++ -*-
//
// Package:    IsoMuAnalyzer
// Class:      IsoMuAnalyzer
// 
/**\class IsoMuAnalyzer IsoMuAnalyzer.cc Workspace/IsoMuAnalyzer/src/IsoMuAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  "Thomas Danielson"
//         Created:  Thu May  8 12:05:03 CDT 2008
// $Id$
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include <string>
#include <cstdlib>
#include <utility>
#include <vector>
#include <map>
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"

#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/TrackCandidate/interface/TrackCandidate.h>
#include <DataFormats/TrackCandidate/interface/TrackCandidateCollection.h>
#include <DataFormats/TrackingRecHit/interface/TrackingRecHit.h>
#include <DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h>
#include "DataFormats/TrajectoryState/interface/TrackCharge.h"
#include "DataFormats/MuonSeed/interface/L3MuonTrajectorySeed.h"
#include "DataFormats/MuonSeed/interface/L3MuonTrajectorySeedCollection.h"

#include "DataFormats/HLTReco/interface/ModuleTiming.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include <DataFormats/MuonDetId/interface/RPCDetId.h>
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"

#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/PatternTools/interface/TrajectoryStateClosestToBeamLineBuilder.h"
#include "TrackingTools/TrajectoryParametrization/interface/CurvilinearTrajectoryError.h"
#include "TrackingTools/GeomPropagators/interface/StateOnTrackerBound.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"

#include "RecoMuon/TrackingTools/interface/MuonErrorMatrix.h"
#include "RecoMuon/Records/interface/MuonRecoGeometryRecord.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
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
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "HepMC/GenEvent.h"
#include "TLorentzVector.h"

#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"

// L1 Muons
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"
// Higher-level muons
#include <DataFormats/TrackReco/interface/Track.h>
#include "DataFormats/MuonSeed/interface/L2MuonTrajectorySeed.h"
#include "DataFormats/MuonSeed/interface/L2MuonTrajectorySeedCollection.h"
#include "DataFormats/MuonReco/interface/MuonTrackLinks.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

// so that we can use the pt sorting method, we use this:
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "HepPDT/ParticleDataTable.hh"

#include "RecoMuon/MuonXRay/interface/IDconverttoBinNum.h"
#include "RecoMuon/MuonXRay/interface/MotherSearch.h"
#include "RecoMuon/MuonXRay/interface/DQMHelp"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"

//inclusion of muon isolation quantities as seen by the HLT.
#include "PhysicsTools/IsolationAlgos/interface/IsoDepositExtractor.h"
#include "PhysicsTools/IsolationAlgos/interface/IsoDepositExtractorFactory.h"
#include "RecoMuon/L3MuonIsolationProducer/src/L3NominalEfficiencyConfigurator.h"
#include "RecoMuon/MuonIsolation/interface/Cuts.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "RecoMuon/MuonIsolation/interface/Range.h"

#ifdef __CINT__ 
#pragma link C++ class std::map<int,std::vector<double> >++;
#pragma link C++ class std::map<std::string,double>++;
#endif

//
// class decleration
//

class IsoMuAnalyzer : public edm::EDAnalyzer {
public:
  explicit IsoMuAnalyzer(const edm::ParameterSet&);
  ~IsoMuAnalyzer();
  
  
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  // ----------member data ---------------------------
  TFile *theFile; // self-explanatory
  TTree *MuTrigData; // reco only.  What you'd fill for actual data.
  TTree *MuTrigMC; // same, but also includes associations, sim muon information.

  //
  // constants, enums and typedefs
  //

  // Members of both MuTrigData and MuTrigMC
  // Event-level information
  int RunNumber;
  int EventNumber;  
  // Overall MuonHLT execution time, and execution time for each module.
  double totalMuonHLTTime;
  std::map<std::string,double> *muonDigiModuleTimes;
  std::map<std::string,double> *muonLocalRecModuleTimes;
  std::map<std::string,double> *muonL2RecModuleTimes;
  std::map<std::string,double> *muonL3RecModuleTimes;
  std::map<std::string,double> *muonL2IsoModuleTimes;
  std::map<std::string,double> *muonL3IsoModuleTimes;
  std::map<std::string,double> *trackerDigiModuleTimes;
  std::map<std::string,double> *trackerRecModuleTimes;
  std::map<std::string,double> *caloDigiModuleTimes;
  std::map<std::string,double> *caloRecModuleTimes;
  // Trigger information for the single muon iso and non-iso paths
  // Until they aren't, the level 1 are identical for both...
  int l1SingleMuNonIsoTriggered;
  int l2SingleMuNonIsoTriggered;
  int l3SingleMuNonIsoTriggered;
  int l1SingleMuIsoTriggered;
  int l2SingleMuIsoPreTriggered;
  int l2SingleMuIsoTriggered;
  int l3SingleMuIsoPreTriggered;
  int l3SingleMuIsoTriggered;
  int l1DiMuNonIsoTriggered;
  int l2DiMuNonIsoTriggered;
  int l3DiMuNonIsoTriggered;
  int l1DiMuIsoTriggered;
  int l2DiMuIsoPreTriggered;
  int l2DiMuIsoTriggered;
  int l3DiMuIsoPreTriggered;
  int l3DiMuIsoTriggered;
  std::vector<int> *triggerDecisions;
  std::vector<std::string> *triggerNames;
  // HLT and L1 muon information
  int nL1; 
  int nL2;
  int nL3;
  int nL3TracksFromL2;
  int nL3Cands;
  int nL3Seeds;
  // basic L3 muon quantities
  std::vector<double> *l3P;
  std::vector<double> *l3Px;
  std::vector<double> *l3Py;
  std::vector<double> *l3Pz;
  std::vector<double> *l3Pt;
  std::vector<double> *l3PtError;
  std::vector<double> *l3Pt90;
  std::vector<double> *l3Eta;
  std::vector<double> *l3EtaError;
  std::vector<double> *l3Phi;
  std::vector<double> *l3PhiError;
  std::vector<double> *l3D0;
  std::vector<double> *l3D0Error;
  std::vector<int> *l3NHits;
  std::vector<double> *l3Charge;
  std::vector<double> *l3Chi2;
  std::vector<double> *l3Ndof;
  std::map<int,std::vector<int> > *l3DetIds;
  std::map<int,std::vector<int> > *l3SubdetIds;
  std::map<int,std::vector<int> > *l3Component;  
  std::map<int,std::vector<int> > *l3RecHitsStatus;
  std::map<int,std::vector<double> > *l3RecHitsX;
  std::map<int,std::vector<double> > *l3RecHitsY;
  std::map<int,std::vector<double> > *l3RecHitsZ;
  std::map<int, int> *l3NMuHits;
  std::map<int,std::vector<int> > *l3MuStationNumber;
  // L3 muon isolation quantities
  std::vector<double> *l3CalIsoDeposit;
  std::vector<double> *l3TrackIsoDeposit;
  std::vector<double> *l3IsoTrackDR;
  std::vector<double> *l3IsoTrackDRMinDelPt;
  // L3 track fitting and error matrix: note that phi is already in there...
  std::vector<double> *l3Dsz;
  std::vector<double> *l3DszError;
  std::vector<double> *l3Dxy;
  std::vector<double> *l3DxyError;
  std::vector<double> *l3Lambda;
  std::vector<double> *l3LambdaError;
  std::vector<double> *l3Qoverp;
  std::vector<double> *l3QoverpError;
  std::vector<reco::TrackBase::CovarianceMatrix> *l3ErrorMatrix;

  // L2 <-> L3 interface

  // L3 seeding from L2
  std::vector<int> *indexL2SeedingL3;
  std::vector<int> *indexL3SeededFromL2;
  std::vector<int> *l2SeedsL3;

  // Error matrix that's rescaled for the OI algos
  std::map<int,std::vector<double> > *muonErrorMatrix;

  // L3 tracks within window described by muonErrorMatrix

  std::vector<double> *l3TrackP;
  std::vector<double> *l3TrackPx;
  std::vector<double> *l3TrackPy;
  std::vector<double> *l3TrackPz;
  std::vector<double> *l3TrackPt;
  std::vector<double> *l3TrackPtError;
  std::vector<double> *l3TrackEta;
  std::vector<double> *l3TrackEtaError;
  std::vector<double> *l3TrackPhi;
  std::vector<double> *l3TrackPhiError;
  std::vector<double> *l3TrackD0;
  std::vector<double> *l3TrackD0Error;
  std::vector<int> *l3TrackNHits;
  std::vector<double> *l3TrackCharge;
  std::vector<double> *l3TrackChi2;
  std::vector<double> *l3TrackNdof;
  std::map<int,std::vector<int> > *l3TrackDetIds;
  std::map<int,std::vector<int> > *l3TrackSubdetIds;
  std::map<int,std::vector<int> > *l3TrackRecHitsStatus;
  std::map<int,std::vector<double> > *l3TrackRecHitsX;
  std::map<int,std::vector<double> > *l3TrackRecHitsY;
  std::map<int,std::vector<double> > *l3TrackRecHitsZ;
  // L3 track fitting and error matrix: note that phi is already in there...
  std::vector<double> *l3TrackDsz;
  std::vector<double> *l3TrackDszError;
  std::vector<double> *l3TrackDxy;
  std::vector<double> *l3TrackDxyError;
  std::vector<double> *l3TrackLambda;
  std::vector<double> *l3TrackLambdaError;
  std::vector<double> *l3TrackQoverp;
  std::vector<double> *l3TrackQoverpError;
  std::vector<reco::TrackBase::CovarianceMatrix> *l3TrackErrorMatrix;
  
  // basic L2 muon quantities
  std::vector<double> *l2P;
  std::vector<double> *l2Px;
  std::vector<double> *l2Py;
  std::vector<double> *l2Pz;
  std::vector<double> *l2Pt;
  std::vector<double> *l2PtError;
  std::vector<double> *l2Pt90;
  std::vector<double> *l2Eta;
  std::vector<double> *l2EtaError;
  std::vector<double> *l2Phi;
  std::vector<double> *l2PhiError;
  std::vector<double> *l2D0;
  std::vector<double> *l2D0Error;
  std::vector<int> *l2NHits;
  std::vector<double> *l2Charge;
  std::vector<double> *l2Chi2;
  std::vector<double> *l2Ndof;
  std::vector<int> *l2NSeeds;
  std::map<int,std::vector<int> > *l2DetIds;
  std::map<int,std::vector<int> > *l2SubdetIds;
  std::map<int,std::vector<int> > *l2Component;
  std::map<int,std::vector<int> > *l2RecHitsStatus;
  std::map<int,std::vector<double> > *l2RecHitsX;
  std::map<int,std::vector<double> > *l2RecHitsY;
  std::map<int,std::vector<double> > *l2RecHitsZ;
  std::map<int, int> *l2NMuHits;
  std::map<int,std::vector<int> > *l2MuStationNumber;
  // L2 muon isolation quantities
  std::vector<double> *l2CalIsoDeposit;
  // L2 track fitting and error matrix: note that phi is already in there...
  std::vector<double> *l2Dsz;
  std::vector<double> *l2DszError;
  std::vector<double> *l2Dxy;
  std::vector<double> *l2DxyError;
  std::vector<double> *l2Lambda;
  std::vector<double> *l2LambdaError;
  std::vector<double> *l2Qoverp;
  std::vector<double> *l2QoverpError;
  std::vector<reco::TrackBase::CovarianceMatrix> *l2ErrorMatrix;
  // L2 seeding from L1.  Not quite ready to add these yet
  std::vector<int> *indexL1SeedingL2;
  std::vector<int> *indexL2SeededFromL1;
  std::vector<int> *l1SeedsL2;

  // L1 quantities
  std::vector<double> *l1P;
  std::vector<double> *l1Pt;
  std::vector<double> *l1Eta;
  std::vector<double> *l1Phi;
  std::vector<int> *l1Quality;
  std::vector<int> *l1IsIso;
  std::vector<int> *l1IsMip;
  std::vector<int> *l1IsForward;
  std::vector<int> *l1IsRPC;

  // EXCLUSIVE TO MuTrigMC
  std::vector<int> *l3IsAssociated;
  std::vector<int> *l3ParentID;
  std::vector<int> *l3MotherBinNumber;
  std::vector<double> *l3AssociationVar;
  std::vector<double> *l3AssociatedSimMuonPt;
  std::vector<double> *l3AssociatedSimMuonEta;
  std::vector<double> *l3AssociatedSimMuonPhi;
  std::vector<int> *l3AssociatedSimMuonNHits;
  std::map<int,std::vector<int> > *l3AssociatedSimMuonDetIds;
  std::map<int,int> *l3AssociatedSimMuonNMuHits;
  std::map<int,std::vector<int> > *l3AssociatedSimMuonMuStationNumber;

  std::vector<double> *l3AssociatedSimMuonDsz;
  std::vector<double> *l3AssociatedSimMuonDxy;
  std::vector<double> *l3AssociatedSimMuonLambda;
  std::vector<double> *l3AssociatedSimMuonQoverP;

  std::vector<int> *l3TrackIsAssociated;
  std::vector<int> *l3TrackParentID;
  std::vector<int> *l3TrackMotherBinNumber;
  std::vector<double> *l3TrackAssociationVar;
  std::vector<double> *l3TrackAssociatedSimMuonPt;
  std::vector<double> *l3TrackAssociatedSimMuonEta;
  std::vector<double> *l3TrackAssociatedSimMuonPhi;
  std::vector<int> *l3TrackAssociatedSimMuonNHits;
  std::map<int,std::vector<int> > *l3TrackAssociatedSimMuonDetIds;

  std::vector<double> *l3TrackAssociatedSimMuonDsz;
  std::vector<double> *l3TrackAssociatedSimMuonDxy;
  std::vector<double> *l3TrackAssociatedSimMuonLambda;
  std::vector<double> *l3TrackAssociatedSimMuonQoverP;

  std::vector<int> *l2IsAssociated;
  std::vector<int> *l2ParentID;
  std::vector<int> *l2MotherBinNumber;
  std::vector<double> *l2AssociationVar;
  std::vector<double> *l2AssociatedSimMuonPt;
  std::vector<double> *l2AssociatedSimMuonEta;
  std::vector<double> *l2AssociatedSimMuonPhi;
  std::vector<int> *l2AssociatedSimMuonNHits;
  std::map<int,std::vector<int> > *l2AssociatedSimMuonDetIds;
  std::map<int,int> *l2AssociatedSimMuonNMuHits;
  std::map<int,std::vector<int> > *l2AssociatedSimMuonMuStationNumber;

  std::vector<double> *l2AssociatedSimMuonDsz;
  std::vector<double> *l2AssociatedSimMuonDxy;
  std::vector<double> *l2AssociatedSimMuonLambda;
  std::vector<double> *l2AssociatedSimMuonQoverP;

  int nSimMuon;
  std::vector<int> *simMuonParentID;
  std::vector<int> *simMuonMotherBinNumber;
  std::vector<double> *simMuonPt;
  std::vector<double> *simMuonEta;
  std::vector<double> *simMuonPhi;
  std::vector<int> *simMuonNHits;
  std::map<int,std::vector<int> > *simMuonDetIds;
  std::map<int,int> *simMuonNMuHits;
  std::map<int,std::vector<int> > *simMuonMuStationNumber;

  std::vector<double> *simMuonDsz;
  std::vector<double> *simMuonDxy;
  std::vector<double> *simMuonLambda;
  std::vector<double> *simMuonQoverP;

  // SimToReco Associations
  std::vector<int> *simToL3Associated;
  std::vector<double> *simToL3AssociationVar;
  std::vector<int> *simToL3RecoIndex;
  std::vector<int> *simToTkAssociated;
  std::vector<double> *simToTkAssociationVar;
  std::vector<int> *simToTkRecoIndex;
  std::vector<int> *simToL2Associated;
  std::vector<double> *simToL2AssociationVar;
  std::vector<int> *simToL2RecoIndex;

  // And that ends the tiems going into the tree itself.  Now we just need the things 
  // for putting the items into the tree.
  
  // Here's where we get the execution time for the modules.

  edm::InputTag theTimerLabel;
  std::vector<std::string> theMuonDigiModules;
  std::vector<std::string> theTrackerDigiModules;
  std::vector<std::string> theTrackerRecModules;
  std::vector<std::string> theTrackerTrackModules;
  std::vector<std::string> theCaloDigiModules;
  std::vector<std::string> theCaloRecModules;
  std::vector<std::string> theMuonLocalRecModules;
  std::vector<std::string> theMuonL2RecModules;
  std::vector<std::string> theMuonL2IsoModules;
  std::vector<std::string> theMuonL3RecModules;
  std::vector<std::string> theMuonL3IsoModules;
  
  std::string theCategory;
  std::string singleMuIsoTriggerName;
  std::string singleMuNonIsoTriggerName;
  std::string diMuIsoTriggerName;
  std::string diMuNonIsoTriggerName; 

  edm::InputTag l1Label;
  edm::InputTag l2Label;
  edm::InputTag l3Label;
  edm::InputTag trackLabel;
  edm::InputTag candLabel;
  edm::InputTag triggerResults_;
  edm::InputTag trackingParticleLabel;
  // To get the seeds for the L2 muons
  edm::InputTag l2SeedCollectionLabel;
  // To get the seeds for the L3 muons
  edm::InputTag l3SeedCollectionLabel;

  // The all-important error matrix
  edm::ParameterSet errorMatrixPset;
  MuonErrorMatrix * theErrorMatrix;
  std::string thePropagatorName;
  edm::ParameterSet muonServiceParams;
  bool theAdjustAtIp;

  std::string l2AssocLabel;
  std::string l3AssocLabel;
  std::string tkAssocLabel;
  edm::ESHandle<TrackAssociatorBase> l2Associator;
  edm::ESHandle<TrackAssociatorBase> l3Associator;
  edm::ESHandle<TrackAssociatorBase> tkAssociator;

  edm::InputTag theLinkLabel;

  edm::ESHandle<MagneticField> field;
  edm::ESHandle<HepPDT::ParticleDataTable> fPDGTable;
  edm::InputTag bsSrc;

  IDconverttoBinNum wantMotherBin;

  //Needed for the ISO quantities
  reco::isodeposit::IsoDepositExtractor* caloDepositExtractor;
  reco::isodeposit::IsoDepositExtractor* trackDepositExtractor;
  std::string calExtractorName;
  std::string L3IsoTrackCollectionName;
  edm::ParameterSet calExtractorPSet;
  edm::ParameterSet trackExtractorPSet;
  edm::ParameterSet trackCutsPSet;
  edm::InputTag CaloTowerCollectionLabel;
  // This track collection is the tracks for the L3iso track cut
  edm::InputTag inputTrackCollection;
  muonisolation::Cuts L2IsoCalCuts;
  muonisolation::Cuts L3IsoTrackCuts;

};

//
// static data member definitions
//

double pt90(const reco::TrackRef & tk, const edm::Event & ev){
  std::string prov=ev.getProvenance(tk.id()).moduleLabel();
  std::string instance = ev.getProvenance(tk.id()).productInstanceName();

  double nsigma_Pt=0;

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

  return ptLx;
}

//
// constructors and destructor
//
IsoMuAnalyzer::IsoMuAnalyzer(const edm::ParameterSet& iConfig):
  wantMotherBin(iConfig.getParameter<edm::ParameterSet>("IDconverttoBinNum"))
{

  edm::LogInfo("IsoMuAnalyzer") << "into the constructor.";

  l1Label = iConfig.getParameter<edm::InputTag>("l1MuonLabel");
  l2Label = iConfig.getParameter<edm::InputTag>("l2MuonLabel");
  l3Label = iConfig.getParameter<edm::InputTag>("l3MuonLabel");
  trackLabel = iConfig.getParameter<edm::InputTag>("trackLabel");
  candLabel = iConfig.getParameter<edm::InputTag>("candLabel");
  l3SeedCollectionLabel = iConfig.getParameter<edm::InputTag>("l3SeedLabel");
  triggerResults_ = iConfig.getParameter<edm::InputTag>("triggerResults_");
  trackingParticleLabel = iConfig.getParameter<edm::InputTag>("trackingParticleLabel");
  bsSrc = iConfig.getParameter<edm::InputTag>("beamSpotLabel");

  edm::LogInfo("IsoMuAnalyzer") << "got the first block of things.  Now getting errorMatrixPset.";
  errorMatrixPset = iConfig.getParameter<edm::ParameterSet>("matrixPset");
  edm::LogInfo("IsoMuAnalyzer") << "Initializing the matrix rescaler.";
  theErrorMatrix = new MuonErrorMatrix(errorMatrixPset);
  edm::LogInfo("IsoMuAnalyzer") << "Getting muonServiceParameters.";
  muonServiceParams = iConfig.getParameter<edm::ParameterSet>("ServiceParameters");
  edm::LogInfo("IsoMuAnalyzer") << "Getting thePropagatorName.";
  thePropagatorName = iConfig.getParameter<std::string>("propagatorName");
  edm::LogInfo("IsoMuAnalyzer") << "Getting theAdjustAtIP.";
  theAdjustAtIp = errorMatrixPset.getParameter<bool>("atIP");
  edm::LogInfo("IsoMuAnalyzer") << "And onwards...";

  // Get our module timing things

  theTimerLabel=iConfig.getParameter<edm::InputTag>("TimerLabel"); 
  edm::ParameterSet ModulesForTiming =  iConfig.getUntrackedParameter<edm::ParameterSet>("TimingModules");
  theMuonDigiModules=ModulesForTiming.getUntrackedParameter<std::vector<std::string> >("MuonDigiModules");
  theMuonLocalRecModules=ModulesForTiming.getUntrackedParameter<std::vector<std::string> >("MuonLocalRecModules");
  theMuonL2RecModules=ModulesForTiming.getUntrackedParameter<std::vector<std::string> >("MuonL2RecModules");
  theMuonL3RecModules=ModulesForTiming.getUntrackedParameter<std::vector<std::string> >("MuonL3RecModules");
  theMuonL2IsoModules=ModulesForTiming.getUntrackedParameter<std::vector<std::string> >("MuonL2IsoModules");
  theMuonL3IsoModules=ModulesForTiming.getUntrackedParameter<std::vector<std::string> >("MuonL3IsoModules");
  theTrackerDigiModules=ModulesForTiming.getUntrackedParameter<std::vector<std::string> >("TrackerDigiModules");
  theTrackerRecModules=ModulesForTiming.getUntrackedParameter<std::vector<std::string> >("TrackerRecModules");
  theCaloDigiModules=ModulesForTiming.getUntrackedParameter<std::vector<std::string> >("CaloDigiModules");
  theCaloRecModules=ModulesForTiming.getUntrackedParameter<std::vector<std::string> >("CaloRecModules");
  
  // Get our trigger names from the config

  singleMuIsoTriggerName = iConfig.getParameter<std::string>("singleMuIsoTriggerName");
  singleMuNonIsoTriggerName = iConfig.getParameter<std::string>("singleMuNonIsoTriggerName");
  diMuIsoTriggerName = iConfig.getParameter<std::string>("diMuIsoTriggerName");
  diMuNonIsoTriggerName = iConfig.getParameter<std::string>("diMuNonIsoTriggerName"); 

  l2AssocLabel = iConfig.getParameter<std::string>("l2AssociatorName");
  l3AssocLabel = iConfig.getParameter<std::string>("l3AssociatorName");
  tkAssocLabel = iConfig.getParameter<std::string>("tkAssociatorName");

  theLinkLabel = iConfig.getParameter<edm::InputTag>("linkLabel");
  
  // Read in the PSets
  calExtractorPSet = iConfig.getParameter<edm::ParameterSet>("calExtractorPSet");
  trackExtractorPSet = iConfig.getParameter<edm::ParameterSet>("trackExtractorPSet");
  trackCutsPSet = iConfig.getParameter<edm::ParameterSet>("trackCutsPSet");
  // Read in the ISO quantities from L2 (CAL) and L3 (track)
  calExtractorName = calExtractorPSet.getParameter<std::string>("calName");
  L3IsoTrackCollectionName = trackExtractorPSet.getParameter<std::string>("ComponentName");
  // we have our L3 extractor PSet.  Use that to get the pixel tracks.
  inputTrackCollection = trackExtractorPSet.getParameter<edm::InputTag>("inputTrackCollection");
  // Create Extractors to read in the deposits
  caloDepositExtractor = IsoDepositExtractorFactory::get()->create( calExtractorName, calExtractorPSet);
  trackDepositExtractor = IsoDepositExtractorFactory::get()->create(L3IsoTrackCollectionName, trackExtractorPSet);

  // start in on the L3 cuts
  std::string L3IsoTrackCutsName = trackCutsPSet.getParameter<std::string>("L3IsoTrackCutsName");
  if (L3IsoTrackCutsName == "SimpleCuts") {
    L2IsoCalCuts = muonisolation::Cuts(calExtractorPSet);
    L3IsoTrackCuts = muonisolation::Cuts(trackCutsPSet);
  }

  // Initialize things so that we have an address for root to write things to
  triggerDecisions = 0;
  triggerNames = 0;
 
  muonDigiModuleTimes = 0;
  muonLocalRecModuleTimes = 0;
  muonL2RecModuleTimes = 0;
  muonL3RecModuleTimes = 0;
  muonL2IsoModuleTimes = 0;
  muonL3IsoModuleTimes = 0;
  trackerDigiModuleTimes = 0;
  trackerRecModuleTimes = 0;
  caloDigiModuleTimes = 0;
  caloRecModuleTimes = 0;
  
  l3P = 0;
  l3Pt = 0;
  l3Px = 0;
  l3Py = 0;
  l3Pz = 0;
  l3PtError = 0;
  l3Pt90 = 0;
  l3Eta = 0;
  l3EtaError = 0;
  l3Phi = 0;
  l3PhiError = 0;
  l3D0 = 0;
  l3D0Error = 0;
  l3NHits = 0;
  l3Charge = 0;
  l3Chi2 = 0;
  l3Ndof = 0;
  l3DetIds = new std::map<int,std::vector<int> >;
  l3SubdetIds = new std::map<int,std::vector<int> >;
  l3NMuHits = new std::map<int,int>;
  l3MuStationNumber = new std::map<int, std::vector <int> >;
  l3Component = new std::map<int,std::vector<int> >;
  l3RecHitsStatus = new std::map<int,std::vector<int> >;
  l3RecHitsX = new std::map<int,std::vector<double> >;
  l3RecHitsY = new std::map<int,std::vector<double> >;
  l3RecHitsZ = new std::map<int,std::vector<double> >;
  
  l3CalIsoDeposit = 0;
  l3TrackIsoDeposit = 0;
  l3IsoTrackDR = 0;
  l3IsoTrackDRMinDelPt = 0;
  
  l3Dsz = 0;
  l3DszError = 0;
  l3Dxy = 0;
  l3DxyError = 0;
  l3Lambda = 0;
  l3LambdaError = 0;
  l3Qoverp = 0;
  l3QoverpError = 0;
  l3ErrorMatrix = 0;
  
  indexL2SeedingL3 = 0;
  indexL3SeededFromL2 = 0;
  l2SeedsL3 = 0;

  muonErrorMatrix = new std::map<int,std::vector<double> >;

  l3IsAssociated = 0;
  l3ParentID = 0;
  l3MotherBinNumber = 0;
  l3AssociationVar = 0;
  l3AssociatedSimMuonPt = 0;
  l3AssociatedSimMuonEta = 0;
  l3AssociatedSimMuonPhi = 0;
  l3AssociatedSimMuonNHits = 0;
  l3AssociatedSimMuonDetIds = new std::map<int,std::vector<int> >;
  l3AssociatedSimMuonNMuHits = new std::map<int, int>;
  l3AssociatedSimMuonMuStationNumber = new std::map<int, std::vector<int> >;

  l3AssociatedSimMuonDsz = 0;
  l3AssociatedSimMuonDxy = 0;
  l3AssociatedSimMuonLambda = 0;
  l3AssociatedSimMuonQoverP = 0;

  l3TrackP = 0;
  l3TrackPt = 0;
  l3TrackPx = 0;
  l3TrackPy = 0;
  l3TrackPz = 0;
  l3TrackPtError = 0;
  l3TrackEta = 0;
  l3TrackEtaError = 0;
  l3TrackPhi = 0;
  l3TrackPhiError = 0;
  l3TrackD0 = 0;
  l3TrackD0Error = 0;
  l3TrackNHits = 0;
  l3TrackCharge = 0;
  l3TrackChi2 = 0;
  l3TrackNdof = 0;
  l3TrackDetIds = new std::map<int,std::vector<int> >;
  l3TrackSubdetIds = new std::map<int,std::vector<int> >;
  l3TrackRecHitsStatus = new std::map<int,std::vector<int> >;
  l3TrackRecHitsX = new std::map<int,std::vector<double> >;
  l3TrackRecHitsY = new std::map<int,std::vector<double> >;
  l3TrackRecHitsZ = new std::map<int,std::vector<double> >;

  l3TrackDsz = 0;
  l3TrackDszError = 0;
  l3TrackDxy = 0;
  l3TrackDxyError = 0;
  l3TrackLambda = 0;
  l3TrackLambdaError = 0;
  l3TrackQoverp = 0;
  l3TrackQoverpError = 0;
  l3TrackErrorMatrix = 0;

  l3TrackIsAssociated = 0;
  l3TrackParentID = 0;
  l3TrackMotherBinNumber = 0;
  l3TrackAssociationVar = 0;
  l3TrackAssociatedSimMuonPt = 0;
  l3TrackAssociatedSimMuonEta = 0;
  l3TrackAssociatedSimMuonPhi = 0;
  l3TrackAssociatedSimMuonNHits = 0;
  l3TrackAssociatedSimMuonDetIds = new std::map<int,std::vector<int> >;

  l3TrackAssociatedSimMuonDsz = 0;
  l3TrackAssociatedSimMuonDxy = 0;
  l3TrackAssociatedSimMuonLambda = 0;
  l3TrackAssociatedSimMuonQoverP = 0;

  l2P = 0;
  l2Px = 0;
  l2Py = 0;
  l2Pz = 0;
  l2Pt = 0;
  l2PtError = 0;
  l2Pt90 = 0;
  l2Eta = 0;
  l2EtaError = 0;
  l2Phi = 0;
  l2PhiError = 0;
  l2D0 = 0;
  l2D0Error = 0;
  l2NHits = 0;
  l2Charge = 0;
  l2Chi2 = 0;
  l2Ndof = 0;
  l2NSeeds = 0;
  l2DetIds = new std::map<int,std::vector<int> >;
  l2SubdetIds = new std::map<int,std::vector<int> >;
  l2Component = new std::map<int,std::vector<int> >;
  l2NMuHits = new std::map<int,int>;
  l2MuStationNumber = new std::map<int, std::vector <int> >;
  l2RecHitsStatus = new std::map<int,std::vector<int> >;
  l2RecHitsX = new std::map<int,std::vector<double> >;
  l2RecHitsY = new std::map<int,std::vector<double> >;
  l2RecHitsZ = new std::map<int,std::vector<double> >;

  l2CalIsoDeposit = 0;

  l2Dsz = 0;
  l2DszError = 0;
  l2Dxy = 0;
  l2DxyError = 0;
  l2Lambda = 0;
  l2LambdaError = 0;
  l2Qoverp = 0;
  l2QoverpError = 0;
  l2ErrorMatrix = 0;

  indexL1SeedingL2 = 0;
  indexL2SeededFromL1 = 0;
  l1SeedsL2 = 0;
  
  l2IsAssociated = 0;
  l2ParentID = 0;
  l2MotherBinNumber = 0;
  l2AssociationVar = 0;
  l2AssociatedSimMuonPt = 0;
  l2AssociatedSimMuonEta = 0;
  l2AssociatedSimMuonPhi = 0;
  l2AssociatedSimMuonNHits = 0;
  l2AssociatedSimMuonDetIds = new std::map<int,std::vector<int> >;
  l2AssociatedSimMuonNMuHits = new std::map<int, int>;
  l2AssociatedSimMuonMuStationNumber = new std::map<int, std::vector<int> >;

  l2AssociatedSimMuonDsz = 0;
  l2AssociatedSimMuonDxy = 0;
  l2AssociatedSimMuonLambda = 0;
  l2AssociatedSimMuonQoverP = 0;

  simMuonParentID = 0;
  simMuonMotherBinNumber = 0;
  simMuonPt = 0;
  simMuonEta = 0;
  simMuonPhi = 0;
  simMuonNHits = 0;
  simMuonDetIds = new std::map<int,std::vector<int> >;
  simMuonNMuHits = new std::map<int, int>;
  simMuonMuStationNumber = new std::map<int, std::vector<int> >;
    

  simMuonDsz = 0;
  simMuonDxy = 0;
  simMuonLambda = 0;
  simMuonQoverP = 0;

  simToL3Associated = 0;
  simToL3AssociationVar = 0;
  simToL3RecoIndex = 0;
  simToTkAssociated = 0;
  simToTkAssociationVar = 0;
  simToTkRecoIndex = 0;
  simToL2Associated = 0;
  simToL2AssociationVar = 0;
  simToL2RecoIndex = 0;

  l1P = 0;
  l1Pt = 0;
  l1Eta = 0;
  l1Phi = 0;
  l1Quality = 0;
  l1IsIso = 0;
  l1IsMip = 0;
  l1IsForward = 0;
  l1IsRPC = 0;

}


IsoMuAnalyzer::~IsoMuAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void IsoMuAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  edm::LogInfo("IsoMuAnalyzer") << "Beginning of the loop.";

  //get the mag field and the beamspot
  iSetup.get<IdealMagneticFieldRecord>().get(field);
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByLabel(bsSrc,recoBeamSpotHandle);
  reco::BeamSpot bs = *recoBeamSpotHandle;      

  edm::Handle<EventTime> evtTime;
  iEvent.getByLabel(theTimerLabel, evtTime); 

  //get the collection of muons at all levels
  edm::Handle<l1extra::L1MuonParticleCollection> l1Muons;
  edm::Handle<reco::TrackCollection> l2Muons;
  edm::Handle<reco::TrackCollection> l3Muons;
  edm::Handle<reco::TrackCollection> l3MuonTracks;
  edm::Handle<TrackCandidateCollection> l3MuonCands;
  edm::Handle<L3MuonTrajectorySeedCollection> l3Seeds; 

  iEvent.getByLabel(l1Label,l1Muons);
  iEvent.getByLabel(l2Label,l2Muons);
  iEvent.getByLabel(l3Label,l3Muons);
  iEvent.getByLabel(trackLabel,l3MuonTracks);
  iEvent.getByLabel(candLabel,l3MuonCands);
  iEvent.getByLabel(l3SeedCollectionLabel,l3Seeds);

  edm::Handle<edm::View<reco::Track> > l2MuonsForAssociation;
  edm::Handle<edm::View<reco::Track> > l3MuonsForAssociation;
  edm::Handle<edm::View<reco::Track> > l3TracksForAssociation;
  iEvent.getByLabel(l2Label,l2MuonsForAssociation);
  iEvent.getByLabel(l3Label,l3MuonsForAssociation);
  iEvent.getByLabel(trackLabel,l3TracksForAssociation);

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
  iSetup.get<TrackAssociatorRecord>().get(l2AssocLabel,l2Associator);
  iSetup.get<TrackAssociatorRecord>().get(l3AssocLabel,l3Associator);
  //  iSetup.get<TrackAssociatorRecord>().get(k3AssocLabel,l3Associator);

  //associate RecoToSim
  reco::RecoToSimCollection l2RecSimColl = l2Associator->associateRecoToSim(l2MuonsForAssociation, TPtracks, &iEvent);
  reco::RecoToSimCollection l3RecSimColl = l3Associator->associateRecoToSim(l3MuonsForAssociation, TPtracks, &iEvent);
  reco::RecoToSimCollection tkRecSimColl = l3Associator->associateRecoToSim(l3MuonsForAssociation, TPtracks, &iEvent);

  //associate SimToReco
  reco::SimToRecoCollection l2SimRecColl = l2Associator->associateSimToReco(l2MuonsForAssociation, TPtracks, &iEvent);
  reco::SimToRecoCollection l3SimRecColl = l3Associator->associateSimToReco(l3MuonsForAssociation, TPtracks, &iEvent);
  reco::SimToRecoCollection tkSimRecColl = l3Associator->associateSimToReco(l3MuonsForAssociation, TPtracks, &iEvent);

  // Event-level information: run, event, and Trigger Table
  EventNumber = iEvent.id().event();
  RunNumber = iEvent.id().run();   

  // Begin entering muon quantities.  Number of muons at each level.
  // These will be the loop limits over L1, L2, and L3.
  nL1 = l1Muons->size();
  nL2 = l2Muons->size();
  nL3 = l3Muons->size();
  if (!l3MuonTracks.failedToGet()) {
    nL3TracksFromL2 = l3MuonTracks->size();
  }
  else {
    nL3TracksFromL2 = -1;
    edm::LogInfo("IsoMuAnalyzer") << "no l3 tracks";
  }
  if (!l3MuonCands.failedToGet()) {
    nL3Cands = l3MuonCands->size();
  }
  else {
    nL3Cands = -1;
    edm::LogInfo("IsoMuAnalyzer") << "we have no cands";
  }
  if (!l3Seeds.failedToGet()) { // do we have L3 seeds?
    nL3Seeds = l3Seeds->size();
  }
  else {
    nL3Seeds = -1;
    edm::LogInfo("IsoMuAnalyzer") << "We have no seeds";
  }

  //  edm::LogInfo("IsoMuAnalyzer") << "How many L1, L2, L3 do we have? " << nL1 << " " << nL2 << " " << nL3;

  //ISO variables go here
  caloDepositExtractor->fillVetos(iEvent,iSetup,*l2Muons);
  trackDepositExtractor->fillVetos(iEvent,iSetup,*l3Muons);

  reco::IsoDeposit::Vetos trackVetos;
  typedef std::vector< std::pair<reco::TrackRef,reco::IsoDeposit> > MuonsWithDeposits;
  MuonsWithDeposits muonsWithDeposits;

  // get hold of the TriggerResults objects
  edm::Handle<TriggerResults> triggerResults;
  iEvent.getByLabel(triggerResults_,triggerResults);
  TriggerNames namesOfTriggers(*triggerResults);

  // Right, that's the setup done.  Now let's put in our execution times.

  unsigned nTotalModules = evtTime->size();
  totalMuonHLTTime = 0;
  for(unsigned int i = 0; i != nTotalModules; ++i){
    std::string module_name = evtTime->name(i);
    for ( unsigned int j = 0; j != theMuonDigiModules.size(); ++j ) {
      if ( theMuonDigiModules[j] == module_name) {
	totalMuonHLTTime+=evtTime->time(i);
	//	(*muonDigiModuleTimes).push_back(evtTime->time(i));
	(*muonDigiModuleTimes).insert(std::make_pair(theMuonDigiModules[j],evtTime->time(i)));
      }
    }
    for ( unsigned int j = 0; j != theMuonLocalRecModules.size(); ++j ) {
      if ( theMuonLocalRecModules[j] == module_name) {
        totalMuonHLTTime+=evtTime->time(i);
	//        (*muonLocalRecModuleTimes).push_back(evtTime->time(i));
        (*muonLocalRecModuleTimes).insert(std::make_pair(theMuonLocalRecModules[j],evtTime->time(i)));
      }
    }
    for ( unsigned int j = 0; j != theMuonL2RecModules.size(); ++j ) {
      if ( theMuonL2RecModules[j] == module_name) {
        totalMuonHLTTime+=evtTime->time(i);
	//        (*muonL2RecModuleTimes).push_back(evtTime->time(i));
	(*muonL2RecModuleTimes).insert(std::make_pair(theMuonL2RecModules[j],evtTime->time(i)));
      }
    }
    for ( unsigned int j = 0; j != theMuonL3RecModules.size(); ++j ) {
      if ( theMuonL3RecModules[j] == module_name) {
        totalMuonHLTTime+=evtTime->time(i);
	//        (*muonL3RecModuleTimes).push_back(evtTime->time(i));
        (*muonL3RecModuleTimes).insert(std::make_pair(theMuonL3RecModules[j],evtTime->time(i)));
      }
    }
    for ( unsigned int j = 0; j != theMuonL2IsoModules.size(); ++j ) {
      if ( theMuonL2IsoModules[j] == module_name) {
        totalMuonHLTTime+=evtTime->time(i);
	//        (*muonL2IsoModuleTimes).push_back(evtTime->time(i));
        (*muonL2IsoModuleTimes).insert(std::make_pair(theMuonL2IsoModules[j],evtTime->time(i)));
      }
    }
    for ( unsigned int j = 0; j != theMuonL3IsoModules.size(); ++j ) {
      if ( theMuonL3IsoModules[j] == module_name) {
        totalMuonHLTTime+=evtTime->time(i);
	//        (*muonL3IsoModuleTimes).push_back(evtTime->time(i));
        (*muonL3IsoModuleTimes).insert(std::make_pair(theMuonL3IsoModules[j],evtTime->time(i)));
      }
    }
    for ( unsigned int j = 0; j != theTrackerDigiModules.size(); ++j ) {
      if ( theTrackerDigiModules[j] == module_name) {
        totalMuonHLTTime+=evtTime->time(i);
	//        (*trackerDigiModuleTimes).push_back(evtTime->time(i));
        (*trackerDigiModuleTimes).insert(std::make_pair(theTrackerDigiModules[j],evtTime->time(i)));
      }
    }
    for ( unsigned int j = 0; j != theTrackerRecModules.size(); ++j ) {
      if ( theTrackerRecModules[j] == module_name) {
        totalMuonHLTTime+=evtTime->time(i);
	//        (*trackerRecModuleTimes).push_back(evtTime->time(i));
        (*trackerRecModuleTimes).insert(std::make_pair(theTrackerRecModules[j],evtTime->time(i)));
      }
    }
    for ( unsigned int j = 0; j != theCaloDigiModules.size(); ++j ) {
      if ( theCaloDigiModules[j] == module_name) {
        totalMuonHLTTime+=evtTime->time(i);
	//        (*caloDigiModuleTimes).push_back(evtTime->time(i));
	(*caloDigiModuleTimes).insert(std::make_pair(theCaloDigiModules[j],evtTime->time(i)));
      }
    }
    for ( unsigned int j = 0; j != theCaloRecModules.size(); ++j ) {
      if ( theCaloRecModules[j] == module_name) {
        totalMuonHLTTime+=evtTime->time(i);
	//        (*caloRecModuleTimes).push_back(evtTime->time(i));
        (*caloRecModuleTimes).insert(std::make_pair(theCaloRecModules[j],evtTime->time(i)));
      }
    }
  }

  unsigned int indexSingleMuIso = namesOfTriggers.triggerIndex(singleMuIsoTriggerName );
  unsigned int indexSingleMuNoIso = namesOfTriggers.triggerIndex(singleMuNonIsoTriggerName);
  unsigned int indexDiMuIso = namesOfTriggers.triggerIndex(diMuIsoTriggerName );
  unsigned int indexDiMuNoIso = namesOfTriggers.triggerIndex(diMuNonIsoTriggerName); 
  //  edm::LogInfo("IsoMuAnalyzer") << "Trying to get the trigger decisions.";
  //  edm::LogInfo("IsoMuAnalyzer") << "Indexes: " << indexSingleMuIso << " " << indexSingleMuNoIso << " " << indexDiMuIso << " " << indexDiMuNoIso;
  if (triggerResults->index(indexSingleMuNoIso) > 6) l1SingleMuNonIsoTriggered = 1;
  else l1SingleMuNonIsoTriggered = 0;
  if (triggerResults->index(indexSingleMuNoIso) > 20) l2SingleMuNonIsoTriggered = 1;
  else l2SingleMuNonIsoTriggered = 0;
  if (triggerResults->state(indexSingleMuNoIso) == 1) l3SingleMuNonIsoTriggered = 1;
  else l3SingleMuNonIsoTriggered = 0;
  if (triggerResults->index(indexSingleMuIso) > 6) l1SingleMuIsoTriggered = 1;
  else l1SingleMuIsoTriggered = 0;
  if (triggerResults->index(indexSingleMuIso) > 20) l2SingleMuIsoPreTriggered = 1;
  else l2SingleMuIsoPreTriggered = 0;
  if (triggerResults->index(indexSingleMuIso) > 34) l2SingleMuIsoTriggered = 1;
  else l2SingleMuIsoTriggered = 0;
  if (triggerResults->index(indexSingleMuIso) > 44) l3SingleMuIsoPreTriggered = 1;
  else l3SingleMuIsoPreTriggered = 0;
  if (triggerResults->state(indexSingleMuIso) == 1) l3SingleMuIsoTriggered = 1;
  else l3SingleMuIsoTriggered = 0;



  if (triggerResults->index(indexDiMuNoIso) > 6) l1DiMuNonIsoTriggered = 1;
  else l1DiMuNonIsoTriggered = 0;
  if (triggerResults->index(indexDiMuNoIso) > 20) l2DiMuNonIsoTriggered = 1;
  else l2DiMuNonIsoTriggered = 0;
  if (triggerResults->state(indexDiMuNoIso) == 1) l3DiMuNonIsoTriggered = 1;
  else l3DiMuNonIsoTriggered = 0;
  if (triggerResults->index(indexDiMuIso) > 6) l1DiMuIsoTriggered = 1;
  else l1DiMuIsoTriggered = 0;
  if (triggerResults->index(indexDiMuIso) > 20) l2DiMuIsoPreTriggered = 1;
  else l2DiMuIsoPreTriggered = 0;
  if (triggerResults->index(indexDiMuIso) > 34) l2DiMuIsoTriggered = 1;
  else l2DiMuIsoTriggered = 0;
  if (triggerResults->index(indexDiMuIso) > 44) l3DiMuIsoPreTriggered = 1;
  else l3DiMuIsoPreTriggered = 0;
  if (triggerResults->state(indexDiMuIso) == 1) l3DiMuIsoTriggered = 1;
  else l3DiMuIsoTriggered = 0;

  for (unsigned int i = 0; i < triggerResults->size(); i++) {
    (*triggerDecisions).push_back(triggerResults->state(i));
    (*triggerNames).push_back(namesOfTriggers.triggerName(i));
  }

  edm::Handle<reco::MuonTrackLinksCollection> l3ToL2Links;
  iEvent.getByLabel(theLinkLabel, l3ToL2Links);

  for(int iL3 = 0;iL3!=nL3;++iL3){
    reco::TrackRef refL3(l3Muons, iL3);
    
    // Fill in basic information of l3Muons
    (*l3P).push_back(refL3->p());
    (*l3Px).push_back(refL3->px());
    (*l3Py).push_back(refL3->py());
    (*l3Pz).push_back(refL3->pz());
    (*l3Pt).push_back(refL3->pt());
    (*l3PtError).push_back(refL3->ptError());
    (*l3Pt90).push_back(pt90(refL3,iEvent));
    (*l3Eta).push_back(refL3->eta());
    (*l3EtaError).push_back(refL3->etaError());
    (*l3Phi).push_back(refL3->phi());
    (*l3PhiError).push_back(refL3->phiError());
    (*l3D0).push_back(refL3->d0());
    (*l3D0Error).push_back(refL3->d0Error());
    (*l3NHits).push_back(refL3->recHitsSize());
    (*l3Charge).push_back(refL3->charge());
    (*l3Chi2).push_back(refL3->chi2());
    (*l3Ndof).push_back(refL3->ndof());
    // Fill in the track fitting parameters (with phi filled in above)
    (*l3Dsz).push_back(refL3->dsz());
    (*l3DszError).push_back(refL3->dszError());
    (*l3Dxy).push_back(refL3->dxy());
    (*l3DxyError).push_back(refL3->dxyError());
    (*l3Lambda).push_back(refL3->lambda());
    (*l3LambdaError).push_back(refL3->lambdaError());
    (*l3Qoverp).push_back(refL3->qoverp());
    (*l3QoverpError).push_back(refL3->qoverpError());
    (*l3ErrorMatrix).push_back(refL3->covariance());
    
    std::vector<int> *idsForThisL3 = new std::vector<int>;
    std::vector<int> *subidsForThisL3 = new std::vector<int>;
    std::vector<int> *detsForThisL3 = new std::vector<int>;
    std::vector<int> *statusForThisL3 = new std::vector<int>;
    std::vector<double> *xForThisL3 = new std::vector<double>;
    std::vector<double> *yForThisL3 = new std::vector<double>;
    std::vector<double> *zForThisL3 = new std::vector<double>;
    std::vector<int> *stationsForThisL3 = new std::vector<int>;
    int nMuHitsForThisL3 = 0;

    edm::ESHandle<TransientTrackingRecHitBuilder> trackBuilder;
    edm::ESHandle<TransientTrackingRecHitBuilder> muonBuilder;
    std::string trackBuilderName = "WithTrackAngle";
    std::string muonBuilderName = "MuonRecHitBuilder";
    iSetup.get<TransientRecHitRecord>().get(trackBuilderName,trackBuilder);
    iSetup.get<TransientRecHitRecord>().get(muonBuilderName,muonBuilder);

    for (trackingRecHit_iterator l3Hit = refL3->recHitsBegin(); l3Hit != refL3->recHitsEnd(); ++l3Hit) {
      (*idsForThisL3).push_back((*l3Hit)->geographicalId().rawId());
      (*subidsForThisL3).push_back((*l3Hit)->geographicalId().subdetId());
      (*detsForThisL3).push_back((*l3Hit)->geographicalId().det());
      (*statusForThisL3).push_back((*l3Hit)->type());
      if ((*l3Hit)->geographicalId().det() == 1) { // Tracker
	TransientTrackingRecHit::RecHitPointer globL3 = trackBuilder->build(&**l3Hit);
	(*xForThisL3).push_back(globL3->globalPosition().x());
	(*yForThisL3).push_back(globL3->globalPosition().y());
	(*zForThisL3).push_back(globL3->globalPosition().z());
      }
      else if ((*l3Hit)->geographicalId().det() == 2) { // Muon System
	nMuHitsForThisL3++;
	TransientTrackingRecHit::RecHitPointer globL3 = muonBuilder->build(&**l3Hit);
	(*xForThisL3).push_back(globL3->globalPosition().x());
	(*yForThisL3).push_back(globL3->globalPosition().y());
	(*zForThisL3).push_back(globL3->globalPosition().z());
	// number station goes here
	if ( (*l3Hit)->geographicalId().subdetId() == 1) { // DT hit 
	  const DTChamberId& id = DTChamberId((*l3Hit)->geographicalId());
	  (*stationsForThisL3).push_back(id.station());
	}
        if ( (*l3Hit)->geographicalId().subdetId() == 2) { // CSC hit
          const CSCDetId& id = CSCDetId((*l3Hit)->geographicalId());
          (*stationsForThisL3).push_back(id.station());
        }
        if ( (*l3Hit)->geographicalId().subdetId() == 3) { // RPC hit
          const RPCDetId& id = RPCDetId((*l3Hit)->geographicalId());
          (*stationsForThisL3).push_back(id.station());
        }
      }
    }
    
    (*l3DetIds).insert(std::make_pair(iL3,*idsForThisL3));    
    (*l3SubdetIds).insert(std::make_pair(iL3,*subidsForThisL3));
    (*l3Component).insert(std::make_pair(iL3,*detsForThisL3));
    (*l3RecHitsStatus).insert(std::make_pair(iL3,*statusForThisL3));
    (*l3NMuHits).insert(std::make_pair(iL3,nMuHitsForThisL3));
    (*l3MuStationNumber).insert(std::make_pair(iL3,*stationsForThisL3));
    (*l3RecHitsX).insert(std::make_pair(iL3,*xForThisL3));
    (*l3RecHitsY).insert(std::make_pair(iL3,*yForThisL3));
    (*l3RecHitsZ).insert(std::make_pair(iL3,*zForThisL3));

    idsForThisL3->clear();
    subidsForThisL3->clear();
    detsForThisL3->clear();
    statusForThisL3->clear();
    stationsForThisL3->clear();
    xForThisL3->clear();
    yForThisL3->clear();
    zForThisL3->clear();
    nMuHitsForThisL3 = 0;

    // Find the correct L2 muon using track links.  Start with borrowed code.
    bool correctTrackLink = false;
    uint matchingLink=0;
    for(uint iLink=0;iLink!=l3ToL2Links->size();iLink++){
      // L3 is treated as a global track.
      if ((*l3ToL2Links)[iLink].globalTrack() == refL3){
        //this L3muon is coming from this L2 track
        correctTrackLink =true;
        matchingLink=iLink;
      }//same track ref
    }//loop over track links
    if(correctTrackLink) {
      reco::TrackRef refL2FromTrackLinks = (*l3ToL2Links)[matchingLink].standAloneTrack();
      // Get match indices here
      for(int iL2 = 0;iL2!=nL2;++iL2){
	reco::TrackRef refL2ToMatch(l2Muons, iL2);
        if (refL2FromTrackLinks == refL2ToMatch) {
          // Fill these so we know which L2 muon seeded this L3 muon
          (*indexL2SeedingL3).push_back(iL2);
          (*indexL3SeededFromL2).push_back(iL3);
        }
      }
    }

    // get the CAL deposits associated with this muon
    reco::IsoDeposit calDeposit = caloDepositExtractor->deposit(iEvent, iSetup, *refL3);
    // cutting for the L3 muon cal isolation (just getting the cuts and vetos)
    muonisolation::Cuts::CutSpec calo_cuts_here = L2IsoCalCuts(refL3->eta());
    // and deposit for the L3 muon cal isolation
    double conesize = calo_cuts_here.conesize;
    (*l3CalIsoDeposit).push_back(calDeposit.depositWithin(conesize));
    // now for the L3 tracking things
    reco::IsoDeposit trackDeposit = trackDepositExtractor->deposit(iEvent, iSetup, *refL3);
    reco::IsoDeposit::const_iterator pixEnd = trackDeposit.end();
    double dRMinDelPt = 1;
    double dRMin = 1;
    double DelPtMin = 999;
    for (reco::IsoDeposit::const_iterator pix = trackDeposit.begin(); pix != pixEnd; ++pix) {
      double tempDelphi = pix.phi() - refL3->phi();
      if (tempDelphi > TMath::Pi()) tempDelphi = (2 * TMath::Pi()) - tempDelphi;
      double tempdR = sqrt(((pix.eta() - refL3->eta()) * (pix.eta() - refL3->eta())) + (tempDelphi * tempDelphi));
      double tempPt = pix.value();
      if (tempdR < dRMin) dRMin = tempdR;
      if (fabs(refL3->pt() - tempPt) < DelPtMin) {
	dRMinDelPt = tempdR;
	DelPtMin = fabs(refL3->pt() - tempPt);
      }
    }
    (*l3IsoTrackDR).push_back(dRMin);
    (*l3IsoTrackDRMinDelPt).push_back(dRMinDelPt);
    trackVetos.push_back(trackDeposit.veto());
    const muonisolation::Cuts::CutSpec & trackCut = L3IsoTrackCuts(refL3->eta());
    // get the Tracking deposits for our muon
    (*l3TrackIsoDeposit).push_back(trackDeposit.depositWithin(trackCut.conesize, trackVetos));

    bool associated = false;
    // With the detector-level things filled, time to start doing the associations to sim
    int sim_index = 0;
    for (reco::RecoToSimCollection::const_iterator findRefL3 = l3RecSimColl.begin(); findRefL3 != l3RecSimColl.end(); ++findRefL3) {
      const edm::RefToBase<reco::Track> & l3RecSimMatch = findRefL3->key;
      if (l3RecSimMatch->pt() == refL3->pt()) {
	associated = true;
	const std::vector<std::pair<TrackingParticleRef,double> > & tp = findRefL3->val;
	const TrackingParticleRef & trp = tp.begin()->first;
	
	(*l3AssociationVar).push_back(tp.begin()->second);
       
	int particle_ID = trp->pdgId();
	//	int myBin = wantMotherBin.GetBinNum(particle_ID);
	
	if(abs(particle_ID) == 13){
	  // put in the associated pt,eta,phi
	  (*l3AssociatedSimMuonPt).push_back(trp->pt());
	  (*l3AssociatedSimMuonEta).push_back(trp->eta());
	  (*l3AssociatedSimMuonPhi).push_back(trp->phi());
	  if (fabs(trp->phi() - refL3->phi()) > 1) {
	    //Note: keeping this in.  This happens sometimes when the associator used is the 
	    //steppingHelixPropagatorAny
	    edm::LogInfo("IsoMuAnalyzer") << "Something's gone wrong here. First our indexes";
	    edm::LogInfo("IsoMuAnalyzer") << "iL3, sim_index = " << iL3 <<" " << sim_index;
	    edm::LogInfo("IsoMuAnalyzer") << "What about recSimMatch vs trp phi?" << l3RecSimMatch->phi() <<" " << trp->phi();
	  }
	  // put in the detIDs for this sim muon
	  std::vector<int> *idsForSimL3 = new std::vector<int>;
	  std::vector<int> *stationsForSim = new std::vector<int>;
	  int simHitCounter = 0;
	  int simMuHitCounter = 0;
	  for (PSimHitContainer::const_iterator l3SimHit = trp->pSimHit_begin(); l3SimHit != trp->pSimHit_end(); ++l3SimHit) {
	    (*idsForSimL3).push_back((*l3SimHit).detUnitId());
	    DetId theDetUnitId(l3SimHit->detUnitId());
	    int detector = theDetUnitId.det();
	    int subdetector = theDetUnitId.subdetId();
	    if (detector == 2) { //Muon system
	      simMuHitCounter ++;
	      if (subdetector == 1) { //DT
		const DTChamberId& id = DTChamberId(l3SimHit->detUnitId());
		(*stationsForSim).push_back(id.station());
	      }
	      if (subdetector == 2) { //CSC
		const CSCDetId& id=CSCDetId(l3SimHit->detUnitId());
		(*stationsForSim).push_back(id.station());
	      }
	      if (subdetector == 3) { //RPC
		const RPCDetId& id = RPCDetId(l3SimHit->detUnitId());
		(*stationsForSim).push_back(id.station());
	      }
	    }
	    simHitCounter++;
	  }
	  (*l3AssociatedSimMuonNHits).push_back(simHitCounter);
	  (*l3AssociatedSimMuonDetIds).insert(std::make_pair(sim_index,*idsForSimL3));
	  (*l3AssociatedSimMuonMuStationNumber).insert(std::make_pair(sim_index,*stationsForSim));
	  (*l3AssociatedSimMuonNMuHits).insert(std::make_pair(sim_index,simMuHitCounter));
	  idsForSimL3->clear();
	  stationsForSim->clear();
	  sim_index++;
	  //---------------------- MOTHERHOOD --------------------------------
	  //find the parent of tracking particle
	  for(TrackingParticle::g4t_iterator isimtk = trp->g4Track_begin();isimtk!=trp->g4Track_end();isimtk++)
	    {
	      LogDebug(theCategory)<<"I am here 1";
	      if(isimtk->type()==13||isimtk->type()==-13)
		{
		  // This is the sim track for this tracking particle.  Time to put in the parameters
		  FreeTrajectoryState 
		    ftsAtProduction(GlobalPoint(trp->vertex().x(),trp->vertex().y(),trp->vertex().z()),
				    GlobalVector(isimtk->momentum().x(),isimtk->momentum().y(),isimtk->momentum().z()),
				    TrackCharge(trp->charge()),
				    field.product());
		  TrajectoryStateClosestToBeamLineBuilder tscblBuilder;
		  TrajectoryStateClosestToBeamLine tsAtClosestApproach = tscblBuilder(ftsAtProduction,bs);//as in TrackProducerAlgorithm
		  GlobalPoint v1 = tsAtClosestApproach.trackStateAtPCA().position();
		  GlobalVector p = tsAtClosestApproach.trackStateAtPCA().momentum();
		  GlobalPoint v(v1.x()-bs.x0(),v1.y()-bs.y0(),v1.z()-bs.z0());
		  
		  double qoverpSim = tsAtClosestApproach.trackStateAtPCA().charge()/p.mag();
		  double lambdaSim = M_PI/2-p.theta();
		  double dxySim    = (-v.x()*sin(p.phi())+v.y()*cos(p.phi()));
		  double dzSim     = v.z() - (v.x()*p.x()+v.y()*p.y())/p.perp() * p.z()/p.perp();
		  
		  (*l3AssociatedSimMuonDsz).push_back(dzSim);
		  (*l3AssociatedSimMuonDxy).push_back(dxySim);
		  (*l3AssociatedSimMuonLambda).push_back(lambdaSim);
		  (*l3AssociatedSimMuonQoverP).push_back(qoverpSim);

		  //calculate mother hood
		  MotherSearch mother(&*isimtk, SimTk, SimVtx, hepmc);
		  //FIXME, use reco::Particle mother.mother();
		  //                double pt,eta,phi;
		  //                int parentID;
		  //                int motherBinNumber;
		  
		  if (mother.IsValid()){
		    if (mother.SimIsValid()){
		      (*l3ParentID).push_back(mother.Sim_mother->type());
		      (*l3MotherBinNumber).push_back(wantMotherBin.GetBinNum(mother.Sim_mother->type()));
		    }
		    else {
		      (*l3ParentID).push_back(mother.Gen_mother->pdg_id());
		      (*l3MotherBinNumber).push_back(wantMotherBin.GetBinNum(mother.Gen_mother->pdg_id()));
		    }
		    //do it once per tracking particle once it succeeds
		    break;
		  }
		  else{
		    // This handles cases when we have an associated sim muon, but "tricky" muon without 
		    // valid parent (e.g. singleMu)
		    (*l3ParentID).push_back(0);
		    (*l3MotherBinNumber).push_back(-999);
		    edm::LogError(theCategory)<<"tricky muon from TrackingParticle.";
		  }
		}//sim track is a muon
	      else{
		(*l3ParentID).push_back(0);
		(*l3MotherBinNumber).push_back(-999);
		edm::LogError(theCategory)<<"the sim track attached to the tracking particle is not a muon.";
	      }
	    }//loop over SimTrack of tracking particle
	}//muon associated
	else{
	  //a reco muon is associated to something else than a muon
	  edm::LogError(theCategory)<<"a reconstructed muon is associated to: "<<particle_ID;
	  (*l3ParentID).push_back(0);
	  (*l3MotherBinNumber).push_back(-999);
	}
      }//track has an association
      else{
	//this track was not associated.
	edm::LogError(theCategory)<<"a reconstructed muon is not associated.";
      }
    }
    if (associated) {
      //      std::cout << "Associated..." << std::endl;
      (*l3IsAssociated).push_back(1);
    }
    else {
      (*l3IsAssociated).push_back(0);
      (*l3AssociationVar).push_back(-999);
      (*l3AssociatedSimMuonPt).push_back(-999);
      (*l3AssociatedSimMuonEta).push_back(-999);
      (*l3AssociatedSimMuonPhi).push_back(-999);
      (*l3AssociatedSimMuonNHits).push_back(-999);
      (*l3AssociatedSimMuonQoverP).push_back(-999);
      (*l3AssociatedSimMuonLambda).push_back(-999);
      (*l3AssociatedSimMuonDxy).push_back(-999);
      (*l3AssociatedSimMuonDsz).push_back(-999);
      (*l3ParentID).push_back(0);
      (*l3MotherBinNumber).push_back(-999);
    }
  } //loop over l3Muons

  // loop now over hltL3TkTracksFromL2

  for(int iTk = 0;iTk < nL3TracksFromL2;++iTk){
    reco::TrackRef refTk(l3MuonTracks, iTk);
    
    // Fill in basic information of l3Muons
    (*l3TrackP).push_back(refTk->p());
    (*l3TrackPx).push_back(refTk->px());
    (*l3TrackPy).push_back(refTk->py());
    (*l3TrackPz).push_back(refTk->pz());
    (*l3TrackPt).push_back(refTk->pt());
    (*l3TrackPtError).push_back(refTk->ptError());
    (*l3TrackEta).push_back(refTk->eta());
    (*l3TrackEtaError).push_back(refTk->etaError());
    (*l3TrackPhi).push_back(refTk->phi());
    (*l3TrackPhiError).push_back(refTk->phiError());
    (*l3TrackD0).push_back(refTk->d0());
    (*l3TrackD0Error).push_back(refTk->d0Error());
    (*l3TrackNHits).push_back(refTk->recHitsSize());
    (*l3TrackCharge).push_back(refTk->charge());
    (*l3TrackChi2).push_back(refTk->chi2());
    (*l3TrackNdof).push_back(refTk->ndof());
    // Fill in the track fitting parameters (with phi filled in above)
    (*l3TrackDsz).push_back(refTk->dsz());
    (*l3TrackDszError).push_back(refTk->dszError());
    (*l3TrackDxy).push_back(refTk->dxy());
    (*l3TrackDxyError).push_back(refTk->dxyError());
    (*l3TrackLambda).push_back(refTk->lambda());
    (*l3TrackLambdaError).push_back(refTk->lambdaError());
    (*l3TrackQoverp).push_back(refTk->qoverp());
    (*l3TrackQoverpError).push_back(refTk->qoverpError());
    (*l3TrackErrorMatrix).push_back(refTk->covariance());
    
    std::vector<int> *idsForThisTk = new std::vector<int>;
    std::vector<int> *subidsForThisTk = new std::vector<int>;
    std::vector<int> *statusForThisTk = new std::vector<int>;
    std::vector<double> *xForThisTk = new std::vector<double>;
    std::vector<double> *yForThisTk = new std::vector<double>;
    std::vector<double> *zForThisTk = new std::vector<double>;
    
    edm::ESHandle<TransientTrackingRecHitBuilder> trackBuilder;
    std::string trackBuilderName = "WithTrackAngle";
    iSetup.get<TransientRecHitRecord>().get(trackBuilderName,trackBuilder);
    
    edm::LogInfo("IsoMuAnalyzer") << "iterating over tracking rechits";
    for (trackingRecHit_iterator tkHit = refTk->recHitsBegin(); tkHit != refTk->recHitsEnd(); ++tkHit) {     
      if ((*tkHit)->isValid()) {
	(*idsForThisTk).push_back((*tkHit)->geographicalId().rawId());
	(*subidsForThisTk).push_back((*tkHit)->geographicalId().subdetId());
	(*statusForThisTk).push_back((*tkHit)->type());
	if ((*tkHit)->geographicalId().det() == 1) {
	  TransientTrackingRecHit::RecHitPointer globTk = trackBuilder->build(&**tkHit);
	  (*xForThisTk).push_back(globTk->globalPosition().x());
	  (*yForThisTk).push_back(globTk->globalPosition().y());
	  (*zForThisTk).push_back(globTk->globalPosition().z());
	}
	else edm::LogInfo("IsoMuAnalyzer") << "rechit found in detector " << (*tkHit)->geographicalId().det();
      }
    }
    edm::LogInfo("IsoMuAnalyzer") << "iteration finished";
    
    (*l3TrackDetIds).insert(std::make_pair(iTk,*idsForThisTk));
    (*l3TrackSubdetIds).insert(std::make_pair(iTk,*subidsForThisTk));
    (*l3TrackRecHitsStatus).insert(std::make_pair(iTk,*statusForThisTk));
    (*l3TrackRecHitsX).insert(std::make_pair(iTk,*xForThisTk));
    (*l3TrackRecHitsY).insert(std::make_pair(iTk,*yForThisTk));
    (*l3TrackRecHitsZ).insert(std::make_pair(iTk,*zForThisTk));
    
    idsForThisTk->clear();
    subidsForThisTk->clear();
    statusForThisTk->clear();
    xForThisTk->clear();
    yForThisTk->clear();
    zForThisTk->clear();
    
    bool associated = false;
    // With the detector-level things filled, time to start doing the associations to sim
    int sim_index = 0;
    for (reco::RecoToSimCollection::const_iterator findRefTk = tkRecSimColl.begin(); findRefTk != tkRecSimColl.end(); ++findRefTk) {
      const edm::RefToBase<reco::Track> & tkRecSimMatch = findRefTk->key;
      if (tkRecSimMatch->pt() == refTk->pt()) {
	associated = true;
	const std::vector<std::pair<TrackingParticleRef,double> > & tp = findRefTk->val;
	const TrackingParticleRef & trp = tp.begin()->first;
	
	(*l3TrackAssociationVar).push_back(tp.begin()->second);
	
	int particle_ID = trp->pdgId();
	//	int myBin = wantMotherBin.GetBinNum(particle_ID);
	
	if(abs(particle_ID) == 13){
	  // put in the associated pt,eta,phi
	  (*l3TrackAssociatedSimMuonPt).push_back(trp->pt());
	  (*l3TrackAssociatedSimMuonEta).push_back(trp->eta());
	  (*l3TrackAssociatedSimMuonPhi).push_back(trp->phi());
	  if (fabs(trp->phi() - refTk->phi()) > 1) {
	    //Note: keeping this in.  This happens sometimes when the associator used is the
	    //steppingHelixPropagatorAny
	    edm::LogInfo("IsoMuAnalyzer") << "Something's gone wrong here. First our indexes";
	    edm::LogInfo("IsoMuAnalyzer") << "iTk, sim_index = " << iTk <<" " << sim_index;
	    edm::LogInfo("IsoMuAnalyzer") << "What about recSimMatch vs trp phi?" << tkRecSimMatch->phi() <<" " << trp->phi();
	  }
	  // put in the detIDs for this sim muon
	  std::vector<int> *idsForSimTk = new std::vector<int>;
	  int simHitCounter = 0;
	  for (PSimHitContainer::const_iterator tkSimHit = trp->pSimHit_begin(); tkSimHit != trp->pSimHit_end(); ++tkSimHit) {
	    (*idsForSimTk).push_back((*tkSimHit).detUnitId());
	    DetId theDetUnitId(tkSimHit->detUnitId());
	    simHitCounter++;
	  }
	  (*l3TrackAssociatedSimMuonNHits).push_back(simHitCounter);
	  (*l3TrackAssociatedSimMuonDetIds).insert(std::make_pair(sim_index,*idsForSimTk));
	  idsForSimTk->clear();
	  sim_index++;
	  //---------------------- MOTHERHOOD --------------------------------
	  //find the parent of tracking particle
	  for(TrackingParticle::g4t_iterator isimtk = trp->g4Track_begin();isimtk!=trp->g4Track_end();isimtk++)
	    {
	      LogDebug(theCategory)<<"I am here 1";
	      if(isimtk->type()==13||isimtk->type()==-13)
		{
		  // This is the sim track for this tracking particle.  Time to put in the parameters
		  FreeTrajectoryState
		    ftsAtProduction(GlobalPoint(trp->vertex().x(),trp->vertex().y(),trp->vertex().z()),
				    GlobalVector(isimtk->momentum().x(),isimtk->momentum().y(),isimtk->momentum().z()),
				    TrackCharge(trp->charge()),
				    field.product());
		  TrajectoryStateClosestToBeamLineBuilder tscblBuilder;
		  TrajectoryStateClosestToBeamLine tsAtClosestApproach = tscblBuilder(ftsAtProduction,bs);//as in TrackProducerAlgorithm
		  GlobalPoint v1 = tsAtClosestApproach.trackStateAtPCA().position();
		  GlobalVector p = tsAtClosestApproach.trackStateAtPCA().momentum();
		  GlobalPoint v(v1.x()-bs.x0(),v1.y()-bs.y0(),v1.z()-bs.z0());
		  
		  double qoverpSim = tsAtClosestApproach.trackStateAtPCA().charge()/p.mag();
		  double lambdaSim = M_PI/2-p.theta();
		  double dxySim    = (-v.x()*sin(p.phi())+v.y()*cos(p.phi()));
		  double dzSim     = v.z() - (v.x()*p.x()+v.y()*p.y())/p.perp() * p.z()/p.perp();
		  
		  (*l3TrackAssociatedSimMuonDsz).push_back(dzSim);
		  (*l3TrackAssociatedSimMuonDxy).push_back(dxySim);
		  (*l3TrackAssociatedSimMuonLambda).push_back(lambdaSim);
		  (*l3TrackAssociatedSimMuonQoverP).push_back(qoverpSim);
		  
		  //calculate mother hood
		  MotherSearch mother(&*isimtk, SimTk, SimVtx, hepmc);
		  //FIXME, use reco::Particle mother.mother();
		  //                double pt,eta,phi;
		  //                int parentID;
		  //                int motherBinNumber;
		  
		  if (mother.IsValid()){
		    if (mother.SimIsValid()){
		      (*l3TrackParentID).push_back(mother.Sim_mother->type());
		      (*l3TrackMotherBinNumber).push_back(wantMotherBin.GetBinNum(mother.Sim_mother->type()));
		    }
		    else {
		      (*l3TrackParentID).push_back(mother.Gen_mother->pdg_id());
		      (*l3TrackMotherBinNumber).push_back(wantMotherBin.GetBinNum(mother.Gen_mother->pdg_id()));
		    }
		    //do it once per tracking particle once it succeeds
		    break;
		  }
		  else{
		    // This handles cases when we have an associated sim muon, but "tricky" muon without
		    // valid parent (e.g. singleMu)
		    (*l3TrackParentID).push_back(0);
		    (*l3TrackMotherBinNumber).push_back(-999);
		    edm::LogError(theCategory)<<"tricky muon from TrackingParticle.";
		  }
		}//sim track is a muon
	      else{
		(*l3TrackParentID).push_back(0);
		(*l3TrackMotherBinNumber).push_back(-999);
		edm::LogError(theCategory)<<"the sim track attached to the tracking particle is not a muon.";
	      }
	    }//loop over SimTrack of tracking particle
	}//muon associated
	else{
	  //a reco muon is associated to something else than a muon
	  edm::LogError(theCategory)<<"a reconstructed muon is associated to: "<<particle_ID;
	  (*l3TrackParentID).push_back(0);
	  (*l3TrackMotherBinNumber).push_back(-999);
	}
      }//track has an association
      else{
	//this track was not associated.
	edm::LogError(theCategory)<<"a reconstructed muon is not associated.";
      }
    }
    if (associated) {
      //      std::cout << "Associated..." << std::endl;
      (*l3TrackIsAssociated).push_back(1);
    }
    else {
      (*l3TrackIsAssociated).push_back(0);
      (*l3TrackAssociationVar).push_back(-999);
      (*l3TrackAssociatedSimMuonPt).push_back(-999);
      (*l3TrackAssociatedSimMuonEta).push_back(-999);
      (*l3TrackAssociatedSimMuonPhi).push_back(-999);
      (*l3TrackAssociatedSimMuonNHits).push_back(-999);
      (*l3TrackAssociatedSimMuonQoverP).push_back(-999);
      (*l3TrackAssociatedSimMuonLambda).push_back(-999);
      (*l3TrackAssociatedSimMuonDxy).push_back(-999);
      (*l3TrackAssociatedSimMuonDsz).push_back(-999);
      (*l3TrackParentID).push_back(0);
      (*l3TrackMotherBinNumber).push_back(-999);
    }
  } // loop over hltL3TkTracksFromL2

  for(int iL2 = 0;iL2!=nL2;++iL2){
    reco::TrackRef refL2(l2Muons, iL2);
    
    // Fill in basic information of l2Muons
    (*l2P).push_back(refL2->p());
    (*l2Px).push_back(refL2->px());
    (*l2Py).push_back(refL2->py());
    (*l2Pz).push_back(refL2->pz());
    (*l2Pt).push_back(refL2->pt());
    (*l2PtError).push_back(refL2->ptError());
    (*l2Pt90).push_back(pt90(refL2,iEvent));
    (*l2Eta).push_back(refL2->eta());
    (*l2EtaError).push_back(refL2->etaError());
    (*l2Phi).push_back(refL2->phi());
    (*l2PhiError).push_back(refL2->phiError());
    (*l2D0).push_back(refL2->d0());
    (*l2D0Error).push_back(refL2->d0Error());
    (*l2NHits).push_back(refL2->recHitsSize());
    (*l2Charge).push_back(refL2->charge());
    (*l2Chi2).push_back(refL2->chi2());
    (*l2Ndof).push_back(refL2->ndof());
    // Fill in the track fitting parameters (with phi filled in above)
    (*l2Dsz).push_back(refL2->dsz());
    (*l2DszError).push_back(refL2->dszError());
    (*l2Dxy).push_back(refL2->dxy());
    (*l2DxyError).push_back(refL2->dxyError());
    (*l2Lambda).push_back(refL2->lambda());
    (*l2LambdaError).push_back(refL2->lambdaError());
    (*l2Qoverp).push_back(refL2->qoverp());
    (*l2QoverpError).push_back(refL2->qoverpError());
    (*l2ErrorMatrix).push_back(refL2->covariance());
    // Filling in THE muon error matrix
    //    edm::LogInfo("IsoMuAnalyzer") << "Trying to fill the muon error matrix.";

    std::vector<double> *matrixValuesForThisL2 = new std::vector<double>;
    std::vector<int> *idsForThisL2 = new std::vector<int>;
    std::vector<int> *subidsForThisL2 = new std::vector<int>;
    std::vector<int> *detsForThisL2 = new std::vector<int>;
    std::vector<int> *statusForThisL2 = new std::vector<int>;
    std::vector<double> *xForThisL2 = new std::vector<double>;
    std::vector<double> *yForThisL2 = new std::vector<double>;
    std::vector<double> *zForThisL2 = new std::vector<double>;
    std::vector<int> *stationsForThisL2 = new std::vector<int>;
    int nMuHitsForThisL2 = 0;

    // We first start with checking for a valid FTS from the IP: i.e. do we have a valid L2 to seed an L3?  
    // Then we rescale the errorMatrix at the IP.

    TrajectoryStateTransform transform; 

    // To get our error matrix things, we need a MuonServiceProxy.  Declared here.

    MuonServiceProxy *proxy = new MuonServiceProxy(muonServiceParams);
    proxy->update(iSetup);

    // And this we need to do because the f.t.s. doesn't take TrackRefs
    const reco::Track *tk = const_cast<const reco::Track *>(&*refL2);

    FreeTrajectoryState fts = transform.initialFreeState(*tk,&*proxy->magneticField());
   
    //rescale the error at IP, but only if I have to.  And I have to.
    // matrixForThisL2 is the rescaling matrix to be applied.
    AlgebraicSymMatrix55 matrixForThisL2 = theErrorMatrix->get(GlobalVector(refL2->px(),refL2->py(),refL2->pz())).matrix();

    if (theErrorMatrix && theAdjustAtIp){ 
      CurvilinearTrajectoryError oMat = fts.curvilinearError();
      CurvilinearTrajectoryError sfMat = theErrorMatrix->get(fts.momentum());//FIXME with position    
      MuonErrorMatrix::multiply(oMat, sfMat);
      fts = FreeTrajectoryState(fts.parameters(),oMat); 
    }
    
    edm::ESHandle<Propagator> testProp = proxy->propagator(thePropagatorName);

    if (testProp.isValid()) {
      
      StateOnTrackerBound onBounds(testProp.product());
      
      TrajectoryStateOnSurface outer = onBounds(fts);
      
      if (theErrorMatrix && !theAdjustAtIp){ 
	CurvilinearTrajectoryError oMat = outer.curvilinearError();
	CurvilinearTrajectoryError sfMat = theErrorMatrix->get(outer.globalMomentum());//FIXME with position
	MuonErrorMatrix::multiply(oMat, sfMat);   
	outer = TrajectoryStateOnSurface(outer.globalParameters(),
					 oMat,
					 outer.surface(),
					 outer.surfaceSide(),
					 outer.weight());
      }
      
      if (outer.isValid()) {
	AlgebraicSymMatrix55 matrixForThisL2 = outer.curvilinearError().matrix();
	for (int i = 0; i < 5; i++) {
	  for (int j = 0; j < 5; j++) {
	    double temp = matrixForThisL2(i,j);
	    (*matrixValuesForThisL2).push_back(temp);
	  }
	}
      }
      (*muonErrorMatrix).insert(std::make_pair(iL2,*matrixValuesForThisL2));
    }
    matrixValuesForThisL2->clear();



    edm::ESHandle<TransientTrackingRecHitBuilder> muonBuilder;
    std::string muonBuilderName = "MuonRecHitBuilder";
    iSetup.get<TransientRecHitRecord>().get(muonBuilderName,muonBuilder);

    for (trackingRecHit_iterator l2Hit = refL2->recHitsBegin(); l2Hit != refL2->recHitsEnd(); ++l2Hit) {
      (*idsForThisL2).push_back((*l2Hit)->geographicalId().rawId());
      (*subidsForThisL2).push_back((*l2Hit)->geographicalId().subdetId());
      (*detsForThisL2).push_back((*l2Hit)->geographicalId().det());
      (*statusForThisL2).push_back((*l2Hit)->type());
      if ((*l2Hit)->geographicalId().det() == 2) { // Muon System
	nMuHitsForThisL2++;
	TransientTrackingRecHit::RecHitPointer globL2 = muonBuilder->build(&**l2Hit);
        (*xForThisL2).push_back(globL2->globalPosition().x());
        (*yForThisL2).push_back(globL2->globalPosition().y());
        (*zForThisL2).push_back(globL2->globalPosition().z());
	if ( (*l2Hit)->geographicalId().subdetId() == 1) { // DT hit
          const DTChamberId& id = DTChamberId((*l2Hit)->geographicalId());
          (*stationsForThisL2).push_back(id.station());
        }
        if ( (*l2Hit)->geographicalId().subdetId() == 2) { // CSC hit
          const CSCDetId& id = CSCDetId((*l2Hit)->geographicalId());
          (*stationsForThisL2).push_back(id.station());
        }
        if ( (*l2Hit)->geographicalId().subdetId() == 3) { // RPC hit
          const RPCDetId& id = RPCDetId((*l2Hit)->geographicalId());
          (*stationsForThisL2).push_back(id.station());
        }
      }
      else edm::LogError(theCategory)<<"Hits for L2 muon outside muon system.";
    }

    (*l2DetIds).insert(std::make_pair(iL2,*idsForThisL2));
    (*l2SubdetIds).insert(std::make_pair(iL2,*subidsForThisL2));
    (*l2Component).insert(std::make_pair(iL2,*detsForThisL2));
    (*l2RecHitsStatus).insert(std::make_pair(iL2,*statusForThisL2));
    (*l2NMuHits).insert(std::make_pair(iL2,nMuHitsForThisL2));
    (*l2MuStationNumber).insert(std::make_pair(iL2,*stationsForThisL2));
    (*l2RecHitsX).insert(std::make_pair(iL2,*xForThisL2));
    (*l2RecHitsY).insert(std::make_pair(iL2,*yForThisL2));
    (*l2RecHitsZ).insert(std::make_pair(iL2,*zForThisL2));

    int l3_seed_counter = 0;
    if (!l3Seeds.failedToGet()) { // did we get our collection?
      for(L3MuonTrajectorySeedCollection::const_iterator l3Seed = l3Seeds->begin(); l3Seed != l3Seeds->end(); ++l3Seed) {
	reco::TrackRef tmpL2 = l3Seed->l2Track();
	if (tmpL2 == refL2) l3_seed_counter++;
      }
    }
    else l3_seed_counter = -1;
    (*l2NSeeds).push_back(l3_seed_counter);

    idsForThisL2->clear();
    subidsForThisL2->clear();
    detsForThisL2->clear();
    statusForThisL2->clear();
    stationsForThisL2->clear();
    xForThisL2->clear();
    yForThisL2->clear();
    zForThisL2->clear();
    nMuHitsForThisL2 = 0;

    // Determine whether this L2 muon seeds L3
    int seedsL3 = 0;
    for (unsigned int i = 0; i < indexL2SeedingL3->size(); i++) {
      if (iL2 == indexL2SeedingL3->at(i)) seedsL3 = 1;
    }
    (*l2SeedsL3).push_back(seedsL3);

    // Find the correct L1 muon using TrajectorySeeds.
    edm::Ref<L2MuonTrajectorySeedCollection> l2SeedRef = refL2->seedRef().castTo<edm::Ref<L2MuonTrajectorySeedCollection> >();
    l1extra::L1MuonParticleRef l1Ref = l2SeedRef->l1Particle();
    int iL1 = 0;
    for(l1extra::L1MuonParticleCollection::const_iterator itL1 = l1Muons->begin(); itL1 != l1Muons->end(); ++itL1) {
      if (itL1->pt() == l1Ref->pt() && itL1->eta() == l1Ref->eta() && itL1->phi() == l1Ref->phi()) {
	(*indexL1SeedingL2).push_back(iL1);
	(*indexL2SeededFromL1).push_back(iL2);
      }
      iL1++;
    }
    
    // get the Calorimeter isolation deposits for this L2 muon
    reco::IsoDeposit calDeposit = caloDepositExtractor->deposit(iEvent, iSetup, *refL2);
    // cutting for the L2 muon isolation
    muonisolation::Cuts::CutSpec calo_cuts_here = L2IsoCalCuts(refL2->eta());
    // and deposit for the L2 muon isolation
    double conesize = calo_cuts_here.conesize;
    (*l2CalIsoDeposit).push_back(calDeposit.depositWithin(conesize));

    // With the detector-level things filled, time to start doing the associations to sim
    bool associated = false;
    // With the detector-level things filled, time to start doing the associations to sim
    for (reco::RecoToSimCollection::const_iterator findRefL2 = l2RecSimColl.begin(); findRefL2 != l2RecSimColl.end(); ++findRefL2) {
      const edm::RefToBase<reco::Track> & l2RecSimMatch = findRefL2->key;
      int sim_index = 0;
      if (l2RecSimMatch->pt() == refL2->pt()) {
	associated = true;
	const std::vector<std::pair<TrackingParticleRef,double> > & tp = findRefL2->val;
	const TrackingParticleRef & trp = tp.begin()->first;

        (*l2AssociationVar).push_back(tp.begin()->second);

	int particle_ID = trp->pdgId();
	//	int myBin = wantMotherBin.GetBinNum(particle_ID);

	if(abs(particle_ID) == 13){
          // put in the associated pt,eta,phi
          (*l2AssociatedSimMuonPt).push_back(trp->pt());
          (*l2AssociatedSimMuonEta).push_back(trp->eta());
          (*l2AssociatedSimMuonPhi).push_back(trp->phi());
          // put in the detIDs for this sim muon
	  std::vector<int> *idsForSimL2 = new std::vector<int>;
	  std::vector<int> *stationsForSim = new std::vector<int>;
          int simHitCounter = 0;
          int simMuHitCounter = 0;
          for (PSimHitContainer::const_iterator l2SimHit = trp->pSimHit_begin(); l2SimHit != trp->pSimHit_end(); ++l2SimHit) {
            (*idsForSimL2).push_back((*l2SimHit).detUnitId());
	    DetId theDetUnitId(l2SimHit->detUnitId());
            int detector = theDetUnitId.det();
            int subdetector = theDetUnitId.subdetId();
            if (detector == 2) { //Muon system
              simMuHitCounter ++;
              if (subdetector == 1) { //DT
                const DTChamberId& id = DTChamberId(l2SimHit->detUnitId());
                (*stationsForSim).push_back(id.station());
              }
              if (subdetector == 2) { //CSC
                const CSCDetId& id=CSCDetId(l2SimHit->detUnitId());
                (*stationsForSim).push_back(id.station());
              }
              if (subdetector == 3) { //RPC
                const RPCDetId& id = RPCDetId(l2SimHit->detUnitId());
                (*stationsForSim).push_back(id.station());
              }
            }
            simHitCounter++;
          }
          (*l2AssociatedSimMuonNHits).push_back(simHitCounter);
          (*l2AssociatedSimMuonDetIds).insert(std::make_pair(sim_index,*idsForSimL2));
	  (*l2AssociatedSimMuonMuStationNumber).insert(std::make_pair(sim_index,*stationsForSim));
          (*l2AssociatedSimMuonNMuHits).insert(std::make_pair(sim_index,simMuHitCounter));
          idsForSimL2->clear();
	  stationsForSim->clear();
          sim_index++;
	  //---------------------- MOTHERHOOD --------------------------------
	  //find the parent of tracking particle
	  for(TrackingParticle::g4t_iterator isimtk = trp->g4Track_begin();isimtk!=trp->g4Track_end();isimtk++)
	    {
	      LogDebug(theCategory)<<"I am here 1";
	      if(isimtk->type()==13||isimtk->type()==-13)
		{
		  // This is the sim track for this tracking particle.  Time to put in the parameters
		  FreeTrajectoryState 
		    ftsAtProduction(GlobalPoint(trp->vertex().x(),trp->vertex().y(),trp->vertex().z()),
				    GlobalVector(isimtk->momentum().x(),isimtk->momentum().y(),isimtk->momentum().z()),
				    TrackCharge(trp->charge()),
				    field.product());
		  TrajectoryStateClosestToBeamLineBuilder tscblBuilder;
		  TrajectoryStateClosestToBeamLine tsAtClosestApproach = tscblBuilder(ftsAtProduction,bs);//as in TrackProducerAlgorithm
		  GlobalPoint v1 = tsAtClosestApproach.trackStateAtPCA().position();
		  GlobalVector p = tsAtClosestApproach.trackStateAtPCA().momentum();
		  GlobalPoint v(v1.x()-bs.x0(),v1.y()-bs.y0(),v1.z()-bs.z0());
		  
		  double qoverpSim = tsAtClosestApproach.trackStateAtPCA().charge()/p.mag();
		  double lambdaSim = M_PI/2-p.theta();
		  double dxySim    = (-v.x()*sin(p.phi())+v.y()*cos(p.phi()));
		  double dzSim     = v.z() - (v.x()*p.x()+v.y()*p.y())/p.perp() * p.z()/p.perp();
		  
		  (*l2AssociatedSimMuonDsz).push_back(dzSim);
		  (*l2AssociatedSimMuonDxy).push_back(dxySim);
		  (*l2AssociatedSimMuonLambda).push_back(lambdaSim);
		  (*l2AssociatedSimMuonQoverP).push_back(qoverpSim);
		  //calculate mother hood
		  MotherSearch mother(&*isimtk, SimTk, SimVtx, hepmc);
		  //FIXME, use reco::Particle mother.mother();
		  //                double pt,eta,phi;
		  //                int parentID;
		  //                int motherBinNumber;
		  if (mother.IsValid()){
		    if (mother.SimIsValid()){
		      (*l2ParentID).push_back(mother.Sim_mother->type());
		      (*l2MotherBinNumber).push_back(wantMotherBin.GetBinNum(mother.Sim_mother->type()));
		    }
		    else {
		      (*l2ParentID).push_back(mother.Gen_mother->pdg_id());
		      (*l2MotherBinNumber).push_back(wantMotherBin.GetBinNum(mother.Gen_mother->pdg_id()));
		    }
		    //do it once per tracking particle once it succeed
		    break;
		  }
		  else{
		    (*l2ParentID).push_back(0);
		    (*l2MotherBinNumber).push_back(-999);
		    edm::LogError(theCategory)<<"tricky muon from TrackingParticle.";
		  }
		}//sim track is a muon
	      else{
		edm::LogError(theCategory)<<"the sim track attached to the tracking particle is not a muon.";
		(*l2ParentID).push_back(0);
		(*l2MotherBinNumber).push_back(-999);
	      }
	    }//loop over SimTrack of tracking particle
	}//muon associated
	else{
	  //a reco muon is associated to something else than a muon
	  edm::LogError(theCategory)<<"a reconstructed muon is associated to: "<<particle_ID;
	  (*l2ParentID).push_back(0);
	  (*l2MotherBinNumber).push_back(-999);
	}
      }//track has an association
      else{
	//this track was not associated.
	edm::LogError(theCategory)<<"a reconstructed muon is not associated.";
      }
    }
    if (associated) (*l2IsAssociated).push_back(1);
    else {
      (*l2IsAssociated).push_back(0);
      (*l2AssociationVar).push_back(-999);
      (*l2AssociatedSimMuonPt).push_back(-999);
      (*l2AssociatedSimMuonEta).push_back(-999);
      (*l2AssociatedSimMuonPhi).push_back(-999);
      (*l2AssociatedSimMuonNHits).push_back(-999);
      (*l2AssociatedSimMuonQoverP).push_back(-999);
      (*l2AssociatedSimMuonLambda).push_back(-999);
      (*l2AssociatedSimMuonDxy).push_back(-999);
      (*l2AssociatedSimMuonDsz).push_back(-999);
      (*l2ParentID).push_back(0);
      (*l2MotherBinNumber).push_back(-999);
    }
  } //loop over l2Muons

  int iL1 = 0;
  for(l1extra::L1MuonParticleCollection::const_iterator itL1 = l1Muons->begin(); itL1 != l1Muons->end(); ++itL1) {
    (*l1P).push_back(itL1->p());
    (*l1Pt).push_back(itL1->pt());
    (*l1Eta).push_back(itL1->eta());
    (*l1Phi).push_back(itL1->phi());
    (*l1Quality).push_back(itL1->gmtMuonCand().quality());
    if (itL1->isIsolated()) (*l1IsIso).push_back(1);
    else (*l1IsIso).push_back(0);
    if (itL1->isMip()) (*l1IsMip).push_back(1);
    else (*l1IsMip).push_back(0);
    if (itL1->isForward()) (*l1IsForward).push_back(1);
    else (*l1IsForward).push_back(0);
    if (itL1->isRPC()) (*l1IsRPC).push_back(1);
    else (*l1IsRPC).push_back(0);
    int seedsL2 = 0;
    for (unsigned int i = 0; i < indexL1SeedingL2->size(); i++) {
      if (iL1 == indexL1SeedingL2->at(i)) seedsL2 = 1;
    }
    (*l1SeedsL2).push_back(seedsL2);
    iL1++;
  } // loop over l1Muons

  // Loop over all tracking particles

  int sim_index = 0;
  nSimMuon = 0;
  //  for (TrackingParticleCollection::const_iterator trp = (*TPtracks).begin();
  //       trp != (*TPtracks).end(); ++trp) {

  //  edm::LogInfo("IsoMuAnalyzer") << "total number of sim particles = " << (*TPtracks).size();
  
  for (unsigned int iSim = 0; iSim != (*TPtracks).size(); iSim++) {
    
    TrackingParticleRef trp(TPtracks, iSim);
    int particle_ID = trp->pdgId();
    //    if (abs(particle_ID) != 13) edm::LogInfo("IsoMuAnalyzer") << "we have a non-muon in the collection.";
    (*simMuonPt).push_back(trp->pt());
    (*simMuonEta).push_back(trp->eta());
    (*simMuonPhi).push_back(trp->phi());
    
    std::vector<int> *idsForSim = new std::vector<int>;
    std::vector<int> *stationsForSim = new std::vector<int>;
    int simHitCounter = 0;
    int simMuHitCounter = 0;
    for (PSimHitContainer::const_iterator simHit = trp->pSimHit_begin(); simHit != trp->pSimHit_end(); ++simHit) {
      (*idsForSim).push_back((*simHit).detUnitId());
      DetId theDetUnitId(simHit->detUnitId());
      int detector = theDetUnitId.det();
      int subdetector = theDetUnitId.subdetId();
      if (detector == 2) { //Muon system
	simMuHitCounter ++;
	if (subdetector == 1) { //DT 
	  const DTChamberId& id = DTChamberId(simHit->detUnitId());
	  (*stationsForSim).push_back(id.station());
	}
	if (subdetector == 2) { //CSC
	  const CSCDetId& id=CSCDetId(simHit->detUnitId());
	  (*stationsForSim).push_back(id.station());
	}
	if (subdetector == 3) { //RPC
	  const RPCDetId& id = RPCDetId(simHit->detUnitId());
	  (*stationsForSim).push_back(id.station());
	}
      }
      simHitCounter++;
    }
    (*simMuonNHits).push_back(simHitCounter);
    (*simMuonDetIds).insert(std::make_pair(sim_index,*idsForSim));
    (*simMuonMuStationNumber).insert(std::make_pair(sim_index,*stationsForSim));
    (*simMuonNMuHits).insert(std::make_pair(sim_index,simMuHitCounter));
    idsForSim->clear();
    stationsForSim->clear();
    
    for(TrackingParticle::g4t_iterator isimtk = trp->g4Track_begin();isimtk!=trp->g4Track_end();isimtk++)  {
      if(isimtk->type()==13||isimtk->type()==-13) {
	// This is the sim track for this tracking particle.  Time to put in the parameters
	FreeTrajectoryState 
	  ftsAtProduction(GlobalPoint(trp->vertex().x(),trp->vertex().y(),trp->vertex().z()),
			  GlobalVector(isimtk->momentum().x(),isimtk->momentum().y(),isimtk->momentum().z()),
			  TrackCharge(trp->charge()),
			  field.product());
	TrajectoryStateClosestToBeamLineBuilder tscblBuilder;
	TrajectoryStateClosestToBeamLine tsAtClosestApproach = tscblBuilder(ftsAtProduction,bs);//as in TrackProducerAlgorithm
	GlobalPoint v1 = tsAtClosestApproach.trackStateAtPCA().position();
	GlobalVector p = tsAtClosestApproach.trackStateAtPCA().momentum();
	GlobalPoint v(v1.x()-bs.x0(),v1.y()-bs.y0(),v1.z()-bs.z0());
	
	double qoverpSim = tsAtClosestApproach.trackStateAtPCA().charge()/p.mag();
	double lambdaSim = M_PI/2-p.theta();
	double dxySim    = (-v.x()*sin(p.phi())+v.y()*cos(p.phi()));
	double dzSim     = v.z() - (v.x()*p.x()+v.y()*p.y())/p.perp() * p.z()/p.perp();
	
	(*simMuonDsz).push_back(dzSim);
	(*simMuonDxy).push_back(dxySim);
	(*simMuonLambda).push_back(lambdaSim);
	(*simMuonQoverP).push_back(qoverpSim);
	//calculate mother hood
	MotherSearch mother(&*isimtk, SimTk, SimVtx, hepmc);
	if (mother.IsValid()){
	  if (mother.SimIsValid()){
	    (*simMuonParentID).push_back(mother.Sim_mother->type());
	    (*simMuonMotherBinNumber).push_back(wantMotherBin.GetBinNum(mother.Sim_mother->type()));
	  } // simIsValid
	  else {
	  (*simMuonParentID).push_back(mother.Gen_mother->pdg_id());
	  (*simMuonMotherBinNumber).push_back(wantMotherBin.GetBinNum(mother.Gen_mother->pdg_id()));
	  } // gen used otherwise
	} // motherIsValid
	else {
	  (*simMuonParentID).push_back(0);
	  (*simMuonMotherBinNumber).push_back(-999);
	}
	//do it once per tracking particle once it succeed
	break;
      } // sim track is muon
    } //loop over g4_iterator
      // SimToReco associations. First the TrackingParticleRef
      // First look to see if there's a SimToRec match at L3
      /* Temporary commentary here:
	 
      We do an on-the-fly SimToReco association.  There's a small number of events where the result 
      is double-counting: one sim can be associated to 2 reco, if there's a low-pt reco within Delta R 
      of the sim.  The size of this effect is about 1/500 for L2, less for L3.
      
      The MultiTrackValidator uses an AssociatorMap to get around this problem.  I suppose the thing to do
      is to look at what it does, and what I can do to get the same results.  To the LXR, then.
      
      */
    
    std::vector<std::pair<RefToBase<reco::Track>, double> > rt;
    if(l3SimRecColl.find(trp) != l3SimRecColl.end()){
      rt = (std::vector<std::pair<RefToBase<reco::Track>, double> >) l3SimRecColl[trp];
      if (rt.size()!=0) { // Association to L3 successful
	(*simToL3Associated).push_back(1);
	(*simToL3AssociationVar).push_back(rt.begin()->second);
	for (int iL3 = 0; iL3 != nL3; iL3++) {
	  if (rt.begin()->first->pt() == (*l3Pt).at(iL3)) {
	    (*simToL3RecoIndex).push_back(iL3);
	  }
	}
      }
      else { // Something went wrong
	edm::LogInfo("IsoMuAnalyzer")<<"rt.size = 0, but l3SimRec finds trp";
	(*simToL3Associated).push_back(0);
	(*simToL3AssociationVar).push_back(-999);
	(*simToL3RecoIndex).push_back(-999);
      }
    }
    else { // Association to L3 unsuccessful
      (*simToL3Associated).push_back(0);
      (*simToL3AssociationVar).push_back(-999);
      (*simToL3RecoIndex).push_back(-999);
    }
    if(tkSimRecColl.find(trp) != tkSimRecColl.end()){
      rt = (std::vector<std::pair<RefToBase<reco::Track>, double> >) tkSimRecColl[trp];
      if (rt.size()!=0) { // Association to TK successful
        (*simToTkAssociated).push_back(1);
        (*simToTkAssociationVar).push_back(rt.begin()->second);
        for (int iL3 = 0; iL3 != nL3; iL3++) {
          if (rt.begin()->first->pt() == (*l3Pt).at(iL3)) {
            (*simToTkRecoIndex).push_back(iL3);
          }
        }
      }
      else { // Something went wrong
	edm::LogInfo("IsoMuAnalyzer")<<"rt.size = 0, but l3SimRec finds trp";
        (*simToTkAssociated).push_back(0);
        (*simToTkAssociationVar).push_back(-999);
        (*simToTkRecoIndex).push_back(-999);
      }
    }
    else { // Association to L3 unsuccessful
      (*simToTkAssociated).push_back(0);
      (*simToTkAssociationVar).push_back(-999);
      (*simToTkRecoIndex).push_back(-999);
    }
    if(l2SimRecColl.find(trp) != l2SimRecColl.end()){
      rt = (std::vector<std::pair<RefToBase<reco::Track>, double> >) l2SimRecColl[trp];
      if (rt.size()!=0) { // Association to L2 successful
	(*simToL2Associated).push_back(1);
	(*simToL2AssociationVar).push_back(rt.begin()->second);
	for (int iL2 = 0; iL2 != nL2; iL2++) {
	  if (rt.begin()->first->pt() == (*l2Pt).at(iL2)) {
	    (*simToL2RecoIndex).push_back(iL2);
	  }
	}
      }
      else { // Something went wrong
	edm::LogInfo("IsoMuAnalyzer")<<"rt.size = 0, but l2SimRec finds trp";
	(*simToL2Associated).push_back(0);
	(*simToL2AssociationVar).push_back(-999);
	(*simToL2RecoIndex).push_back(-999);
      }
    }
    else { // Association to L2 unsuccessful
      (*simToL2Associated).push_back(0);
      (*simToL2AssociationVar).push_back(-999);
      (*simToL2RecoIndex).push_back(-999);
    }
    
    sim_index++;
    
    //    } //trackingParticle is a muon
  }//loop over all trackingparticles
  
  nSimMuon = sim_index;

  MuTrigData->Fill();
  MuTrigMC->Fill();

  triggerDecisions->clear(); 
  triggerNames->clear();

  muonDigiModuleTimes->clear();
  muonLocalRecModuleTimes->clear();
  muonL2RecModuleTimes->clear();
  muonL3RecModuleTimes->clear();
  muonL2IsoModuleTimes->clear();
  muonL3IsoModuleTimes->clear();
  trackerDigiModuleTimes->clear();
  trackerRecModuleTimes->clear();
  caloDigiModuleTimes->clear();
  caloRecModuleTimes->clear();

  l3P->clear();
  l3Px->clear();
  l3Py->clear();
  l3Pz->clear();
  l3Pt->clear();
  l3PtError->clear();
  l3Pt90->clear();
  l3Eta->clear();
  l3EtaError->clear();
  l3Phi->clear();
  l3PhiError->clear();
  l3D0->clear();
  l3D0Error->clear();
  l3NHits->clear();
  l3Charge->clear();
  l3Chi2->clear();
  l3Ndof->clear();
  l3DetIds->clear();
  l3SubdetIds->clear();
  l3Component->clear();
  l3NMuHits->clear();
  l3MuStationNumber->clear();
  l3RecHitsStatus->clear();
  l3RecHitsX->clear();
  l3RecHitsY->clear();
  l3RecHitsZ->clear();

  l3CalIsoDeposit->clear();
  l3TrackIsoDeposit->clear();
  l3IsoTrackDR->clear();
  l3IsoTrackDRMinDelPt->clear();

  l3Dsz->clear();
  l3DszError->clear();
  l3Dxy->clear();
  l3DxyError->clear();
  l3Lambda->clear();
  l3LambdaError->clear();
  l3Qoverp->clear();
  l3QoverpError->clear();
  l3ErrorMatrix->clear();

  l2SeedsL3->clear();
  indexL2SeedingL3->clear();
  indexL3SeededFromL2->clear();

  muonErrorMatrix->clear();

  l3IsAssociated->clear();
  l3ParentID->clear();
  l3MotherBinNumber->clear();
  l3AssociationVar->clear();
  l3AssociatedSimMuonPt->clear();
  l3AssociatedSimMuonEta->clear();
  l3AssociatedSimMuonPhi->clear();
  l3AssociatedSimMuonNHits->clear();
  l3AssociatedSimMuonDetIds->clear();
  l3AssociatedSimMuonNMuHits->clear();
  l3AssociatedSimMuonMuStationNumber->clear();

  l3AssociatedSimMuonDsz->clear();
  l3AssociatedSimMuonDxy->clear();
  l3AssociatedSimMuonLambda->clear();
  l3AssociatedSimMuonQoverP->clear();

  l3TrackP->clear();
  l3TrackPx->clear();
  l3TrackPy->clear();
  l3TrackPz->clear();
  l3TrackPt->clear();
  l3TrackPtError->clear();
  l3TrackEta->clear();
  l3TrackEtaError->clear();
  l3TrackPhi->clear();
  l3TrackPhiError->clear();
  l3TrackD0->clear();
  l3TrackD0Error->clear();
  l3TrackNHits->clear();
  l3TrackCharge->clear();
  l3TrackChi2->clear();
  l3TrackNdof->clear();
  l3TrackDetIds->clear();
  l3TrackSubdetIds->clear();
  l3TrackRecHitsStatus->clear();
  l3TrackRecHitsX->clear();
  l3TrackRecHitsY->clear();
  l3TrackRecHitsZ->clear();

  l3TrackDsz->clear();
  l3TrackDszError->clear();
  l3TrackDxy->clear();
  l3TrackDxyError->clear();
  l3TrackLambda->clear();
  l3TrackLambdaError->clear();
  l3TrackQoverp->clear();
  l3TrackQoverpError->clear();
  l3TrackErrorMatrix->clear();

  l3TrackIsAssociated->clear();
  l3TrackParentID->clear();
  l3TrackMotherBinNumber->clear();
  l3TrackAssociationVar->clear();
  l3TrackAssociatedSimMuonPt->clear();
  l3TrackAssociatedSimMuonEta->clear();
  l3TrackAssociatedSimMuonPhi->clear();
  l3TrackAssociatedSimMuonNHits->clear();
  l3TrackAssociatedSimMuonDetIds->clear();

  l3TrackAssociatedSimMuonDsz->clear();
  l3TrackAssociatedSimMuonDxy->clear();
  l3TrackAssociatedSimMuonLambda->clear();
  l3TrackAssociatedSimMuonQoverP->clear();

  l2P->clear();
  l2Px->clear();
  l2Py->clear();
  l2Pz->clear();
  l2Pt->clear();
  l2PtError->clear();
  l2Pt90->clear();
  l2Eta->clear();
  l2EtaError->clear();
  l2Phi->clear();
  l2PhiError->clear();
  l2D0->clear();
  l2D0Error->clear();
  l2NHits->clear();
  l2Charge->clear();
  l2Chi2->clear();
  l2Ndof->clear();
  l2DetIds->clear();
  l2SubdetIds->clear();
  l2Component->clear();
  l2NMuHits->clear();
  l2NSeeds->clear();
  l2MuStationNumber->clear();
  l2RecHitsStatus->clear();
  l2RecHitsX->clear();
  l2RecHitsY->clear();
  l2RecHitsZ->clear();

  l2CalIsoDeposit->clear();

  l2Dsz->clear();
  l2DszError->clear();
  l2Dxy->clear();
  l2DxyError->clear();
  l2Lambda->clear();
  l2LambdaError->clear();
  l2Qoverp->clear();
  l2QoverpError->clear();
  l2ErrorMatrix->clear();

  l1SeedsL2->clear();
  indexL1SeedingL2->clear();
  indexL2SeededFromL1->clear();

  l2IsAssociated->clear();
  l2ParentID->clear();
  l2MotherBinNumber->clear();
  l2AssociationVar->clear();
  l2AssociatedSimMuonPt->clear();
  l2AssociatedSimMuonEta->clear();
  l2AssociatedSimMuonPhi->clear();
  l2AssociatedSimMuonNHits->clear();
  l2AssociatedSimMuonDetIds->clear();
  l2AssociatedSimMuonNMuHits->clear();
  l2AssociatedSimMuonMuStationNumber->clear();

  l2AssociatedSimMuonDsz->clear();
  l2AssociatedSimMuonDxy->clear();
  l2AssociatedSimMuonLambda->clear();
  l2AssociatedSimMuonQoverP->clear();

  simMuonParentID->clear();
  simMuonMotherBinNumber->clear();
  simMuonPt->clear();
  simMuonEta->clear();
  simMuonPhi->clear();
  simMuonNHits->clear();
  simMuonDetIds->clear();
  simMuonNMuHits->clear();
  simMuonMuStationNumber->clear();

  simMuonDsz->clear();
  simMuonDxy->clear();
  simMuonLambda->clear();
  simMuonQoverP->clear();

  simToL3Associated->clear();
  simToL3AssociationVar->clear();
  simToL3RecoIndex->clear();
  simToTkAssociated->clear();
  simToTkAssociationVar->clear();
  simToTkRecoIndex->clear();
  simToL2Associated->clear();
  simToL2AssociationVar->clear();
  simToL2RecoIndex->clear();

  l1P->clear();
  l1Pt->clear();
  l1Eta->clear();
  l1Phi->clear();
  l1Quality->clear();
  l1IsIso->clear();
  l1IsMip->clear();
  l1IsForward->clear();
  l1IsRPC->clear();

}


// ------------ method called once each job just before starting event loop  ------------
void 
IsoMuAnalyzer::beginJob(const edm::EventSetup&)
{

  theFile = new TFile("HLTMuonTree.root","recreate");
  theFile->cd();
  
  MuTrigData = new TTree("MuTrigData","MuTrigData");
  MuTrigMC = new TTree("MuTrigMC","MuTrigMC");

  // MuTrigData branches

  //Event-level information
  MuTrigData->Branch("RunNumber",&RunNumber,"RunNumber/I");
  MuTrigData->Branch("EventNumber",&EventNumber,"EventNumber/I");
  // Execution times
  MuTrigData->Branch("totalMuonHLTTime",&totalMuonHLTTime,"totalMuonHLTTime/D");
  MuTrigData->Branch("muonDigiModuleTimes",&muonDigiModuleTimes);  
  MuTrigData->Branch("muonLocalRecModuleTimes",&muonLocalRecModuleTimes);
  MuTrigData->Branch("muonL2RecModuleTimes",&muonL2RecModuleTimes);
  MuTrigData->Branch("muonL3RecModuleTimes",&muonL3RecModuleTimes);
  MuTrigData->Branch("muonL2IsoModuleTimes",&muonL2IsoModuleTimes);
  MuTrigData->Branch("muonL3IsoModuleTimes",&muonL3IsoModuleTimes);
  MuTrigData->Branch("trackerDigiModuleTimes",&trackerDigiModuleTimes);
  MuTrigData->Branch("trackerRecModuleTimes",&trackerRecModuleTimes);
  MuTrigData->Branch("caloDigiModuleTimes",&caloDigiModuleTimes);
  MuTrigData->Branch("caloRecModuleTimes",&caloRecModuleTimes);

  //Trigger information
  MuTrigData->Branch("l1SingleMuNonIsoTriggered",&l1SingleMuNonIsoTriggered,"l1SingleMuNonIsoTriggered/I");
  MuTrigData->Branch("l2SingleMuNonIsoTriggered",&l2SingleMuNonIsoTriggered,"l2SingleMuNonIsoTriggered/I");
  MuTrigData->Branch("l3SingleMuNonIsoTriggered",&l3SingleMuNonIsoTriggered,"l3SingleMuNonIsoTriggered/I");
  MuTrigData->Branch("l1SingleMuIsoTriggered",&l1SingleMuIsoTriggered,"l1SingleMuIsoTriggered/I");
  MuTrigData->Branch("l2SingleMuIsoPreTriggered",&l2SingleMuIsoPreTriggered,"l2SingleMuIsoPreTriggered/I");
  MuTrigData->Branch("l2SingleMuIsoTriggered",&l2SingleMuIsoTriggered,"l2SingleMuIsoTriggered/I");
  MuTrigData->Branch("l3SingleMuIsoPreTriggered",&l3SingleMuIsoPreTriggered,"l3SingleMuIsoPreTriggered/I");
  MuTrigData->Branch("l3SingleMuIsoTriggered",&l3SingleMuIsoTriggered,"l3SingleMuIsoTriggered/I");
  MuTrigData->Branch("l1DiMuNonIsoTriggered",&l1DiMuNonIsoTriggered,"l1DiMuNonIsoTriggered/I");
  MuTrigData->Branch("l2DiMuNonIsoTriggered",&l2DiMuNonIsoTriggered,"l2DiMuNonIsoTriggered/I");
  MuTrigData->Branch("l3DiMuNonIsoTriggered",&l3DiMuNonIsoTriggered,"l3DiMuNonIsoTriggered/I");
  MuTrigData->Branch("l1DiMuIsoTriggered",&l1DiMuIsoTriggered,"l1DiMuIsoTriggered/I");
  MuTrigData->Branch("l2DiMuIsoPreTriggered",&l2DiMuIsoPreTriggered,"l2DiMuIsoPreTriggered/I");
  MuTrigData->Branch("l2DiMuIsoTriggered",&l2DiMuIsoTriggered,"l2DiMuIsoTriggered/I");
  MuTrigData->Branch("l3DiMuIsoPreTriggered",&l3DiMuIsoPreTriggered,"l3DiMuIsoPreTriggered/I");
  MuTrigData->Branch("l3DiMuIsoTriggered",&l3DiMuIsoTriggered,"l3DiMuIsoTriggered/I"); 
  MuTrigData->Branch("triggerDecisions",&triggerDecisions);
  MuTrigData->Branch("triggerNames",&triggerNames);
  // number of muons at each level
  MuTrigData->Branch("nL1",&nL1,"nL1/I");
  MuTrigData->Branch("nL2",&nL2,"nL2/I");
  MuTrigData->Branch("nL3",&nL3,"nL3/I");
  MuTrigData->Branch("nL3TracksFromL2",&nL3TracksFromL2,"nL3TracksFromL2/I");
  MuTrigData->Branch("nL3Cands",&nL3Cands,"nL3Cands/I");
  MuTrigData->Branch("nL3Seeds",&nL3Seeds,"nL3Seeds/I");
  // L3 muon information: the basics
  MuTrigData->Branch("l3P",&l3P);
  MuTrigData->Branch("l3Px",&l3Px);
  MuTrigData->Branch("l3Py",&l3Py);
  MuTrigData->Branch("l3Pz",&l3Pz);
  MuTrigData->Branch("l3Pt",&l3Pt);
  MuTrigData->Branch("l3PtError",&l3PtError);
  MuTrigData->Branch("l3Pt90",&l3Pt90);
  MuTrigData->Branch("l3Eta",&l3Eta);
  MuTrigData->Branch("l3EtaError",&l3EtaError);
  MuTrigData->Branch("l3Phi",&l3Phi);
  MuTrigData->Branch("l3PhiError",&l3PhiError);
  MuTrigData->Branch("l3D0",&l3D0);
  MuTrigData->Branch("l3D0Error",&l3D0Error);
  MuTrigData->Branch("l3NHits",&l3NHits);
  MuTrigData->Branch("l3Charge",&l3Charge);
  MuTrigData->Branch("l3Chi2",&l3Chi2);
  MuTrigData->Branch("l3Ndof",&l3Ndof);
  // L3 Muon Isolation quantities
  MuTrigData->Branch("l3CalIsoDeposit",&l3CalIsoDeposit);
  MuTrigData->Branch("l3TrackIsoDeposit",&l3TrackIsoDeposit);
  MuTrigData->Branch("l3IsoTrackDR",&l3IsoTrackDR);
  MuTrigData->Branch("l3IsoTrackDRMinDelPt",&l3IsoTrackDRMinDelPt);
  // L3 Muon Track fitting parameters (with phi already declared above)
  MuTrigData->Branch("l3Dsz",&l3Dsz);
  MuTrigData->Branch("l3DszError",&l3DszError);
  MuTrigData->Branch("l3Dxy",&l3Dxy);
  MuTrigData->Branch("l3DxyError",&l3DxyError);
  MuTrigData->Branch("l3Lambda",&l3Lambda);
  MuTrigData->Branch("l3LambdaError",&l3LambdaError);
  MuTrigData->Branch("l3Qoverp",&l3Qoverp);
  MuTrigData->Branch("l3QoverpError",&l3QoverpError);
  MuTrigData->Branch("l3ErrorMatrix",&l3ErrorMatrix);
  MuTrigData->Branch("l3DetIds",&l3DetIds);
  MuTrigData->Branch("l3SubdetIds",&l3SubdetIds);
  MuTrigData->Branch("l3Component",&l3Component);
  MuTrigData->Branch("l3NMuHits",&l3NMuHits);
  MuTrigData->Branch("l3MuStationNumber",&l3MuStationNumber);
  MuTrigData->Branch("l3RecHitsStatus",&l3RecHitsStatus);
  MuTrigData->Branch("l3RecHitsX",&l3RecHitsX);
  MuTrigData->Branch("l3RecHitsY",&l3RecHitsY);
  MuTrigData->Branch("l3RecHitsZ",&l3RecHitsZ);
  // Indices for L3<->L2
  MuTrigData->Branch("indexL2SeedingL3",&indexL2SeedingL3);
  MuTrigData->Branch("indexL3SeededFromL2",&indexL3SeededFromL2);
  MuTrigData->Branch("l2SeedsL3",&l2SeedsL3);

  // The muon error matrix
  MuTrigData->Branch("muonErrorMatrix",&muonErrorMatrix);

  // information from hltL3TkTracksFromL2
  MuTrigData->Branch("l3TrackP",&l3TrackP);
  MuTrigData->Branch("l3TrackPx",&l3TrackPx);
  MuTrigData->Branch("l3TrackPy",&l3TrackPy);
  MuTrigData->Branch("l3TrackPz",&l3TrackPz);
  MuTrigData->Branch("l3TrackPt",&l3TrackPt);
  MuTrigData->Branch("l3TrackPtError",&l3TrackPtError);
  MuTrigData->Branch("l3TrackEta",&l3TrackEta);
  MuTrigData->Branch("l3TrackEtaError",&l3TrackEtaError);
  MuTrigData->Branch("l3TrackPhi",&l3TrackPhi);
  MuTrigData->Branch("l3TrackPhiError",&l3TrackPhiError);
  MuTrigData->Branch("l3TrackD0",&l3TrackD0);
  MuTrigData->Branch("l3TrackD0Error",&l3TrackD0Error);
  MuTrigData->Branch("l3TrackNHits",&l3TrackNHits);
  MuTrigData->Branch("l3TrackCharge",&l3TrackCharge);
  MuTrigData->Branch("l3TrackChi2",&l3TrackChi2);
  MuTrigData->Branch("l3TrackNdof",&l3TrackNdof);
  // L3TRACK Muon Track fitting parameters (with phi already declared above)
  MuTrigData->Branch("l3TrackDsz",&l3TrackDsz);
  MuTrigData->Branch("l3TrackDszError",&l3TrackDszError);
  MuTrigData->Branch("l3TrackDxy",&l3TrackDxy);
  MuTrigData->Branch("l3TrackDxyError",&l3TrackDxyError);
  MuTrigData->Branch("l3TrackLambda",&l3TrackLambda);
  MuTrigData->Branch("l3TrackLambdaError",&l3TrackLambdaError);
  MuTrigData->Branch("l3TrackQoverp",&l3TrackQoverp);
  MuTrigData->Branch("l3TrackQoverpError",&l3TrackQoverpError);
  MuTrigData->Branch("l3TrackErrorMatrix",&l3TrackErrorMatrix);
  MuTrigData->Branch("l3TrackDetIds",&l3TrackDetIds);
  MuTrigData->Branch("l3TrackSubdetIds",&l3TrackSubdetIds);
  MuTrigData->Branch("l3TrackRecHitsStatus",&l3TrackRecHitsStatus);
  MuTrigData->Branch("l3TrackRecHitsX",&l3TrackRecHitsX);
  MuTrigData->Branch("l3TrackRecHitsY",&l3TrackRecHitsY);
  MuTrigData->Branch("l3TrackRecHitsZ",&l3TrackRecHitsZ);

  // L2 muon information: the basics
  MuTrigData->Branch("l2P",&l2P);
  MuTrigData->Branch("l2Px",&l2Px);
  MuTrigData->Branch("l2Py",&l2Py);
  MuTrigData->Branch("l2Pz",&l2Pz);
  MuTrigData->Branch("l2Pt",&l2Pt);
  MuTrigData->Branch("l2PtError",&l2PtError);
  MuTrigData->Branch("l2Pt90",&l2Pt90);
  MuTrigData->Branch("l2Eta",&l2Eta);
  MuTrigData->Branch("l2EtaError",&l2EtaError);
  MuTrigData->Branch("l2Phi",&l2Phi);
  MuTrigData->Branch("l2PhiError",&l2PhiError);
  MuTrigData->Branch("l2D0",&l2D0);
  MuTrigData->Branch("l2D0Error",&l2D0Error);
  MuTrigData->Branch("l2NHits",&l2NHits);
  MuTrigData->Branch("l2Charge",&l2Charge);
  MuTrigData->Branch("l2Chi2",&l2Chi2);
  MuTrigData->Branch("l2Ndof",&l2Ndof);
  MuTrigData->Branch("l2NSeeds",&l2NSeeds);
  // L2 Muon Isolation quantities
  MuTrigData->Branch("l2CalIsoDeposit",&l2CalIsoDeposit);
  // L2 Muon Track fitting parameters (with phi already declared above)
  MuTrigData->Branch("l2Dsz",&l2Dsz);
  MuTrigData->Branch("l2DszError",&l2DszError);
  MuTrigData->Branch("l2Dxy",&l2Dxy);
  MuTrigData->Branch("l2DxyError",&l2DxyError);
  MuTrigData->Branch("l2Lambda",&l2Lambda);
  MuTrigData->Branch("l2LambdaError",&l2LambdaError);
  MuTrigData->Branch("l2Qoverp",&l2Qoverp);
  MuTrigData->Branch("l2QoverpError",&l2QoverpError);
  MuTrigData->Branch("l2ErrorMatrix",&l2ErrorMatrix);
  MuTrigData->Branch("l2DetIds",&l2DetIds);
  MuTrigData->Branch("l2SubdetIds",&l2SubdetIds);
  MuTrigData->Branch("l2Component",&l2Component);
  MuTrigData->Branch("l2NMuHits",&l2NMuHits);
  MuTrigData->Branch("l2MuStationNumber",&l2MuStationNumber);
  MuTrigData->Branch("l2RecHitsStatus",&l2RecHitsStatus);
  MuTrigData->Branch("l2RecHitsX",&l2RecHitsX);
  MuTrigData->Branch("l2RecHitsY",&l2RecHitsY);
  MuTrigData->Branch("l2RecHitsZ",&l2RecHitsZ);

  MuTrigData->Branch("l1SeedsL2",&l1SeedsL2);
  MuTrigData->Branch("indexL1SeedingL2",&indexL1SeedingL2);
  MuTrigData->Branch("indexL2SeededFromL1",&indexL2SeededFromL1);

  MuTrigData->Branch("l1P",&l1P);
  MuTrigData->Branch("l1Pt",&l1Pt);
  MuTrigData->Branch("l1Eta",&l1Eta);
  MuTrigData->Branch("l1Phi",&l1Phi);
  MuTrigData->Branch("l1Quality",&l1Quality);
  MuTrigData->Branch("l1IsIso",&l1IsIso);
  MuTrigData->Branch("l1IsMip",&l1IsMip);
  MuTrigData->Branch("l1IsForward",&l1IsForward);
  MuTrigData->Branch("l1IsRPC",&l1IsRPC);


  // MuTrigMC branches

  //Event-level information
  MuTrigMC->Branch("RunNumber",&RunNumber,"RunNumber/I");
  MuTrigMC->Branch("EventNumber",&EventNumber,"EventNumber/I");
  // Execution times
  MuTrigMC->Branch("totalMuonHLTTime",&totalMuonHLTTime,"totalMuonHLTTime/D");
  MuTrigMC->Branch("muonDigiModuleTimes",&muonDigiModuleTimes);
  MuTrigMC->Branch("muonLocalRecModuleTimes",&muonLocalRecModuleTimes);
  MuTrigMC->Branch("muonL2RecModuleTimes",&muonL2RecModuleTimes);
  MuTrigMC->Branch("muonL3RecModuleTimes",&muonL3RecModuleTimes);
  MuTrigMC->Branch("muonL2IsoModuleTimes",&muonL2IsoModuleTimes);
  MuTrigMC->Branch("muonL3IsoModuleTimes",&muonL3IsoModuleTimes);
  MuTrigMC->Branch("trackerDigiModuleTimes",&trackerDigiModuleTimes);
  MuTrigMC->Branch("trackerRecModuleTimes",&trackerRecModuleTimes);
  MuTrigMC->Branch("caloDigiModuleTimes",&caloDigiModuleTimes);
  MuTrigMC->Branch("caloRecModuleTimes",&caloRecModuleTimes);
  //Trigger information
  MuTrigMC->Branch("l1SingleMuNonIsoTriggered",&l1SingleMuNonIsoTriggered,"l1SingleMuNonIsoTriggered/I");
  MuTrigMC->Branch("l2SingleMuNonIsoTriggered",&l2SingleMuNonIsoTriggered,"l2SingleMuNonIsoTriggered/I");
  MuTrigMC->Branch("l3SingleMuNonIsoTriggered",&l3SingleMuNonIsoTriggered,"l3SingleMuNonIsoTriggered/I");
  MuTrigMC->Branch("l1SingleMuIsoTriggered",&l1SingleMuIsoTriggered,"l1SingleMuIsoTriggered/I");
  MuTrigMC->Branch("l2SingleMuIsoPreTriggered",&l2SingleMuIsoPreTriggered,"l2SingleMuIsoPreTriggered/I");
  MuTrigMC->Branch("l2SingleMuIsoTriggered",&l2SingleMuIsoTriggered,"l2SingleMuIsoTriggered/I");
  MuTrigMC->Branch("l3SingleMuIsoPreTriggered",&l3SingleMuIsoPreTriggered,"l3SingleMuIsoPreTriggered/I");
  MuTrigMC->Branch("l3SingleMuIsoTriggered",&l3SingleMuIsoTriggered,"l3SingleMuIsoTriggered/I");
  MuTrigMC->Branch("l1DiMuNonIsoTriggered",&l1DiMuNonIsoTriggered,"l1DiMuNonIsoTriggered/I");
  MuTrigMC->Branch("l2DiMuNonIsoTriggered",&l2DiMuNonIsoTriggered,"l2DiMuNonIsoTriggered/I");
  MuTrigMC->Branch("l3DiMuNonIsoTriggered",&l3DiMuNonIsoTriggered,"l3DiMuNonIsoTriggered/I");
  MuTrigMC->Branch("l1DiMuIsoTriggered",&l1DiMuIsoTriggered,"l1DiMuIsoTriggered/I");
  MuTrigMC->Branch("l2DiMuIsoPreTriggered",&l2DiMuIsoPreTriggered,"l2DiMuIsoPreTriggered/I");
  MuTrigMC->Branch("l2DiMuIsoTriggered",&l2DiMuIsoTriggered,"l2DiMuIsoTriggered/I");
  MuTrigMC->Branch("l3DiMuIsoPreTriggered",&l3DiMuIsoPreTriggered,"l3DiMuIsoPreTriggered/I");
  MuTrigMC->Branch("l3DiMuIsoTriggered",&l3DiMuIsoTriggered,"l3DiMuIsoTriggered/I"); 
  MuTrigMC->Branch("triggerDecisions",&triggerDecisions);
  MuTrigMC->Branch("triggerNames",&triggerNames);
  // number of muons at each level
  MuTrigMC->Branch("nL1",&nL1,"nL1/I");
  MuTrigMC->Branch("nL2",&nL2,"nL2/I");
  MuTrigMC->Branch("nL3",&nL3,"nL3/I");
  MuTrigMC->Branch("nL3TracksFromL2",&nL3TracksFromL2,"nL3TracksFromL2/I");
  MuTrigMC->Branch("nL3Cands",&nL3Cands,"nL3Cands/I");
  MuTrigMC->Branch("nL3Seeds",&nL3Seeds,"nL3Seeds/I");
  // L3 muon information: the basics
  MuTrigMC->Branch("l3P",&l3P);
  MuTrigMC->Branch("l3Px",&l3Px);
  MuTrigMC->Branch("l3Py",&l3Py);
  MuTrigMC->Branch("l3Pz",&l3Pz);
  MuTrigMC->Branch("l3Pt",&l3Pt);
  MuTrigMC->Branch("l3PtError",&l3PtError);
  MuTrigMC->Branch("l3Pt90",&l3Pt90);
  MuTrigMC->Branch("l3Eta",&l3Eta);
  MuTrigMC->Branch("l3EtaError",&l3EtaError);
  MuTrigMC->Branch("l3Phi",&l3Phi);
  MuTrigMC->Branch("l3PhiError",&l3PhiError);
  MuTrigMC->Branch("l3D0",&l3D0);
  MuTrigMC->Branch("l3D0Error",&l3D0Error);
  MuTrigMC->Branch("l3NHits",&l3NHits);
  MuTrigMC->Branch("l3Charge",&l3Charge);
  MuTrigMC->Branch("l3Chi2",&l3Chi2);
  MuTrigMC->Branch("l3Ndof",&l3Ndof);
  // L3 Muon Isolation quantities
  MuTrigMC->Branch("l3CalIsoDeposit",&l3CalIsoDeposit);
  MuTrigMC->Branch("l3TrackIsoDeposit",&l3TrackIsoDeposit);
  MuTrigMC->Branch("l3IsoTrackDR",&l3IsoTrackDR);
  MuTrigMC->Branch("l3IsoTrackDRMinDelPt",&l3IsoTrackDRMinDelPt);
  // L3 Muon Track fitting parameters (with phi already declared above)
  MuTrigMC->Branch("l3Dsz",&l3Dsz);
  MuTrigMC->Branch("l3DszError",&l3DszError);
  MuTrigMC->Branch("l3Dxy",&l3Dxy);
  MuTrigMC->Branch("l3DxyError",&l3DxyError);
  MuTrigMC->Branch("l3Lambda",&l3Lambda);
  MuTrigMC->Branch("l3LambdaError",&l3LambdaError);
  MuTrigMC->Branch("l3Qoverp",&l3Qoverp);
  MuTrigMC->Branch("l3QoverpError",&l3QoverpError);
  MuTrigMC->Branch("l3ErrorMatrix",&l3ErrorMatrix);
  MuTrigMC->Branch("l3DetIds",&l3DetIds);
  MuTrigMC->Branch("l3SubdetIds",&l3SubdetIds);
  MuTrigMC->Branch("l3Component",&l3Component);
  MuTrigMC->Branch("l3NMuHits",&l3NMuHits);
  MuTrigMC->Branch("l3MuStationNumber",&l3MuStationNumber);
  MuTrigMC->Branch("l3RecHitsStatus",&l3RecHitsStatus);
  MuTrigMC->Branch("l3RecHitsX",&l3RecHitsX);
  MuTrigMC->Branch("l3RecHitsY",&l3RecHitsY);
  MuTrigMC->Branch("l3RecHitsZ",&l3RecHitsZ);
  // Indices for L3<->L2
  MuTrigMC->Branch("indexL2SeedingL3",&indexL2SeedingL3);
  MuTrigMC->Branch("indexL3SeededFromL2",&indexL3SeededFromL2);
  MuTrigMC->Branch("l2SeedsL3",&l2SeedsL3);

  // The muon error matrix
  MuTrigMC->Branch("muonErrorMatrix",&muonErrorMatrix);

  // information from hltL3TkTracksFromL2
  MuTrigMC->Branch("l3TrackP",&l3TrackP);
  MuTrigMC->Branch("l3TrackPx",&l3TrackPx);
  MuTrigMC->Branch("l3TrackPy",&l3TrackPy);
  MuTrigMC->Branch("l3TrackPz",&l3TrackPz);
  MuTrigMC->Branch("l3TrackPt",&l3TrackPt);
  MuTrigMC->Branch("l3TrackPtError",&l3TrackPtError);
  MuTrigMC->Branch("l3TrackEta",&l3TrackEta);
  MuTrigMC->Branch("l3TrackEtaError",&l3TrackEtaError);
  MuTrigMC->Branch("l3TrackPhi",&l3TrackPhi);
  MuTrigMC->Branch("l3TrackPhiError",&l3TrackPhiError);
  MuTrigMC->Branch("l3TrackD0",&l3TrackD0);
  MuTrigMC->Branch("l3TrackD0Error",&l3TrackD0Error);
  MuTrigMC->Branch("l3TrackNHits",&l3TrackNHits);
  MuTrigMC->Branch("l3TrackCharge",&l3TrackCharge);
  MuTrigMC->Branch("l3TrackChi2",&l3TrackChi2);
  MuTrigMC->Branch("l3TrackNdof",&l3TrackNdof);
  // L3TRACK Muon Track fitting parameters (with phi already declared above)
  MuTrigMC->Branch("l3TrackDsz",&l3TrackDsz);
  MuTrigMC->Branch("l3TrackDszError",&l3TrackDszError);
  MuTrigMC->Branch("l3TrackDxy",&l3TrackDxy);
  MuTrigMC->Branch("l3TrackDxyError",&l3TrackDxyError);
  MuTrigMC->Branch("l3TrackLambda",&l3TrackLambda);
  MuTrigMC->Branch("l3TrackLambdaError",&l3TrackLambdaError);
  MuTrigMC->Branch("l3TrackQoverp",&l3TrackQoverp);
  MuTrigMC->Branch("l3TrackQoverpError",&l3TrackQoverpError);
  MuTrigMC->Branch("l3TrackErrorMatrix",&l3TrackErrorMatrix);
  MuTrigMC->Branch("l3TrackDetIds",&l3TrackDetIds);
  MuTrigMC->Branch("l3TrackSubdetIds",&l3TrackSubdetIds);
  MuTrigMC->Branch("l3TrackRecHitsStatus",&l3TrackRecHitsStatus);
  MuTrigMC->Branch("l3TrackRecHitsX",&l3TrackRecHitsX);
  MuTrigMC->Branch("l3TrackRecHitsY",&l3TrackRecHitsY);
  MuTrigMC->Branch("l3TrackRecHitsZ",&l3TrackRecHitsZ);

  // L2 muon information: the basics
  MuTrigMC->Branch("l2P",&l2P);
  MuTrigMC->Branch("l2Px",&l2Px);
  MuTrigMC->Branch("l2Py",&l2Py);
  MuTrigMC->Branch("l2Pz",&l2Pz);
  MuTrigMC->Branch("l2Pt",&l2Pt);
  MuTrigMC->Branch("l2PtError",&l2PtError);
  MuTrigMC->Branch("l2Pt90",&l2Pt90);
  MuTrigMC->Branch("l2Eta",&l2Eta);
  MuTrigMC->Branch("l2EtaError",&l2EtaError);
  MuTrigMC->Branch("l2Phi",&l2Phi);
  MuTrigMC->Branch("l2PhiError",&l2PhiError);
  MuTrigMC->Branch("l2D0",&l2D0);
  MuTrigMC->Branch("l2D0Error",&l2D0Error);
  MuTrigMC->Branch("l2NHits",&l2NHits);
  MuTrigMC->Branch("l2Charge",&l2Charge);
  MuTrigMC->Branch("l2Chi2",&l2Chi2);
  MuTrigMC->Branch("l2Ndof",&l2Ndof);
  MuTrigMC->Branch("l2NSeeds",&l2NSeeds);
  MuTrigMC->Branch("l2DetIds",&l2DetIds);
  MuTrigMC->Branch("l2SubdetIds",&l2SubdetIds);
  MuTrigMC->Branch("l2Component",&l2Component);
  MuTrigMC->Branch("l2NMuHits",&l2NMuHits);
  MuTrigMC->Branch("l2MuStationNumber",&l2MuStationNumber);
  MuTrigMC->Branch("l2RecHitsStatus",&l2RecHitsStatus);
  MuTrigMC->Branch("l2RecHitsX",&l2RecHitsX);
  MuTrigMC->Branch("l2RecHitsY",&l2RecHitsY);
  MuTrigMC->Branch("l2RecHitsZ",&l2RecHitsZ);
  // L2 Muon Isolation quantities
  MuTrigMC->Branch("l2CalIsoDeposit",&l2CalIsoDeposit);
  // L2 Muon Track fitting parameters (with phi already declared above)
  MuTrigMC->Branch("l2Dsz",&l2Dsz);
  MuTrigMC->Branch("l2DszError",&l2DszError);
  MuTrigMC->Branch("l2Dxy",&l2Dxy);
  MuTrigMC->Branch("l2DxyError",&l2DxyError);
  MuTrigMC->Branch("l2Lambda",&l2Lambda);
  MuTrigMC->Branch("l2LambdaError",&l2LambdaError);
  MuTrigMC->Branch("l2Qoverp",&l2Qoverp);
  MuTrigMC->Branch("l2QoverpError",&l2QoverpError);
  MuTrigMC->Branch("l2ErrorMatrix",&l2ErrorMatrix);

  MuTrigMC->Branch("l1SeedsL2",&l1SeedsL2);
  MuTrigMC->Branch("indexL1SeedingL2",&indexL1SeedingL2);
  MuTrigMC->Branch("indexL2SeededFromL1",&indexL2SeededFromL1);

  MuTrigMC->Branch("l1P",&l1P);
  MuTrigMC->Branch("l1Pt",&l1Pt);
  MuTrigMC->Branch("l1Eta",&l1Eta);
  MuTrigMC->Branch("l1Phi",&l1Phi);
  MuTrigMC->Branch("l1Quality",&l1Quality);
  MuTrigMC->Branch("l1IsIso",&l1IsIso);
  MuTrigMC->Branch("l1IsMip",&l1IsMip);
  MuTrigMC->Branch("l1IsForward",&l1IsForward);
  MuTrigMC->Branch("l1IsRPC",&l1IsRPC);

  // Specific to MuTrigMC
  MuTrigMC->Branch("l3IsAssociated",&l3IsAssociated);
  MuTrigMC->Branch("l3ParentID",&l3ParentID);
  MuTrigMC->Branch("l3MotherBinNumber",&l3MotherBinNumber);
  MuTrigMC->Branch("l3AssociationVar",&l3AssociationVar);
  MuTrigMC->Branch("l3AssociatedSimMuonPt",&l3AssociatedSimMuonPt);
  MuTrigMC->Branch("l3AssociatedSimMuonEta",&l3AssociatedSimMuonEta);
  MuTrigMC->Branch("l3AssociatedSimMuonPhi",&l3AssociatedSimMuonPhi);
  MuTrigMC->Branch("l3AssociatedSimMuonNHits",&l3AssociatedSimMuonNHits);
  MuTrigMC->Branch("l3AssociatedSimMuonDetIds",&l3AssociatedSimMuonDetIds);
  MuTrigMC->Branch("l3AssociatedSimMuonNMuHits",&l3AssociatedSimMuonNMuHits);
  MuTrigMC->Branch("l3AssociatedSimMuonMuStationNumber",&l3AssociatedSimMuonMuStationNumber);

  MuTrigMC->Branch("l3AssociatedSimMuonDsz",&l3AssociatedSimMuonDsz);
  MuTrigMC->Branch("l3AssociatedSimMuonDxy",&l3AssociatedSimMuonDxy);
  MuTrigMC->Branch("l3AssociatedSimMuonLambda",&l3AssociatedSimMuonLambda);
  MuTrigMC->Branch("l3AssociatedSimMuonQoverP",&l3AssociatedSimMuonQoverP);

  MuTrigMC->Branch("l3TrackIsAssociated",&l3TrackIsAssociated);
  MuTrigMC->Branch("l3TrackParentID",&l3TrackParentID);
  MuTrigMC->Branch("l3TrackMotherBinNumber",&l3TrackMotherBinNumber);
  MuTrigMC->Branch("l3TrackAssociationVar",&l3TrackAssociationVar);
  MuTrigMC->Branch("l3TrackAssociatedSimMuonPt",&l3TrackAssociatedSimMuonPt);
  MuTrigMC->Branch("l3TrackAssociatedSimMuonEta",&l3TrackAssociatedSimMuonEta);
  MuTrigMC->Branch("l3TrackAssociatedSimMuonPhi",&l3TrackAssociatedSimMuonPhi);
  MuTrigMC->Branch("l3TrackAssociatedSimMuonNHits",&l3TrackAssociatedSimMuonNHits);
  MuTrigMC->Branch("l3TrackAssociatedSimMuonDetIds",&l3TrackAssociatedSimMuonDetIds);

  MuTrigMC->Branch("l3TrackAssociatedSimMuonDsz",&l3TrackAssociatedSimMuonDsz);
  MuTrigMC->Branch("l3TrackAssociatedSimMuonDxy",&l3TrackAssociatedSimMuonDxy);
  MuTrigMC->Branch("l3TrackAssociatedSimMuonLambda",&l3TrackAssociatedSimMuonLambda);
  MuTrigMC->Branch("l3TrackAssociatedSimMuonQoverP",&l3TrackAssociatedSimMuonQoverP);

  MuTrigMC->Branch("l2IsAssociated",&l2IsAssociated);
  MuTrigMC->Branch("l2ParentID",&l2ParentID);
  MuTrigMC->Branch("l2MotherBinNumber",&l2MotherBinNumber);
  MuTrigMC->Branch("l2AssociationVar",&l2AssociationVar);
  MuTrigMC->Branch("l2AssociatedSimMuonPt",&l2AssociatedSimMuonPt);
  MuTrigMC->Branch("l2AssociatedSimMuonEta",&l2AssociatedSimMuonEta);
  MuTrigMC->Branch("l2AssociatedSimMuonPhi",&l2AssociatedSimMuonPhi);
  MuTrigMC->Branch("l2AssociatedSimMuonNHits",&l2AssociatedSimMuonNHits);
  MuTrigMC->Branch("l2AssociatedSimMuonDetIds",&l2AssociatedSimMuonDetIds);
  MuTrigMC->Branch("l2AssociatedSimMuonNMuHits",&l2AssociatedSimMuonNMuHits);
  MuTrigMC->Branch("l2AssociatedSimMuonMuStationNumber",&l2AssociatedSimMuonMuStationNumber);

  MuTrigMC->Branch("l2AssociatedSimMuonDsz",&l2AssociatedSimMuonDsz);
  MuTrigMC->Branch("l2AssociatedSimMuonDxy",&l2AssociatedSimMuonDxy);
  MuTrigMC->Branch("l2AssociatedSimMuonLambda",&l2AssociatedSimMuonLambda);
  MuTrigMC->Branch("l2AssociatedSimMuonQoverP",&l2AssociatedSimMuonQoverP);

  MuTrigMC->Branch("nSimMuon",&nSimMuon,"nSimMuon/I");
  MuTrigMC->Branch("simMuonParentID",&simMuonParentID);
  MuTrigMC->Branch("simMuonMotherBinNumber",&simMuonMotherBinNumber);
  MuTrigMC->Branch("simMuonPt",&simMuonPt);
  MuTrigMC->Branch("simMuonEta",&simMuonEta);
  MuTrigMC->Branch("simMuonPhi",&simMuonPhi);
  MuTrigMC->Branch("simMuonNHits",&simMuonNHits);
  MuTrigMC->Branch("simMuonDetIds",&simMuonDetIds);
  MuTrigMC->Branch("simMuonNMuHits",&simMuonNMuHits);
  MuTrigMC->Branch("simMuonMuStationNumber",&simMuonMuStationNumber);

  MuTrigMC->Branch("simMuonDsz",&simMuonDsz);
  MuTrigMC->Branch("simMuonDxy",&simMuonDxy);
  MuTrigMC->Branch("simMuonLambda",&simMuonLambda);
  MuTrigMC->Branch("simMuonQoverP",&simMuonQoverP);

  MuTrigMC->Branch("simToL3Associated",&simToL3Associated);
  MuTrigMC->Branch("simToL3AssociationVar",&simToL3AssociationVar);
  MuTrigMC->Branch("simToL3RecoIndex",&simToL3RecoIndex);
  MuTrigMC->Branch("simToTkAssociated",&simToTkAssociated);
  MuTrigMC->Branch("simToTkAssociationVar",&simToTkAssociationVar);
  MuTrigMC->Branch("simToTkRecoIndex",&simToTkRecoIndex);
  MuTrigMC->Branch("simToL2Associated",&simToL2Associated);
  MuTrigMC->Branch("simToL2AssociationVar",&simToL2AssociationVar);
  MuTrigMC->Branch("simToL2RecoIndex",&simToL2RecoIndex);

  edm::LogInfo("IsoMuAnalyzer")<<"beginJob executed.  Problems not from here.";

}

// ------------ method called once each job just after ending the event loop  ------------
void 
IsoMuAnalyzer::endJob() {

  edm::LogInfo("IsoMuAnalyzer")<<"Starting to write the trees.  Changing to directory";
  theFile->cd();
  edm::LogInfo("IsoMuAnalyzer")<<"Starting to write the trees.  Writing MuTrigData";
  MuTrigData->Write();
  edm::LogInfo("IsoMuAnalyzer")<<"Starting to write the trees.  Writing MuTrigMC";
  MuTrigMC->Write();
  edm::LogInfo("IsoMuAnalyzer")<<"Finished writing trees.  Closing file.";
  theFile->Close();
  edm::LogInfo("IsoMuAnalyzer")<<"All done.  Nothing left to do.";

}

//define this as a plug-in
DEFINE_FWK_MODULE(IsoMuAnalyzer);
