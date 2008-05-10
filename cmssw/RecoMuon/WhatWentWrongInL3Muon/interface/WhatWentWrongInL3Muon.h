#ifndef WhatWentWrongInL3Muon_H
#define WhatWentWrongInL3Muon_H

// -*- C++ -*-
//
// Package:    WhatWentWrongInL3Muon
// Class:      WhatWentWrongInL3Muon
// 
/**\class WhatWentWrongInL3Muon WhatWentWrongInL3Muon.cc RecoMuon/WhatWentWrongInL3Muon/src/WhatWentWrongInL3Muon.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jean-Roch Vlimant
//         Created:  Wed Sep 12 02:24:59 CEST 2007
// $Id: WhatWentWrongInL3Muon.h,v 1.1 2007/11/27 20:36:32 vlimant Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
//
// class decleration
//

#include "DataFormats/MuonReco/interface/MuonTrackLinks.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include <TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h>
#include <TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h>


#include <MagneticField/Records/interface/IdealMagneticFieldRecord.h>
#include <MagneticField/Engine/interface/MagneticField.h>


#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <TH2F.h>
#include <TFile.h>
#include <TMath.h>

#include <sstream>

#include "DQMServices/Core/interface/DaqMonitorBEInterface.h"
#include "DQMServices/Daemon/interface/MonitorDaemon.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

#include <Geometry/Records/interface/TrackerDigiGeometryRecord.h>
#include <Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h>
#include <Geometry/CommonDetUnit/interface/GeomDet.h>

#include <TrackingTools/Records/interface/TrackingComponentsRecord.h>
#include <TrackingTools/GeomPropagators/interface/Propagator.h>

#include <TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h>

#include <DataFormats/GeometrySurface/interface/BoundPlane.h>

class nodePlotter{
public:
  //constructor
  nodePlotter(std::string s):nodeName(s),owner(0){}
  nodePlotter():nodeName("unknown"),owner(0){}
  void setName(std::string s){nodeName=s;}
  void setOwner(nodePlotter * np){owner=np;}
  
  nodePlotter &  registerNode(std::string s) {
    nodePlotter  & inserted =dirMap[s];
    inserted.nodeName=s;
    inserted.owner=this;
    return inserted;}
  
  nodePlotter &  registerNode(edm::InputTag & tag){ return registerNode(encode(tag));}
  
  MonitorElement* registerElement(std::string s, MonitorElement* e) {elementMap[s]=e; return e;}
  MonitorElement* registerElement(MonitorElement* e) {elementMap[dynamic_cast<MonitorElementRootObject*>(e)->getName()]=e; return e;}
  //accessor
  nodePlotter & dir(std::string s) {
    nodePlotter & n=dirMap[s];
      if (n.name()=="unknown"){edm::LogError("nodePlotter")<<s<<" node does not exist in: "<<fullName();}
      return n; }

  nodePlotter & dir(edm::InputTag & tag) {return dir(encode(tag)); }

  MonitorElement * element(std::string s) {
    LogDebug("nodePlotter")<<"accessing: "<<s<<" in:"<<fullName();
    MonitorElement * e = elementMap[s];
    if (!e) {edm::LogError("nodePlotter")<<s<<" element does not exist in: "<<fullName();}
    return e;}
  const std::string & name(){return nodeName;}
  const std::string fullName() {
    if (owner) return (owner->fullName()+"/"+nodeName);
    else return nodeName; }
  std::string encode(edm::InputTag & tag, std::string separator="_"){
    return tag.label()
	+(tag.instance().empty()?std::string() : separator+tag.instance())
      +(tag.process().empty()?std::string() : separator+tag.process());}
  
private:
  std::map<std::string, nodePlotter> dirMap;
  std::string nodeName;
  std::map<std::string, MonitorElement*> elementMap;
  nodePlotter * owner;
};


class WhatWentWrongInL3Muon : public edm::EDAnalyzer {
   public:
      explicit WhatWentWrongInL3Muon(const edm::ParameterSet&);
      ~WhatWentWrongInL3Muon();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      void oneLocal(nodePlotter & node1, reco::RecoToSimCollection & recSimColl, edm::RefToBase<reco::Track> & refL2, FreeTrajectoryState & l2State);
      void one(nodePlotter & node1, edm::RefToBase<reco::Track> & refL2, FreeTrajectoryState & l2State);
      // ----------member data ---------------------------
  std::string category;
  
  DaqMonitorBEInterface* theDQM;
  
  edm::InputTag L2Label;
  std::vector<edm::InputTag> L3Labels;

  std::string theRootFileName;

  nodePlotter plotter;

  edm::InputTag trackingParticleLabel;
  std::string theAssocName;
  edm::ESHandle<TrackAssociatorBase> theAssociator;
  edm::ESHandle<TrackerGeometry> theTrackerGeometry;
  std::string thePropagatorName;
  edm::ESHandle<Propagator> thePropagator;

  bool theTellMe;
};

#endif
