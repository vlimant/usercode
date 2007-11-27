#ifndef  MuonKinematics_H
#define  MuonKinematics_H

// -*- C++ -*-
//
// Package:    MuonKinematics
// Class:      MuonKinematics
// 
/**\class MuonKinematics MuonKinematics.cc RecoMuon/MuonKinematics/src/MuonKinematics.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jean-Roch Vlimant
//         Created:  Thu Nov 15 23:01:46 CET 2007
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <SimDataFormats/HepMCProduct/interface/HepMCProduct.h>

//
// class declaration
//

class MuonKinematics : public edm::EDFilter {
   public:
      explicit MuonKinematics(const edm::ParameterSet&);
      ~MuonKinematics();

   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------
  edm::InputTag HEPsourceLabel;
  double etaMin,etaMax;
};

#endif
