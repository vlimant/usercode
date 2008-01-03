// -*- C++ -*-
//
// Package:    TTbarWithMuon
// Class:      TTbarWithMuon
// 
/**\class TTbarWithMuon TTbarWithMuon.cc RecoTracker/TTbarWithMuon/src/TTbarWithMuon.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jean-Roch Vlimant
//         Created:  Thu Nov  9 15:04:53 CST 2006
// $Id: TTbarWithMuon.h,v 1.2 2007/11/27 19:09:23 vlimant Exp $
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

#include <SimDataFormats/Track/interface/SimTrackContainer.h>
#include <SimDataFormats/Vertex/interface/SimVertexContainer.h>

#include <SimDataFormats/Track/interface/SimTrack.h>
#include <SimDataFormats/Vertex/interface/SimVertex.h>

#include <vector>

//
// class declaration
//

class TTbarWithMuon : public edm::EDFilter {
   public:
      explicit TTbarWithMuon(const edm::ParameterSet&);
      ~TTbarWithMuon();

   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------
  edm::InputTag HEPsourceLabel;
  std::vector<int> B_decay;
  bool W_decay;

  std::map<int,int> counts;
};

