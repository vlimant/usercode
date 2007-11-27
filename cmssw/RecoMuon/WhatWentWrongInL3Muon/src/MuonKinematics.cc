#include "RecoMuon/WhatWentWrongInL3Muon/interface/MuonKinematics.h"

MuonKinematics::MuonKinematics(const edm::ParameterSet& iConfig)
{
  HEPsourceLabel = iConfig.getParameter<edm::InputTag>("HEPsourceLabel");
  etaMin = iConfig.getParameter<double>("etaMin");
  etaMax = iConfig.getParameter<double>("etaMax");
}


MuonKinematics::~MuonKinematics()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
MuonKinematics::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   Handle<HepMCProduct> hepmc;
   iEvent.getByLabel(HEPsourceLabel, hepmc);
   if(!hepmc.isValid())
     {return false;}
   
   const HepMC::GenEvent* evt = hepmc->GetEvent();
   
   for(HepMC::GenEvent::particle_const_iterator pitr = evt->particles_begin(); pitr !=evt->particles_end();++pitr)
     {
       const HepMC::GenParticle * part = (*pitr);
       
       int ipdg = part->pdg_id();
       
       if (abs(ipdg) == 13) //mu
	 {
	   double feta =fabs(part->momentum().eta());
	   if (feta> etaMin &&  feta<etaMax) return true;
	 }
     }
   return false;
}

// ------------ method called once each job just before starting event loop  ------------
void 
MuonKinematics::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonKinematics::endJob() {
}
