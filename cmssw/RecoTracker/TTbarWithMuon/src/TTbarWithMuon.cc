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
// $Id$
//
//


#include <RecoTracker/TTbarWithMuon/interface/TTbarWithMuon.h>

#include <FWCore/MessageLogger/interface/MessageLogger.h>

#include <SimDataFormats/HepMCProduct/interface/HepMCProduct.h>
#include <iostream>


TTbarWithMuon::TTbarWithMuon(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  HEPsourceLabel = iConfig.getParameter<edm::InputTag>("HEPsourceLabel");
  B_decay = iConfig.getParameter< std::vector<int> >("B_decay");
  W_decay = iConfig.getParameter<bool>("W_decay");
}


TTbarWithMuon::~TTbarWithMuon()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
TTbarWithMuon::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;
   //purpose:
   //select ttbar event where a mu is coming 
   //   from W decay 
   //or 
   //   from B meson decay

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
	   if (!part->production_vertex()) continue;
	   if (part->production_vertex()->particles_in_size()==0) continue;
	   const HepMC::GenParticle * mother = *(part->production_vertex()->particles_in_const_begin());
	   if (mother){

	     if (abs(mother->pdg_id()) ==24 && W_decay)// from W decay 
	       { edm::LogInfo("TTbarWithMuon::filter(...)")<<" W decaying into mu";return true;}

	     if (B_decay.size()!=0){/* from b decay */
	       for (std::vector<int>::iterator B_pdg_it=B_decay.begin();B_pdg_it!=B_decay.end();++B_pdg_it)
		 if (abs(mother->pdg_id())==(*B_pdg_it))
		   { edm::LogInfo("TTbarWithMuon::filter(...)")<<" b-quark meson ("<<mother->pdg_id()<<") decaying into mu";return true;}}
	   }}
     }

   return false;
}

// ------------ method called once each job just before starting event loop  ------------
void 
TTbarWithMuon::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TTbarWithMuon::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(TTbarWithMuon);
