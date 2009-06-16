#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include <SimDataFormats/GeneratorProducts/interface/HepMCProduct.h>
#include "HepMC/GenEvent.h"
#include "TLorentzVector.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <DataFormats/Candidate/interface/Particle.h>

class MotherSearch {
public:
  MotherSearch( const SimTrack * isimtk,
		edm::Handle<edm::SimTrackContainer> & SimTk,
		edm::Handle<edm::SimVertexContainer> & SimVtx,
		edm::Handle<edm::HepMCProduct> & hepmc);

  bool IsValid() { return SimIsValid() || GenIsValid();}
  bool SimIsValid() { return (!useGen && Sim_vertex && Sim_mother && simtrack);}
  bool GenIsValid() { return (useGen && Gen_vertex && Gen_mother && simtrack && gentrack); }

  reco::Particle particle();
  reco::Particle mother();

  const SimTrack * simtrack;
  const SimVertex * Sim_vertex;
  const SimTrack * Sim_mother;
  bool useGen;
  const HepMC::GenParticle * gentrack;
  const HepMC::GenVertex * Gen_vertex;
  const HepMC::GenParticle *Gen_mother;
};
