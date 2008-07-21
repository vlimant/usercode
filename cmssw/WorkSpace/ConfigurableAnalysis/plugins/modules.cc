#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/ServiceRegistry/interface/ServiceMaker.h"

DEFINE_SEAL_MODULE();


#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "PhysicsTools/UtilAlgos/interface/ObjectSelector.h"
#include "PhysicsTools/UtilAlgos/interface/SortCollectionSelector.h"
#include "PhysicsTools/Utilities/interface/PtComparator.h"

typedef ObjectSelector<SortCollectionSelector<reco::GenParticleCollection,GreaterByPt<reco::GenParticle> > > GenParticleSorterByPt;

DEFINE_FWK_MODULE(GenParticleSorterByPt);

#include "Workspace/ConfigurableAnalysis/interface/ProcessIdSplitter.h"
DEFINE_EDM_PLUGIN(CachingVariableFactory, ProcessIdSplitter, "ProcessIdSplitter");


#include "Workspace/ConfigurableAnalysis/interface/StringBasedNTupler.h"
DEFINE_EDM_PLUGIN(NTuplerFactory, StringBasedNTupler, "StringBasedNTupler");
