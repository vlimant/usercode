#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/ServiceRegistry/interface/ServiceMaker.h"

DEFINE_SEAL_MODULE();

#include "Workspace/EventSelectors/interface/EventSelectorFactory.h"

#include "Workspace/ConfigurableAnalysis/plugins/JetMetDphiEventSelector.h"
#include "Workspace/ConfigurableAnalysis/plugins/JetJetDphiEventSelector.h"
#include "Workspace/ConfigurableAnalysis/plugins/JetEventSelector.h"
#include "Workspace/ConfigurableAnalysis/plugins/METEventSelector.h"
#include "Workspace/ConfigurableAnalysis/plugins/MuonEventSelector.h"
#include "Workspace/ConfigurableAnalysis/plugins/VariableEventSelector.h"

DEFINE_EDM_PLUGIN(EventSelectorFactory, JetMetDphiEventSelector, "JetMetDphiEventSelector");
DEFINE_EDM_PLUGIN(EventSelectorFactory, JetJetDphiEventSelector, "JetJetDphiEventSelector");
DEFINE_EDM_PLUGIN(EventSelectorFactory, aJetEventSelector, "aJetEventSelector");
DEFINE_EDM_PLUGIN(EventSelectorFactory, JetSEventSelector, "JetSEventSelector");
DEFINE_EDM_PLUGIN(EventSelectorFactory, aMETEventSelector, "aMETEventSelector");
DEFINE_EDM_PLUGIN(EventSelectorFactory, METSEventSelector, "METSEventSelector");
DEFINE_EDM_PLUGIN(EventSelectorFactory, MuonEventSelector, "MuonEventSelector");
DEFINE_EDM_PLUGIN(EventSelectorFactory, MuonSEventSelector, "MuonSEventSelector");
DEFINE_EDM_PLUGIN(EventSelectorFactory, VariableEventSelector, "VariableEventSelector");

#include "Workspace/ConfigurableAnalysis/interface/CachingVariableFactory.h"
#include "Workspace/ConfigurableAnalysis/interface/CachingVariable.h"

namespace configurableAnalysis{
  char Jet[]="pat::Jet";
  char Muon[]="pat::Muon";
  char MET[]="pat::MET";
}
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
typedef ExpressionVariable<pat::Jet,configurableAnalysis::Jet> JetExpressionVariable;
typedef ExpressionVariable<pat::MET,configurableAnalysis::MET> METExpressionVariable;
typedef ExpressionVariable<pat::Muon,configurableAnalysis::Muon> MuonExpressionVariable;

DEFINE_EDM_PLUGIN(CachingVariableFactory, JetExpressionVariable, "JetExpressionVariable");
DEFINE_EDM_PLUGIN(CachingVariableFactory, METExpressionVariable, "METExpressionVariable");
DEFINE_EDM_PLUGIN(CachingVariableFactory, MuonExpressionVariable, "MuonExpressionVariable");

typedef CosDphiVariable<pat::Jet,configurableAnalysis::Jet,pat::Muon,configurableAnalysis::Muon> JetMuonCosDphiVariable;
typedef CosDphiVariable<pat::Jet,configurableAnalysis::Jet,pat::MET,configurableAnalysis::MET> JetMETCosDphiVariable;
typedef CosDphiVariable<pat::Jet,configurableAnalysis::Jet,pat::Jet,configurableAnalysis::Jet> JetJetCosDphiVariable;

DEFINE_EDM_PLUGIN(CachingVariableFactory, JetMuonCosDphiVariable, "JetMuonCosDphiVariable");
DEFINE_EDM_PLUGIN(CachingVariableFactory, JetMETCosDphiVariable, "JetMETCosDphiVariable");
DEFINE_EDM_PLUGIN(CachingVariableFactory, JetJetCosDphiVariable, "JetJetCosDphiVariable");


DEFINE_EDM_PLUGIN(CachingVariableFactory, Power, "Power");
DEFINE_EDM_PLUGIN(CachingVariableFactory, CalculateMHT, "CalculateMHT");
DEFINE_EDM_PLUGIN(CachingVariableFactory, VarSplitter, "VarSplitter");
DEFINE_EDM_PLUGIN(CachingVariableFactory, ProcessIdSplitter, "ProcessIdSplitter");

