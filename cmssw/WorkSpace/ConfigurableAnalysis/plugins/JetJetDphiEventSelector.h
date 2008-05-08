#ifndef Workspace_JetJetDphiEventSelector_h_
#define Workspace_JetJetDphiEventSelector_h_


#include "Workspace/ConfigurableAnalysis/plugins/CosDphiEventSelector.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

typedef CosDphiEventSelector<pat::Jet,pat::Jet> JetJetDphiEventSelector;

#endif
