#ifndef Workspace_JetMetDphiEventSelector_h_
#define Workspace_JetMetDphiEventSelector_h_


#include "Workspace/ConfigurableAnalysis/plugins/CosDphiEventSelector.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

typedef CosDphiEventSelector<pat::Jet,pat::MET> JetMetDphiEventSelector;

#endif
