#include "Workspace/ConfigurableAnalysis/plugins/StringCutEventSelector.h"
#include "Workspace/ConfigurableAnalysis/plugins/StringCutsEventSelector.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

typedef StringCutEventSelector<pat::Jet> aJetEventSelector;
typedef StringCutsEventSelector<pat::Jet> JetSEventSelector;
