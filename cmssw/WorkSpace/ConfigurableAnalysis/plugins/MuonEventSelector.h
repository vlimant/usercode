#include "Workspace/ConfigurableAnalysis/plugins/StringCutEventSelector.h"
#include "Workspace/ConfigurableAnalysis/plugins/StringCutsEventSelector.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

typedef StringCutEventSelector<pat::Muon> MuonEventSelector;
typedef StringCutsEventSelector<pat::Muon> MuonSEventSelector;
