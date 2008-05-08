#include "Workspace/ConfigurableAnalysis/plugins/StringCutEventSelector.h"
#include "Workspace/ConfigurableAnalysis/plugins/StringCutsEventSelector.h"
#include "DataFormats/PatCandidates/interface/MET.h"

typedef StringCutEventSelector<pat::MET> aMETEventSelector;
typedef StringCutsEventSelector<pat::MET> METSEventSelector;
