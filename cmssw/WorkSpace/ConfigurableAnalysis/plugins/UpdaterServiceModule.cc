#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/ServiceRegistry/interface/ServiceMaker.h"

#include "Workspace/ConfigurableAnalysis/interface/UpdaterService.h"
DEFINE_FWK_SERVICE( UpdaterService );
//DEFINE_FWK_SERVICE( VariableHelperService );
