#ifndef ConfigurableAnalysis_CachingVariableFactory_H
#define ConfigurableAnalysis_CachingVariableFactory_H

#include "FWCore/PluginManager/interface/PluginFactory.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Workspace/ConfigurableAnalysis/interface/CachingVariable.h"

typedef edmplugin::PluginFactory< CachingVariable* (std::string , edm::ParameterSet&) > CachingVariableFactory;

#endif
