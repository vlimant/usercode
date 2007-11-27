#include "FWCore/Framework/interface/MakerMacros.h"

#include "RecoMuon/TrackerSeedGenerator/interface/TrackerSeedGeneratorFactory.h"

#include "RecoMuon/TrackerSeedGenerator/interface/TrackerSeedGeneratorBC.h"
#include "RecoMuon/TrackerSeedGenerator/interface/TSGFromOrderedHits.h"
#include "RecoMuon/TrackerSeedGenerator/interface/TSGForRoadSearch.h"
#include "RecoMuon/TrackerSeedGenerator/interface/TSGFromPropagation.h"

#include "RecoMuon/TrackerSeedGenerator/interface/DualByHitFractionTSG.h"
#include "RecoMuon/TrackerSeedGenerator/interface/DualByEtaTSG.h"
#include "RecoMuon/TrackerSeedGenerator/interface/DualByZTSG.h"
#include "RecoMuon/TrackerSeedGenerator/interface/CombinedTSG.h"

DEFINE_EDM_PLUGIN(TrackerSeedGeneratorFactory, TrackerSeedGeneratorBC, "TrackerSeedGeneratorBC");
DEFINE_EDM_PLUGIN(TrackerSeedGeneratorFactory, TSGFromOrderedHits, "TSGFromOrderedHits");
DEFINE_EDM_PLUGIN(TrackerSeedGeneratorFactory, TSGForRoadSearch, "TSGForRoadSearch");
DEFINE_EDM_PLUGIN(TrackerSeedGeneratorFactory, TSGFromPropagation, "TSGFromPropagation");
DEFINE_EDM_PLUGIN(TrackerSeedGeneratorFactory, DualByHitFractionTSG, "DualByHitFractionTSG");
DEFINE_EDM_PLUGIN(TrackerSeedGeneratorFactory, DualByEtaTSG, "DualByEtaTSG");
DEFINE_EDM_PLUGIN(TrackerSeedGeneratorFactory, DualByZTSG, "DualByZTSG");
DEFINE_EDM_PLUGIN(TrackerSeedGeneratorFactory, CombinedTSG, "CombinedTSG");
