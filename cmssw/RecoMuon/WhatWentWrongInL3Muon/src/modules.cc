#include "FWCore/Framework/interface/MakerMacros.h"
//define this as a plug-in

#include "RecoMuon/WhatWentWrongInL3Muon/interface/WhatWentWrongInL3Muon.h"
#include "RecoMuon/WhatWentWrongInL3Muon/interface/MuonKinematics.h"

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(WhatWentWrongInL3Muon);
DEFINE_ANOTHER_FWK_MODULE(MuonKinematics);
