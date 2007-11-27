#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/eventsetupdata_registration_macro.h"

#include "RecoTracker/TrackProducerFromSimTrack/interface/TrackProducerFromSimTrack.h"

//define this as a plug-in
DEFINE_FWK_MODULE(TrackProducerFromSimTrack);
