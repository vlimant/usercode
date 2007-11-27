#include "RecoMuon/TrackerSeedGenerator/interface/CompositeTSG.h"

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "RecoMuon/TrackerSeedGenerator/interface/TrackerSeedGeneratorFactory.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

CompositeTSG::CompositeTSG(const edm::ParameterSet & par){
  theCategory = "CompositeTSG";

  std::string passingOnTo = par.getParameter<std::string>("StateOnTrackerBoundOutPropagator");
  /*
    std::vector<edm::ParameterSet> TSGsPSet = par.getParameter<std::vector< edm::ParameterSet> > ("TSGsPSet");
    
    for (std::vector<edm::ParameterSet>::iterator PSit = TSGsPSet.begin(); PSit!= TSGsPSet.end(); ++PSit){
    //    edm::ParameterSet TSGpset = PSit->getParameter<edm::ParameterSet>("SeedGeneratorParameters");
    edm::ParameterSet & TSGpset = *PSit;
    TSGpset.addParameter<std::string>("StateOnTrackerBoundOutPropagator",passingOnTo);
    std::string SeedGenName = TSGpset.getParameter<std::string>("ComponentName");
    theNames.push_back(SeedGenName);
    theTSGs.push_back(TrackerSeedGeneratorFactory::get()->create(SeedGenName,TSGpset));
    }
  */

  //alternative way of doing things
  std::vector<std::string> PSetNames =  par.getParameter<std::vector<std::string> >("PSetNames");
  for (std::vector<std::string>::iterator nIt = PSetNames.begin();nIt!=PSetNames.end();nIt++){
    edm::ParameterSet TSGpset = par.getParameter<edm::ParameterSet>(*nIt);
    TSGpset.addParameter<std::string>("StateOnTrackerBoundOutPropagator",passingOnTo);
    std::string SeedGenName = TSGpset.getParameter<std::string>("ComponentName");
    theNames.push_back((*nIt)+":"+SeedGenName);
    theTSGs.push_back(TrackerSeedGeneratorFactory::get()->create(SeedGenName,TSGpset));
  }

}

CompositeTSG::~CompositeTSG(){
  //delete the components ?
}


void CompositeTSG::init(const MuonServiceProxy* service){
  theProxyService = service;
  for (uint iTSG=0; iTSG!=theTSGs.size();iTSG++){
    theTSGs[iTSG]->init(service);}
}

void CompositeTSG::setEvent(const edm::Event &event){
  for (uint iTSG=0; iTSG!=theTSGs.size();iTSG++){
    theTSGs[iTSG]->setEvent(event);}
}


/*
  void CompositeTSG::trackerSeeds(const TrackCand & muonTrackCand, const TrackingRegion& region, std::vector<TrajectorySeed> & result){
  edm::LogError(theCategory)<<"I should never been called. throw exception.";
  return ;
}
*/
 
