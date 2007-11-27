#ifndef RecoMuon_TrackerSeedGenerator_H
#define RecoMuon_TrackerSeedGenerator_H

/** \class TrackerSeedGenerator
 *  Generate seed from muon trajectory.
 */

#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
class Trajectory;
class TrackingRegion;
class MuonServiceProxy;

namespace edm {class ParameterSet; class Event; class EventSetup;}

class TrackerSeedGenerator {

public:
  typedef std::vector<TrajectorySeed> BTSeedCollection;  

  TrackerSeedGenerator() : theEvent(0), theProxyService(0) {}

  typedef std::pair<const Trajectory*, reco::TrackRef> TrackCand;

  virtual void init(const MuonServiceProxy *service);
  
  /// destructor
  virtual ~TrackerSeedGenerator() {}

  virtual void trackerSeeds(const TrackCand&, const TrackingRegion&, BTSeedCollection &);
    
  virtual void setEvent(const edm::Event&);

private:

  virtual void run(TrajectorySeedCollection &seeds, 
      const edm::Event &ev, const edm::EventSetup &es, const TrackingRegion& region) {} 

  const edm::Event * theEvent;
  const MuonServiceProxy * theProxyService;
  
};

#endif

