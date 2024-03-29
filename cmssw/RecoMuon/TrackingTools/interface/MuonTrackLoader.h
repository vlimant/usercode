#ifndef RecoMuon_TrackingTools_MuonTrackLoader_H
#define RecoMuon_TrackingTools_MuonTrackLoader_H

/** \class MuonTrackLoader
 *  Class to load the tracks in the event, it provide some common functionalities
 *  both for all the RecoMuon producers.
 *
 *  $Date: 2007/07/13 20:38:03 $
 *  $Revision: 1.20 $
 *  \author R. Bellan - INFN Torino <riccardo.bellan@cern.ch>
 */

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/Common/interface/OrphanHandle.h"
#include "RecoMuon/TrackingTools/interface/MuonCandidate.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

namespace edm {class Event; class EventSetup; class ParameterSet;}

class Trajectory;
class Propagator;
class MuonServiceProxy;
class MuonUpdatorAtVertex;
class TrajectorySmoother;

class MuonTrackLoader {
  public:

    typedef MuonCandidate::TrajectoryContainer TrajectoryContainer;
    typedef MuonCandidate::CandidateContainer CandidateContainer;

    /// Constructor for the STA reco the args must be specify!
    MuonTrackLoader(edm::ParameterSet &parameterSet, const MuonServiceProxy *service =0);

    /// Destructor
    virtual ~MuonTrackLoader();
   
    /// Convert the trajectories into tracks and load the tracks in the event
/*     edm::OrphanHandle<reco::TrackCollection> loadTracks(const TrajectoryContainer&,  */
/*                                                         edm::Event&); */

    edm::OrphanHandle<reco::TrackCollection> loadTracks(const TrajectoryContainer&, 
                                                        edm::Event&,const std::string& = "", 
							bool = true);

    /// Convert the trajectories into tracks and load the tracks in the event
    edm::OrphanHandle<reco::MuonTrackLinksCollection> loadTracks(const CandidateContainer&,
								 edm::Event&); 
  
  private:
 
    /// Build a track at the PCA WITHOUT any vertex constriant
    std::pair<bool,reco::Track> buildTrackAtPCA(const Trajectory& trajectory) const;

    /// Takes a track at the PCA and applies the vertex constriant
    reco::Track buildTrackUpdatedAtPCA(const reco::Track& trackAtPCA) const;

    reco::TrackExtra buildTrackExtra(const Trajectory&) const;

    double computeNDOF(const Trajectory& trajectory) const;
  
    const MuonServiceProxy *theService;

    bool theUpdatingAtVtx;
    MuonUpdatorAtVertex *theUpdatorAtVtx;

    bool theTrajectoryFlag;

    bool theSmoothingStep;
    std::string theSmootherName;
    edm::ESHandle<TrajectorySmoother> theSmoother;

    /// Label for L2SeededTracks
    std::string theL2SeededTkLabel; 
    bool thePutTkTrackFlag;
    bool theSmoothTkTrackFlag;
    bool theAllowNoVtxFlag;
};
#endif
