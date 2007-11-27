#ifndef RecoMuon_TrackingTools_MuonTrajectoryCleaner_H
#define RecoMuon_TrackingTools_MuonTrajectoryCleaner_H

/** \class MuonTrajectoryCleaner
 *  No description available.
 *
 *  $Date: 2006/11/23 02:29:36 $
 *  $Revision: 1.6 $
 *  \author R. Bellan - INFN Torino
 */

#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "RecoMuon/TrackingTools/interface/MuonCandidate.h"
#include <vector>

class MuonTrajectoryCleaner {
 public:
  typedef MuonCandidate::TrajectoryContainer TrajectoryContainer;
  typedef MuonCandidate::CandidateContainer CandidateContainer;
  

  /// Constructor
  MuonTrajectoryCleaner(){};

  /// Destructor
  virtual ~MuonTrajectoryCleaner(){};

  // Operations

  /// Clean the trajectories container, erasing the (worst) clone trajectory
  void clean(TrajectoryContainer &muonTrajectories); //used by reference...

  /// Clean the candidates container, erasing the (worst) clone trajectory
  void clean(CandidateContainer &muonTrajectories); //used by reference...

protected:

private:

};
#endif

