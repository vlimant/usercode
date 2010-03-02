#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TFile.h"
//#include "Math/LorentzVector.h"
#include <vector>
#include <map>
#include <utility>
#include "CMS1.h"
#include "Style.C"

setTDRStyle();

//int ScanTree ( TTree* tree) {
int ScanTree ( TTree* tree, char *fileName, int type) {

  switch (type) {
  case 1:
    cout << "This is ppMuX" << endl;
    break;
  case 2:
    cout << "This is ttbar" << endl;
    break;
  case 3:
    cout << "This is ZMM" << endl;
    break;
  default:
    cout << "Not sure what this is supposed to be." << endl;
  }

  Init(tree);
  TFile *histFile = new TFile(fileName,"RECREATE");

  // Sanity check plots for L1 (Including figures for MuonHLT Document)
  TH1F *level1Number = new TH1F("level1Number","level1Number",7,0,7);
  TH1F *level1Pt = new TH1F("level1Pt","level1Pt",150,0,150);
  TH1F *level1Eta = new TH1F("level1Eta","level1Eta",100,-2.5,2.5);
  TH1F *level1Phi = new TH1F("level1Phi","level1Phi",160,-3.2,3.2);
  TH1F *level1Quality = new TH1F("level1Quality","level1Quality",8,0,8);
  TH1F *level1Quality_afterL1 = new TH1F("level1Quality_afterL1","level1Quality_afterL1",8,0,8);
  TH1F *level1SeedsL2 = new TH1F("level1SeedsL2","level1SeedsL2",2,0,2); 
  // Then required maps for l1 section
  TH2F *level1GenPt = new TH2F("level1GenPt","level1GenPt",100,0,500,150,0,150); // Binned for single mu 0 - 500
  TH2F *level1GenEta = new TH2F("level1GenEta","level1GenEta",100,-2.5,2.5,100,-2.5,2.5);
  TH2F *level1GenPhi = new TH2F("level1GenPhi","level1GenPhi",160,-3.2,3.2,160,-3.2,3.2);
  TH2F *level1EtaQuality = new TH2F("level1EtaQuality","level1EtaQuality",100,-2.5,2.5,8,0,8);
  TH2F *level1PhiQuality = new TH2F("level1PhiQuality","level1PhiQuality",160,-3.2,3.2,8,0,8);
  // Put in the L1->L2 diagnositcs here.
  // Failing to seed L2
  TH1F *level1Pt_FailingL2Seeding =  new TH1F("level1Pt_FailingL2Seeding","level1Pt_FailingL2Seeding",150,0,150);
  TH1F *level1Eta_FailingL2Seeding =  new TH1F("level1Eta_FailingL2Seeding","level1Eta_FailingL2Seeding",100,-2.5,2.5);
  TH1F *level1Phi_FailingL2Seeding =  new TH1F("level1Phi_FailingL2Seeding","level1Phi_FailingL2Seeding",160,-3.2,3.2);
  TH1F *level1Quality_FailingL2Seeding =  new TH1F("level1Quality_FailingL2Seeding","level1Quality_FailingL2Seeding",8,0,8);
  // 2-D map
  TH2F *level1PtQuality_FailingL2Seeding = new TH2F("level1PtQuality_FailingL2Seeding","level1PtQuality_FailingL2Seeding",150,0,150,8,0,8);
  TH2F *level1EtaQuality_FailingL2Seeding = new TH2F("level1EtaQuality_FailingL2Seeding","level1EtaQuality_FailingL2Seeding",100,-2.5,2.5,8,0,8);
  TH2F *level1PhiQuality_FailingL2Seeding = new TH2F("level1PhiQuality_FailingL2Seeding","level1PhiQuality_FailingL2Seeding",160,-3.2,3.2,8,0,8);
  TH2F *level1EtaPhi_FailingL2Seeding = new TH2F("level1EtaPhi_FailingL2Seeding","level1EtaPhi_FailingL2Seeding",100,-2.5,2.5,160,-3.2,3.2);
  // Seeding L2
  TH2F *level1PtQuality_PassingL2Seeding = new TH2F("level1PtQuality_PassingL2Seeding","level1PtQuality_PassingL2Seeding",150,0,150,8,0,8);
  TH2F *level1EtaQuality_PassingL2Seeding = new TH2F("level1EtaQuality_PassingL2Seeding","level1EtaQuality_PassingL2Seeding",100,-2.5,2.5,8,0,8);
  TH2F *level1PhiQuality_PassingL2Seeding = new TH2F("level1PhiQuality_PassingL2Seeding","level1PhiQuality_PassingL2Seeding",160,-3.2,3.2,8,0,8);
  TH2F *level1EtaPhi_PassingL2Seeding = new TH2F("level1EtaPhi_PassingL2Seeding","level1EtaPhi_PassingL2Seeding",100,-2.5,2.5,160,-3.2,3.2);

  // Stations before/after L1 Filter
  TH1F *simStationsBeforeL1 = new TH1F("simStationsBeforeL1","simStationsBeforeL1",5,0,5);
  TH1F *simStationsAfterL1 = new TH1F("simStationsAfterL1","simStationsAfterL1",5,0,5);

  // Sanity check plots for L2 (Including figures for MuonHLT Document)
  // Reco level only
  TH1F *level2Number = new TH1F("level2Number","level2Number",5,0,5);
  TH1F *level2Pt = new TH1F("level2Pt","level2Pt",50,0,100);
  TH1F *level2Eta = new TH1F("level2Eta","level2Eta",90,-3,3);
  TH1F *level2Phi = new TH1F("level2Phi","level2Phi",96,-3.2,3.2);
  TH1F *level2D0 = new TH1F("level2D0","level2D0",80,-0.2,0.2);
  TH1F *level2NHits = new TH1F("level2NHits","level2NHits",70,0,70);
  TH1F *level2CalDeposit_hist = new TH1F("level2CalDeposit_hist","CAL deposit around L2 #mu",100,0,10);
  TH1F *level2CalDeposit_dist_cd = new TH1F("level2CalDeposit_dist_cd","CAL deposit around L2 #mu",100,0,10);
  // Here we want l2 Cal Deposit by parent ID
  TH1F *level2CalDeposit_1 = new TH1F("level2CalDeposit_1","#pi^{+/-}",100,0,10);
  TH1F *level2CalDeposit_2 = new TH1F("level2CalDeposit_2","K",100,0,10);
  TH1F *level2CalDeposit_3 = new TH1F("level2CalDeposit_3","D",100,0,10);
  TH1F *level2CalDeposit_4 = new TH1F("level2CalDeposit_4","B",100,0,10);
  TH1F *level2CalDeposit_5 = new TH1F("level2CalDeposit_5","#Lambda_{b}",100,0,10);
  TH1F *level2CalDeposit_6 = new TH1F("level2CalDeposit_6","J/#Psi",100,0,10);
  TH1F *level2CalDeposit_7 = new TH1F("level2CalDeposit_7","#Upsilon (ns)",100,0,10);
  TH1F *level2CalDeposit_8 = new TH1F("level2CalDeposit_8","#mu/W^{+/-}",100,0,10);
  TH1F *level2CalDeposit_9 = new TH1F("level2CalDeposit_9","#mu/Z",100,0,10);
  TH1F *level2CalDeposit_10 = new TH1F("level2CalDeposit_10","#tau^{+/-}",100,0,10);
  TH1F *level2CalDeposit_11 = new TH1F("level2CalDeposit_11","#mu^{+/-}",100,0,10);
  TH1F *level2CalDeposit_12 = new TH1F("level2CalDeposit_12","other",100,0,10);
  TH1F *level2CalDeposit_13 = new TH1F("level2CalDeposit_13","non-assoc",100,0,10);

  TH1F *level2CalDeposit_1_cd = new TH1F("level2CalDeposit_1_cd","#pi^{+/-}",100,0,10);
  TH1F *level2CalDeposit_2_cd = new TH1F("level2CalDeposit_2_cd","K",100,0,10);
  TH1F *level2CalDeposit_3_cd = new TH1F("level2CalDeposit_3_cd","D",100,0,10);
  TH1F *level2CalDeposit_4_cd = new TH1F("level2CalDeposit_4_cd","B",100,0,10);
  TH1F *level2CalDeposit_5_cd = new TH1F("level2CalDeposit_5_cd","#Lambda_{b}",100,0,10);
  TH1F *level2CalDeposit_6_cd = new TH1F("level2CalDeposit_6_cd","J/#Psi",100,0,10);
  TH1F *level2CalDeposit_7_cd = new TH1F("level2CalDeposit_7_cd","#Uptilon (ns)",100,0,10);
  TH1F *level2CalDeposit_8_cd = new TH1F("level2CalDeposit_8_cd","#mu/W^{+/-}",100,0,10);
  TH1F *level2CalDeposit_9_cd = new TH1F("level2CalDeposit_9_cd","#mu/Z",100,0,10);
  TH1F *level2CalDeposit_10_cd = new TH1F("level2CalDeposit_10_cd","#tau^{+/-}",100,0,10);
  TH1F *level2CalDeposit_11_cd = new TH1F("level2CalDeposit_11_cd","#mu^{+/-}",100,0,10);
  TH1F *level2CalDeposit_12_cd = new TH1F("level2CalDeposit_12_cd","other",100,0,10);
  TH1F *level2CalDeposit_13_cd = new TH1F("level2CalDeposit_13_cd","non-assoc",100,0,10);

  THStack *level2CalDeposit = new THStack("level2CalDeposit","CAL deposit around L2 #mu");
  THStack *level2CalDeposit_cd = new THStack("level2CalDeposit_cd","CAL deposit around L2 #mu");

  TH1F *level2SeedsL3 = new TH1F("level2SeedsL3","level2SeedsL3",2,0,2);
  // Sim level as well
  TH1F *level2PtRes = new TH1F("level2PtRes","level2PtRes",100,-1,1);
  TH1F *level2SimDeltaR = new TH1F("level2SimDeltaR","level2SimDeltaR",100,0,0.4);
  // Scatters
  TH2F *level2CalDepositEta = new TH2F("level2CalDepositEta","level2CalDepositEta",90,-3,3,100,0,10);
  TH2F *level2CalDepositPt = new TH2F("level2CalDepositPt","level2CalDepositPt",50,0,20,100,0,10);

  // Then histograms for L3 diagnositc
  TH1F *level3Number = new TH1F("level3Number","level3Number",5,0,5);
  TH1F *level3Pt = new TH1F("level3Pt","level3Pt",50,0,100);
  TH1F *level3Eta = new TH1F("level3Eta","level3Eta",90,-3,3);
  TH1F *level3Phi = new TH1F("level3Phi","level3Phi",96,-3.2,3.2);
  TH1F *level3D0 = new TH1F("level3D0","level3D0",80,-0.1,0.1);
  TH1F *level3NHits = new TH1F("level3NHits","level3NHits",70,0,70);
  TH1F *level3TrackDeposit_hist = new TH1F("level3TrackDeposit_hist","PIX deposit around L3 #mu",100,0,4);
  TH1F *level3TrackDeposit_dist_cd = new TH1F("level3TrackDeposit_dist_cd","PIX deposit around L3 #mu",100,0,4);
  // Tracking isolation
  TH1F *level3TrackDeposit_1 = new TH1F("level3TrackDeposit_1","#pi^{+/-}",100,0,4);
  TH1F *level3TrackDeposit_2 = new TH1F("level3TrackDeposit_2","K",100,0,4);
  TH1F *level3TrackDeposit_3 = new TH1F("level3TrackDeposit_3","D",100,0,4);
  TH1F *level3TrackDeposit_4 = new TH1F("level3TrackDeposit_4","B",100,0,4);
  TH1F *level3TrackDeposit_5 = new TH1F("level3TrackDeposit_5","#Lambda_{b}",100,0,4);
  TH1F *level3TrackDeposit_6 = new TH1F("level3TrackDeposit_6","J/#Psi",100,0,4);
  TH1F *level3TrackDeposit_7 = new TH1F("level3TrackDeposit_7","#Upsilon (ns)",100,0,4);
  TH1F *level3TrackDeposit_8 = new TH1F("level3TrackDeposit_8","#mu/W^{+/-}",100,0,4);
  TH1F *level3TrackDeposit_9 = new TH1F("level3TrackDeposit_9","#mu/Z",100,0,4);
  TH1F *level3TrackDeposit_10 = new TH1F("level3TrackDeposit_10","#tau^{+/-}",100,0,4);
  TH1F *level3TrackDeposit_11 = new TH1F("level3TrackDeposit_11","#mu^{+/-}",100,0,4);
  TH1F *level3TrackDeposit_12 = new TH1F("level3TrackDeposit_12","other",100,0,4);
  TH1F *level3TrackDeposit_13 = new TH1F("level3TrackDeposit_13","non-assoc",100,0,4);

  TH1F *level3TrackDeposit_1_cd = new TH1F("level3TrackDeposit_1_cd","#pi^{+/-}",100,0,4);
  TH1F *level3TrackDeposit_2_cd = new TH1F("level3TrackDeposit_2_cd","K",100,0,4);
  TH1F *level3TrackDeposit_3_cd = new TH1F("level3TrackDeposit_3_cd","D",100,0,4);
  TH1F *level3TrackDeposit_4_cd = new TH1F("level3TrackDeposit_4_cd","B",100,0,4);
  TH1F *level3TrackDeposit_5_cd = new TH1F("level3TrackDeposit_5_cd","#Lambda_{b}",100,0,4);
  TH1F *level3TrackDeposit_6_cd = new TH1F("level3TrackDeposit_6_cd","J/#Psi",100,0,4);
  TH1F *level3TrackDeposit_7_cd = new TH1F("level3TrackDeposit_7_cd","#Upsilon (ns)",100,0,4);
  TH1F *level3TrackDeposit_8_cd = new TH1F("level3TrackDeposit_8_cd","#mu/W^{+/-}",100,0,4);
  TH1F *level3TrackDeposit_9_cd = new TH1F("level3TrackDeposit_9_cd","#mu/Z",100,0,4);
  TH1F *level3TrackDeposit_10_cd = new TH1F("level3TrackDeposit_10_cd","#tau^{+/-}",100,0,4);
  TH1F *level3TrackDeposit_11_cd = new TH1F("level3TrackDeposit_11_cd","#mu^{+/-}",100,0,4);
  TH1F *level3TrackDeposit_12_cd = new TH1F("level3TrackDeposit_12_cd","other",100,0,4);
  TH1F *level3TrackDeposit_13_cd = new TH1F("level3TrackDeposit_13_cd","non-assoc",100,0,4);

  THStack *level3TrackDeposit = new THStack("level3TrackDeposit","PIX deposit around L3 #mu");
  THStack *level3TrackDeposit_cd = new THStack("level3TrackDeposit_cd","PIX deposit around L3 #mu");

  // Sim level as well
  TH1F *level3PtRes = new TH1F("level3PtRes","level3PtRes",100,-0.2,0.2);
  TH1F *level3SimDeltaR = new TH1F("level3SimDeltaR","level3SimDeltaR",100,0,0.2);
  // Scatters
  TH2F *level3TrackDepositEta = new TH2F("level3TrackDepositEta","level3TrackDepositEta",90,-3,3,100,0,4);
  TH2F *level3TrackDepositPt = new TH2F("level3TrackDepositPt","level3TrackDepositPt",50,0,20,100,0,4);
  
  int nEntries = tree->GetEntries();
  int numberL2 = 0;
  int numberL3 = 0;
  for( int iEntry = 0; iEntry < nEntries; iEntry++) {
    if (iEntry%1 == 0) cout << "Event " << iEntry << endl;
    tree->GetEntry(iEntry);

    //    cout << "Begin by testing triggers" << endl;

    // First we fill in the L1 information
    level1Number->Fill(nL1);
    bool passL1 = false;
    for (int i = 0; i < nL1; i++) {
      //      if ((*l1Pt).at(i) + 0.01 > 7 && (*l1Quality).at(i) >= 4) {
      if ((*l1Pt).at(i) + 0.01 > 0 && (*l1Quality).at(i) >= 4) {
	passL1 = true;
	continue;
      }
    }
    bool passL2 = false;
    bool passL2Iso = false;
    if (passL1) {
      for (int i = 0; i < nL2; i++) {
	//	if ((*l2Pt).at(i) > 7 && abs((*l2Eta).at(i)) < 2.5) {
	if ((*l2Pt).at(i) > 0 && abs((*l2Eta).at(i)) < 2.5) {
	  passL2 = true;
	  if ((*l2CalIsoDeposit).at(i) < etaCalIsoDepositCut((*l2Eta).at(i))) {
	    passL2Iso = true;
	    continue;
	  }
	}
      }
    }
    bool passL3 = false;
    bool passL3pre = false;
    if (passL2) {
      for (int i = 0; i < nL3; i++) {
	//        if ((*l3Pt).at(i) > 9 && fabs((*l3Eta).at(i)) < 2.5 && abs((*l3D0).at(i)) < 0.2) {
	if ((*l3Pt).at(i) > 0 && fabs((*l3Eta).at(i)) < 2.5 && abs((*l3D0).at(i)) < 0.2) {
          passL3 = true;
	  continue;
        }
      }
    }
    if (passL2Iso) {
      for (int i = 0; i < nL3; i++) {
	//        if ((*l3Pt).at(i) > 7 && fabs((*l3Eta).at(i)) < 2.5 && abs((*l3D0).at(i)) < 0.2) {
	if ((*l3Pt).at(i) > 0 && fabs((*l3Eta).at(i)) < 2.5 && abs((*l3D0).at(i)) < 0.2) {
	  passL3pre = true;
	  continue;
        }
      }
    }
    

    //    cout << "trigger check complete.  into loop for l1" << endl;

    for (int i = 0; i < nL1; i++) {
      bool triggersL1 = false;
      if ((*l1Pt).at(i) + 0.01 > 0 && (*l1Quality).at(i) >= 4) triggersL1 = true;      
      level1Pt->Fill((*l1Pt).at(i));
      level1Eta->Fill((*l1Eta).at(i));
      level1Phi->Fill((*l1Phi).at(i));
      level1Quality->Fill((*l1Quality).at(i));
      if (passL1) level1Quality_afterL1->Fill((*l1Quality).at(i));
      level1SeedsL2->Fill((*l1SeedsL2).at(i));      
      level1EtaQuality->Fill((*l1Eta).at(i),(*l1Quality).at(i));
      level1PhiQuality->Fill((*l1Phi).at(i),(*l1Quality).at(i));
      if (passL1) {
	if ((*l1SeedsL2).at(i) == 0) {
	  level1PtQuality_FailingL2Seeding->Fill((*l1Pt).at(i),(*l1Quality).at(i));
	  level1EtaQuality_FailingL2Seeding->Fill((*l1Eta).at(i),(*l1Quality).at(i));
	  level1PhiQuality_FailingL2Seeding->Fill((*l1Phi).at(i),(*l1Quality).at(i));
	  level1EtaPhi_FailingL2Seeding->Fill((*l1Eta).at(i),(*l1Phi).at(i));
	}
	else {
	  level1PtQuality_PassingL2Seeding->Fill((*l1Pt).at(i),(*l1Quality).at(i));
	  level1EtaQuality_PassingL2Seeding->Fill((*l1Eta).at(i),(*l1Quality).at(i));
	  level1PhiQuality_PassingL2Seeding->Fill((*l1Phi).at(i),(*l1Quality).at(i));
	  level1EtaPhi_PassingL2Seeding->Fill((*l1Eta).at(i),(*l1Phi).at(i));
	}
      }
      // for gen-l1 pt, eta scatters, check matching to sim
      if (nSimMuon > 0) {
	//	cout << "next candidate for a fuckup is here." << endl;
	int simIndex = -999;
	if ((*simMuonEta).size() > 0) simIndex = simIndexMatchingL1((*l1Eta).at(i),(*l1Phi).at(i),simMuonEta,simMuonPhi);
	//	cout << "made it past here" << endl;
	if (simIndex > -1) {
	  level1GenPt->Fill((*simMuonPt).at(simIndex),(*l1Pt).at(i));
	  level1GenEta->Fill((*simMuonEta).at(simIndex),(*l1Eta).at(i));
	  level1GenPhi->Fill((*simMuonPhi).at(simIndex),(*l1Phi).at(i));
	}
	//	cout << "and filled scatters successfully" << endl;
	// Sim Muon Station Numbers
	//	cout << (*simMuonMuStationNumber).begin().first << endl;
	for (map<int,vector<int> >::const_iterator fuckIt = (*simMuonMuStationNumber).begin(); fuckIt != (*simMuonMuStationNumber).end(); fuckIt++) {
	  //	  cout << "here we go into the station loop." << endl;
	  bool firstRun = false;
	  if (i == 0) firstRun = true;
	  int nStations = 0;
	  bool hitFirst = false;
	  bool hitSecond = false;
	  bool hitThird = false;
	  bool hitFourth = false;
	  for (int j = 0; j < fuckIt->second.size(); j++) { // over the hit station numbers
	    if (!hitFirst && fuckIt->second.at(j) == 1) {
	      nStations++;
	      hitFirst = true;
	    }
	    if (!hitSecond && fuckIt->second.at(j) == 2) {
	      nStations++;
	      hitSecond = true;
	    }
	    if (!hitThird && fuckIt->second.at(j) == 3) {
	      nStations++;
	      hitThird = true;
	    }
	    if (!hitFourth && fuckIt->second.at(j) == 4) {
	      nStations++;
	      hitFourth = true;
	    }
	  }
	  if (nStations > 0) {
	    if (firstRun) simStationsBeforeL1->Fill(nStations); // Do it only once, since you were muppet enough to put a loop within a loop.
	    if (fuckIt->first == simIndex && triggersL1) {
	      simStationsAfterL1->Fill(nStations);
	    } // sim muon that goes with triggering L1
	  } // nStations > 0
	} // loop over simMuonStationNumbers
	//	cout << "guessing there's a fuckup here." << endl;
	//	(*simMuonPt).erase((*simMuonPt).begin() + simIndex);
	//	(*simMuonEta).erase((*simMuonEta).begin() + simIndex);
	//	(*simMuonEta).at(simIndex) = -999.0;
	//	(*simMuonPhi).erase((*simMuonPhi).begin() + simIndex);
	//	cout << "but it made it here." << endl;
      } // nSimMuon > 0
    } // loop over L1

    //    cout << "loop over l1 complete.  something you did (the erase) fucks the code up later on." << endl;
    
    level2Number->Fill(nL2);
    for (int i = 0; i < nL2; i++) {
      // 1-D histograms
      level2Pt->Fill((*l2Pt).at(i));
      level2Eta->Fill((*l2Eta).at(i));
      level2Phi->Fill((*l2Phi).at(i));
      level2D0->Fill((*l2D0).at(i));
      level2NHits->Fill((*l2NHits).at(i));      
      level2CalDeposit_hist->Fill((*l2CalIsoDeposit).at(i));
      if (passL2) {
	numberL2++;
	if (type == 1) { // fill ppMuX parents: 1, 2, 3, 4, 5, 6, 10, 12, 13
	  if ((*l2IsAssociated).at(i) && (*l2AssociationVar).at(i) > -0.2) {
	    if ((*l2MotherBinNumber).at(i) == 1) level2CalDeposit_1->Fill((*l2CalIsoDeposit).at(i));
            else if ((*l2MotherBinNumber).at(i) == 2) level2CalDeposit_2->Fill((*l2CalIsoDeposit).at(i));
            else if ((*l2MotherBinNumber).at(i) == 3) level2CalDeposit_3->Fill((*l2CalIsoDeposit).at(i));
            else if ((*l2MotherBinNumber).at(i) == 4) level2CalDeposit_4->Fill((*l2CalIsoDeposit).at(i));
            else if ((*l2MotherBinNumber).at(i) == 5) level2CalDeposit_5->Fill((*l2CalIsoDeposit).at(i));
            else if ((*l2MotherBinNumber).at(i) == 6) level2CalDeposit_6->Fill((*l2CalIsoDeposit).at(i));
            else if ((*l2MotherBinNumber).at(i) == 10) level2CalDeposit_10->Fill((*l2CalIsoDeposit).at(i));
            else level2CalDeposit_12->Fill((*l2CalIsoDeposit).at(i));
	  }
	  else {
	    level2CalDeposit_13->Fill((*l2CalIsoDeposit).at(i));
	  }
	}  
	if (type == 2) { // fill ttbar parents: 3, 4, 8, 10, 11, 12, 13
	  if ((*l2IsAssociated).at(i) && (*l2AssociationVar).at(i) > -0.2) {
            if ((*l2MotherBinNumber).at(i) == 3) level2CalDeposit_3->Fill((*l2CalIsoDeposit).at(i));
            else if ((*l2MotherBinNumber).at(i) == 4) level2CalDeposit_4->Fill((*l2CalIsoDeposit).at(i));
            else if ((*l2MotherBinNumber).at(i) == 8) level2CalDeposit_8->Fill((*l2CalIsoDeposit).at(i));
            else if ((*l2MotherBinNumber).at(i) == 10) level2CalDeposit_10->Fill((*l2CalIsoDeposit).at(i));
            else if ((*l2MotherBinNumber).at(i) == 11) level2CalDeposit_8->Fill((*l2CalIsoDeposit).at(i));
            else level2CalDeposit_12->Fill((*l2CalIsoDeposit).at(i));
	  }
	  else {
	    level2CalDeposit_13->Fill((*l2CalIsoDeposit).at(i));
	  }
	}
	if (type == 3) { // fill ZMM parents: 9, 11, 12, 13
	  if ((*l2IsAssociated).at(i) && (*l2AssociationVar).at(i) > -0.2) {
	    if ((*l2MotherBinNumber).at(i) == 9) level2CalDeposit_9->Fill((*l2CalIsoDeposit).at(i));
	    else if ((*l2MotherBinNumber).at(i) == 11) level2CalDeposit_9->Fill((*l2CalIsoDeposit).at(i));
	    else level2CalDeposit_12->Fill((*l2CalIsoDeposit).at(i));
	  }
	  else level2CalDeposit_13->Fill((*l2CalIsoDeposit).at(i));
	}
      }
      level2SeedsL3->Fill((*l2SeedsL3).at(i));
      level2SimDeltaR->Fill(-1 * (*l2AssociationVar).at(i));
      if ((*l2AssociationVar).at(i) > -0.2) {
	level2PtRes->Fill(((*l2Pt).at(i) - (*l2AssociatedSimMuonPt).at(i))/(*l2AssociatedSimMuonPt).at(i));
      }
      // Fill the scatters
      if (passL1) {
	level2CalDepositPt->Fill((*l2Pt).at(i),(*l2CalIsoDeposit).at(i));
      }
      if (passL2) { 
	level2CalDepositEta->Fill((*l2Eta).at(i),(*l2CalIsoDeposit).at(i));
      }
    }

    //    cout << "The fuckup's not at L2, though" << endl;

    /*    if (iEntry > 115000) {
	  cout << "Diagnostic time" << endl;
	  cout << "nL3 = " << nL3 << endl;
	  }
    */

    level3Number->Fill(nL3);
    for (int i = 0; i < nL3; i++) {
      /*      if (iEntry > 115000) {
	      cout << "in the loop: i = " << i << endl;
	      cout << "variable dump: pt = " << (*l3Pt).at(i) <<endl;
	      cout << "variable dump: eta = " << (*l3Eta).at(i) <<endl;
	      cout << "variable dump: phi = " << (*l3Phi).at(i) <<endl;
	      cout << "variable dump: d0 = " << (*l3D0).at(i) <<endl;
	      cout << "variable dump: NHits = " << (*l3NHits).at(i) <<endl;
	      cout << "variable dump: trackDeposit = " << (*l3TrackIsoDeposit).at(i) <<endl;
	      cout << "variable dump: isAssociated = " << (*l3IsAssociated).at(i) << endl;
	      cout << "variable dump: AssociationVar = " << (*l3AssociationVar).at(i) <<endl;
	      }
      */
      level3Pt->Fill((*l3Pt).at(i));
      level3Eta->Fill((*l3Eta).at(i));
      level3Phi->Fill((*l3Phi).at(i));
      level3D0->Fill((*l3D0).at(i));
      level3NHits->Fill((*l3NHits).at(i));
      level3TrackDeposit_hist->Fill((*l3TrackIsoDeposit).at(i));
      if (passL3pre) {
	numberL3++;
        if (type == 1) { // fill ppMuX parents: 1, 2, 3, 4, 5, 6, 10, 12, 13
	  if ((*l3AssociationVar).at(i) > -0.1) {
	    if ((*l3MotherBinNumber).at(i) == 1) level3TrackDeposit_1->Fill((*l3TrackIsoDeposit).at(i));
            else if ((*l3MotherBinNumber).at(i) == 2) level3TrackDeposit_2->Fill((*l3TrackIsoDeposit).at(i));
            else if ((*l3MotherBinNumber).at(i) == 3) level3TrackDeposit_3->Fill((*l3TrackIsoDeposit).at(i));
            else if ((*l3MotherBinNumber).at(i) == 4) level3TrackDeposit_4->Fill((*l3TrackIsoDeposit).at(i));
            else if ((*l3MotherBinNumber).at(i) == 5) level3TrackDeposit_5->Fill((*l3TrackIsoDeposit).at(i));
            else if ((*l3MotherBinNumber).at(i) == 6) level3TrackDeposit_6->Fill((*l3TrackIsoDeposit).at(i));
            else if ((*l3MotherBinNumber).at(i) == 10) level3TrackDeposit_10->Fill((*l3TrackIsoDeposit).at(i));
	    else level3TrackDeposit_12->Fill((*l3TrackIsoDeposit).at(i));
          }
          else {
            level3TrackDeposit_13->Fill((*l3TrackIsoDeposit).at(i));
          }
        }
        if (type == 2) { // fill ttbar parents: 3, 4, 8, 10, 11, 12, 13
          if ((*l3AssociationVar).at(i) > -0.1) {
            if ((*l3MotherBinNumber).at(i) == 3) level3TrackDeposit_3->Fill((*l3TrackIsoDeposit).at(i));
            else if ((*l3MotherBinNumber).at(i) == 4) level3TrackDeposit_4->Fill((*l3TrackIsoDeposit).at(i));
            else if ((*l3MotherBinNumber).at(i) == 8) level3TrackDeposit_8->Fill((*l3TrackIsoDeposit).at(i));
            else if ((*l3MotherBinNumber).at(i) == 10) level3TrackDeposit_10->Fill((*l3TrackIsoDeposit).at(i));
            else if ((*l3MotherBinNumber).at(i) == 11) level3TrackDeposit_8->Fill((*l3TrackIsoDeposit).at(i));
            else level3TrackDeposit_12->Fill((*l3TrackIsoDeposit).at(i));
          }
          else {
            level3TrackDeposit_13->Fill((*l3TrackIsoDeposit).at(i));
          }
        }
        if (type == 3) { // fill ZMM parents: 9, 11, 12, 13
          if ((*l3AssociationVar).at(i) > -0.1) {
            if ((*l3MotherBinNumber).at(i) == 9) level3TrackDeposit_9->Fill((*l3TrackIsoDeposit).at(i));
            else if ((*l3MotherBinNumber).at(i) == 11) level3TrackDeposit_9->Fill((*l3TrackIsoDeposit).at(i));
            else level3TrackDeposit_12->Fill((*l3TrackIsoDeposit).at(i));
          }
          else {
            level3TrackDeposit_13->Fill((*l3TrackIsoDeposit).at(i));
          }
        }
      }
      level3SimDeltaR->Fill(-1 * (*l3AssociationVar).at(i));
      if ((*l3AssociationVar).at(i) > -0.1) {
        level3PtRes->Fill(((*l3Pt).at(i) - (*l3AssociatedSimMuonPt).at(i))/(*l3AssociatedSimMuonPt).at(i));
      }
      // Fill the scatters
      if (passL2Iso) {
	level3TrackDepositPt->Fill((*l3Pt).at(i),(*l3TrackIsoDeposit).at(i));
      }
      if (passL3pre) {
	level3TrackDepositEta->Fill((*l3Eta).at(i),(*l3TrackIsoDeposit).at(i));
      }
    }
    //    cout << "and we get through the L3 loop as well.  odd." << endl;
  }

  // Here we want to just normalize to the number of entries.
  make_cd_from_hist(level2CalDeposit_hist,level2CalDeposit_dist_cd,level2CalDeposit_hist->GetEntries());

  // For the individual contributions, we want to normalize to 1 by using the fraction of L2
  // muons that fall in each category, which sums up to 1.
  make_cd_from_hist(level2CalDeposit_1, level2CalDeposit_1_cd, numberL2);
  make_cd_from_hist(level2CalDeposit_2, level2CalDeposit_2_cd, numberL2);
  make_cd_from_hist(level2CalDeposit_3, level2CalDeposit_3_cd, numberL2);
  make_cd_from_hist(level2CalDeposit_4, level2CalDeposit_4_cd, numberL2);
  make_cd_from_hist(level2CalDeposit_5, level2CalDeposit_5_cd, numberL2);
  make_cd_from_hist(level2CalDeposit_6, level2CalDeposit_6_cd, numberL2);
  make_cd_from_hist(level2CalDeposit_7, level2CalDeposit_7_cd, numberL2);
  make_cd_from_hist(level2CalDeposit_8, level2CalDeposit_8_cd, numberL2);
  make_cd_from_hist(level2CalDeposit_9, level2CalDeposit_9_cd, numberL2);
  make_cd_from_hist(level2CalDeposit_10, level2CalDeposit_10_cd, numberL2);
  make_cd_from_hist(level2CalDeposit_11, level2CalDeposit_11_cd, numberL2);
  make_cd_from_hist(level2CalDeposit_12, level2CalDeposit_12_cd, numberL2);
  make_cd_from_hist(level2CalDeposit_13, level2CalDeposit_13_cd, numberL2);

  make_cd_from_hist(level3TrackDeposit_hist,level3TrackDeposit_dist_cd,level3TrackDeposit_hist->GetEntries());

  make_cd_from_hist(level3TrackDeposit_1, level3TrackDeposit_1_cd, numberL3);
  make_cd_from_hist(level3TrackDeposit_2, level3TrackDeposit_2_cd, numberL3);
  make_cd_from_hist(level3TrackDeposit_3, level3TrackDeposit_3_cd, numberL3);
  make_cd_from_hist(level3TrackDeposit_4, level3TrackDeposit_4_cd, numberL3);
  make_cd_from_hist(level3TrackDeposit_5, level3TrackDeposit_5_cd, numberL3);
  make_cd_from_hist(level3TrackDeposit_6, level3TrackDeposit_6_cd, numberL3);
  make_cd_from_hist(level3TrackDeposit_7, level3TrackDeposit_7_cd, numberL3);
  make_cd_from_hist(level3TrackDeposit_8, level3TrackDeposit_8_cd, numberL3);
  make_cd_from_hist(level3TrackDeposit_9, level3TrackDeposit_9_cd, numberL3);
  make_cd_from_hist(level3TrackDeposit_10, level3TrackDeposit_10_cd, numberL3);
  make_cd_from_hist(level3TrackDeposit_11, level3TrackDeposit_11_cd, numberL3);
  make_cd_from_hist(level3TrackDeposit_12, level3TrackDeposit_12_cd, numberL3);
  make_cd_from_hist(level3TrackDeposit_13, level3TrackDeposit_13_cd, numberL3);

  if (type == 1 ) { // 1, 2, 3, 4, 5, 6, 10, 12, 13
    level2CalDeposit->Add(level2CalDeposit_1);
    level2CalDeposit->Add(level2CalDeposit_2);
    level2CalDeposit->Add(level2CalDeposit_3);
    level2CalDeposit->Add(level2CalDeposit_4);
    level2CalDeposit->Add(level2CalDeposit_5);
    level2CalDeposit->Add(level2CalDeposit_6);
    level2CalDeposit->Add(level2CalDeposit_10);
    level2CalDeposit->Add(level2CalDeposit_12);
    level2CalDeposit->Add(level2CalDeposit_13);

    level2CalDeposit_cd->Add(level2CalDeposit_1_cd);
    level2CalDeposit_cd->Add(level2CalDeposit_2_cd);
    level2CalDeposit_cd->Add(level2CalDeposit_3_cd);
    level2CalDeposit_cd->Add(level2CalDeposit_4_cd);
    level2CalDeposit_cd->Add(level2CalDeposit_5_cd);
    level2CalDeposit_cd->Add(level2CalDeposit_6_cd);
    level2CalDeposit_cd->Add(level2CalDeposit_10_cd);
    level2CalDeposit_cd->Add(level2CalDeposit_12_cd);
    level2CalDeposit_cd->Add(level2CalDeposit_13_cd);

    level3TrackDeposit->Add(level3TrackDeposit_1);
    level3TrackDeposit->Add(level3TrackDeposit_2);
    level3TrackDeposit->Add(level3TrackDeposit_3);
    level3TrackDeposit->Add(level3TrackDeposit_4);
    level3TrackDeposit->Add(level3TrackDeposit_5);
    level3TrackDeposit->Add(level3TrackDeposit_6);
    level3TrackDeposit->Add(level3TrackDeposit_10);
    level3TrackDeposit->Add(level3TrackDeposit_12);
    level3TrackDeposit->Add(level3TrackDeposit_13);

    level3TrackDeposit_cd->Add(level3TrackDeposit_1_cd);
    level3TrackDeposit_cd->Add(level3TrackDeposit_2_cd);
    level3TrackDeposit_cd->Add(level3TrackDeposit_3_cd);
    level3TrackDeposit_cd->Add(level3TrackDeposit_4_cd);
    level3TrackDeposit_cd->Add(level3TrackDeposit_5_cd);
    level3TrackDeposit_cd->Add(level3TrackDeposit_6_cd);
    level3TrackDeposit_cd->Add(level3TrackDeposit_10_cd);
    level3TrackDeposit_cd->Add(level3TrackDeposit_12_cd);
    level3TrackDeposit_cd->Add(level3TrackDeposit_13_cd);
  }
  
  if (type == 2 ) { // 3, 4, 8, 10, 12, 13
    level2CalDeposit->Add(level2CalDeposit_3);
    level2CalDeposit->Add(level2CalDeposit_4);
    level2CalDeposit->Add(level2CalDeposit_8);
    level2CalDeposit->Add(level2CalDeposit_10);
    level2CalDeposit->Add(level2CalDeposit_12);
    level2CalDeposit->Add(level2CalDeposit_13);

    level2CalDeposit_cd->Add(level2CalDeposit_3_cd);
    level2CalDeposit_cd->Add(level2CalDeposit_4_cd);
    level2CalDeposit_cd->Add(level2CalDeposit_8_cd);
    level2CalDeposit_cd->Add(level2CalDeposit_10_cd);
    level2CalDeposit_cd->Add(level2CalDeposit_12_cd);
    level2CalDeposit_cd->Add(level2CalDeposit_13_cd);

    level3TrackDeposit->Add(level3TrackDeposit_3);
    level3TrackDeposit->Add(level3TrackDeposit_4);
    level3TrackDeposit->Add(level3TrackDeposit_8);
    level3TrackDeposit->Add(level3TrackDeposit_10);
    level3TrackDeposit->Add(level3TrackDeposit_12);
    level3TrackDeposit->Add(level3TrackDeposit_13);

    level3TrackDeposit_cd->Add(level3TrackDeposit_3_cd);
    level3TrackDeposit_cd->Add(level3TrackDeposit_4_cd);
    level3TrackDeposit_cd->Add(level3TrackDeposit_8_cd);
    level3TrackDeposit_cd->Add(level3TrackDeposit_10_cd);
    level3TrackDeposit_cd->Add(level3TrackDeposit_12_cd);
    level3TrackDeposit_cd->Add(level3TrackDeposit_13_cd);
  }
  
  if (type == 3 ) { // 9, 12, 13
    level2CalDeposit->Add(level2CalDeposit_9);
    level2CalDeposit->Add(level2CalDeposit_12);
    level2CalDeposit->Add(level2CalDeposit_13);

    level2CalDeposit_cd->Add(level2CalDeposit_9_cd);
    level2CalDeposit_cd->Add(level2CalDeposit_12_cd);
    level2CalDeposit_cd->Add(level2CalDeposit_13_cd);

    level3TrackDeposit->Add(level3TrackDeposit_9);
    level3TrackDeposit->Add(level3TrackDeposit_12);
    level3TrackDeposit->Add(level3TrackDeposit_13);

    level3TrackDeposit_cd->Add(level3TrackDeposit_9_cd);
    level3TrackDeposit_cd->Add(level3TrackDeposit_12_cd);
    level3TrackDeposit_cd->Add(level3TrackDeposit_13_cd);
  }

  gDirectory->Append(level2CalDeposit);
  gDirectory->Append(level2CalDeposit_cd);
  gDirectory->Append(level3TrackDeposit);
  gDirectory->Append(level3TrackDeposit_cd);

  histFile->Write("",TObject::kOverwrite);
  return 0;
}

int simIndexMatchingL1 (double l1RecoEta, double l1RecoPhi, vector<double> *simEta, vector<double> *simPhi) {
  int index = -1;
  double DeltaR = 999.9;

  //  cout << "testing this method...looping over sim" << endl;
  //  cout << "simEta size" << (*simEta).size() << endl;

  for (int i = 0; i < (*simEta).size(); i++) {
    //    cout << "eta phi at index" << (*simEta).at(i) <<" "<< (*simPhi).at(i) << endl;
    double temp_deleta = l1RecoEta - (*simEta).at(i);
    double temp_delphi = l1RecoPhi - (*simPhi).at(i);
    double temp_DeltaR = sqrt((temp_deleta * temp_deleta) + (temp_delphi * temp_delphi));
    if (temp_DeltaR < DeltaR  && temp_DeltaR < 1) {
      index = i;
      DeltaR = temp_DeltaR;
    }
  }
  //  cout << "DeltaR = " << DeltaR << endl;
  //  cout << "returning" << endl;
  return index;
}

double etaCalIsoDepositCut (double eta) {
  if (abs(eta) < 0.435) return 4.0;
  else if (abs(eta) < 0.1305) return 3.7;
  else if (abs(eta) < 0.2175) return 4.0;
  else if (abs(eta) < 0.3045) return 3.5;
  else if (abs(eta) < 0.3915) return 3.4;
  else if (abs(eta) < 0.4785) return 3.4;
  else if (abs(eta) < 0.5655) return 3.2;
  else if (abs(eta) < 0.6525) return 3.4;
  else if (abs(eta) < 0.7395) return 3.1;
  else if (abs(eta) < 0.8265) return 2.9;
  else if (abs(eta) < 0.9135) return 2.9;
  else if (abs(eta) < 1.0005) return 2.7;
  else if (abs(eta) < 1.0875) return 3.1;
  else if (abs(eta) < 1.1745) return 3.0;
  else if (abs(eta) < 1.2615) return 2.4;
  else if (abs(eta) < 1.3485) return 2.1;
  else if (abs(eta) < 1.4355) return 2.0;
  else if (abs(eta) < 1.5225) return 2.3;
  else if (abs(eta) < 1.6095) return 2.2;
  else if (abs(eta) < 1.6965) return 2.4;
  else if (abs(eta) < 1.785) return 2.5;
  else if (abs(eta) < 1.88) return 2.5;
  else if (abs(eta) < 1.9865) return 2.6;
  else if (abs(eta) < 2.1075) return 2.9;
  else if (abs(eta) < 2.247) return 3.1;
  else if (abs(eta) < 2.411) return 2.9;
  else return 2.9;   
}

// This one takes a hist and makes a cd.  Doesn't have to be same number of bins, but the number of bins should be the same
void make_cd_from_hist(TH1F* hist, TH1F* cd, double norm_factor) {

  int nBins1 = hist->GetNbinsX();
  int nBins2 = cd->GetNbinsX();
  if (hist->GetBinWidth(1) != cd->GetBinWidth(1)) {
    cout << "This is wrong." << endl;
    cout << "width 1, 2 = " << hist->GetBinWidth(1) << " " << cd->GetBinWidth(1) << endl;
    return;
  }

  double content = 0;
  for(int i=1; i<=nBins2; ++i){
    for (int j=i; j<=nBins1 + 1; ++j) { // include overflow
      content+=hist->GetBinContent(j);
    }
    cd->SetBinContent(i, content); //Integrate.
    content = 0; // Reset it, you muppet.
  }

  //  double norm_factor = hist->GetEntries();
   
  if (norm_factor != 0) {
    cd->Scale(1/norm_factor);
  }

  return;
}

