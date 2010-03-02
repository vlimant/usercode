#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
//#include "Math/LorentzVector.h"
#include <vector>
#include <map>
#include <utility>
#include "CMS1.h"
#include "Style.C"

setTDRStyle();

//int ScanTree ( TTree* tree) {
int ScanTree ( TTree* tree, char *fileName) {

  //  double pi = 3.14159265358979323846;
  
  Init(tree);
  TFile *histFile = new TFile(fileName,"RECREATE");

  // Testing to see what L1 really does
  //  int nPtBins = 33;
  //  Double_t pt_Edges[34] = {0, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42.5, 45, 47.5, 50, 55, 60, 70, 80, 100};
  
  int nPtBins = 23;
  double pt_Edges[24] = {0, 5, 5.25, 5.5, 5.75, 6, 6.25, 6.5, 6.75, 7, 7.25, 7.5, 7.75, 8, 8.5, 9, 10, 11, 12, 14, 16, 20, 25, 30};
  TH1F *level1Pt =  new TH1F("level1Pt","level1Pt",150,0,150);
  TH1F *level1Eta = new TH1F("level1Eta","level1Eta",100,-2.5,2.5);
  TH1F *level1Phi = new TH1F("level1Phi","level1Phi",144,-3.18,3.1);
  TH1F *level1Quality = new TH1F("level1Quality","level1Quality",8,0,8);
  TH1F *level1Quality_afterL1 = new TH1F("level1Quality_afterL1","level1Quality_afterL1",8,0,8);

  TH1F *level1LowPtNum = new TH1F("level1LowPtNum","level1LowPtNum",nPtBins,pt_Edges);
  TH1F *level1PtNum = new TH1F("level1PtNum","level1PtNum",100,0,500);
  TH1F *level1EtaNum = new TH1F("level1EtaNum","level1EtaNum",36,-2.4,2.4);
  TH1F *level1PhiNum = new TH1F("level1PhiNum","level1PhiNum",32,-3.2,3.2);

  TH1F *level1LowPtDenom = new TH1F("level1LowPtDenom","level1LowPtDenom",nPtBins,pt_Edges);
  TH1F *level1PtDenom = new TH1F("level1PtDenom","level1PtDenom",100,0,500);
  TH1F *level1EtaDenom = new TH1F("level1EtaDenom","level1EtaDenom",36,-2.4,2.4);
  TH1F *level1PhiDenom = new TH1F("level1PhiDenom","level1PhiDenom",32,-3.2,3.2);

  TH1F *level1LowPtEff = new TH1F("level1LowPtEff","level1LowPtEff",nPtBins,pt_Edges);
  TH1F *level1PtEff = new TH1F("level1PtEff","level1PtEff",100,0,500);
  TH1F *level1EtaEff = new TH1F("level1EtaEff","level1EtaEff",36,-2.4,2.4);
  TH1F *level1PhiEff = new TH1F("level1PhiEff","level1PhiEff",32,-3.2,3.2);

  // efficiency maps
  TH2F *level1PtEtaNum = new TH2F("level1PtEtaNum","level1PtEtaNum",100,0,500,36,-2.4,2.4);
  TH2F *level1PtEtaDenom = new TH2F("level1PtEtaDenom","level1PtEtaDenom",100,0,500,36,-2.4,2.4);
  TH2F *level1PtEtaEff = new TH2F("level1PtEtaEff","level1PtEtaEff",100,0,500,36,-2.4,2.4);

  // This turned out to work so well that I'm going to use this script for scatters and resolutions too
  TH2F *level1PtEta = new TH2F("level1PtEta","level1PtEta",150,0,150,100,-2.5,2.5);
  TH2F *level1EtaQuality = new TH2F("level1EtaQuality","level1EtaQuality",100,-2.5,2.5,8,0,8);
  TH2F *level1PhiQuality = new TH2F("level1PhiQuality","level1PhiQuality",144,-3.18,3.1,8,0,8);
  
  TH2F *level1GenPt = new TH2F("level1GenPt","level1GenPt",100,0,200,150,0,150); // Binned for single mu 0 - 500
  TH2F *level1GenPt_afterL1 = new TH2F("level1GenPt_afterL1","level1GenPt_afterL1",100,0,200,150,0,150); // Binned for single mu 0 - 500
  TH2F *level1GenForwardPt = new TH2F("level1GenForwardPt","level1GenForwardPt",100,0,200,150,0,150); // Binned for single mu 0 - 500
  TH2F *level1GenEta = new TH2F("level1GenEta","level1GenEta",100,-2.5,2.5,100,-2.5,2.5);
  TH2F *level1GenPhi = new TH2F("level1GenPhi","level1GenPhi",144,-3.18,3.1,160,-3.2,3.2);

  TH1F *level1SimDeltaR = new TH1F("level1SimDeltaR","level1SimDeltaR",100,0,0.4);

  TH1F *level1ptRes = new TH1F("level1ptRes","level1ptRes",100,-1,1);  // binned exclusively for TeV.  100,-1,1 otherwise
  TH2F *level1ptRes_vs_simpT = new TH2F("level1ptRes_vs_simpT","level1ptRes_vs_simpT",100,0,200,100,-1,1);
  TH2F *level1ptRes_vs_simEta = new TH2F("level1ptRes_vs_simEta","level1ptRes_vs_simEta",100,-2.5,2.5,100,-1,1);

  // Stations before/after L1 Filter
  TH1F *simStationsBeforeL1 = new TH1F("simStationsBeforeL1","simStationsBeforeL1",5,0,5);
  TH1F *simStationsAfterL1 = new TH1F("simStationsAfterL1","simStationsAfterL1",5,0,5);

  TH2F *simStationsPt = new TH2F("simStationsPt","simStationsPt",100,0,200,5,0,5);
  TH2F *simStationsEta = new TH2F("simStationsEta","simStationsEta",48,-2.4,2.4,5,0,5);

  TH2F *simHitsPt = new TH2F("simHitsPt","simHitsPt",100,0,200,60,0,60);
  TH2F *simHitsEta = new TH2F("simHitsEta","simHitsEta",48,-2.4,2.4,60,0,60);

  int nEntries = tree->GetEntries();
  for( int iEntry = 0; iEntry < nEntries; iEntry++) {
    if (iEntry%1000 == 0) cout << "Event " << iEntry << endl;
    tree->GetEntry(iEntry);

    // First we fill in the L1 information
    int l1PassCount = 0; // keep track so we don't double-count cases where > 2 l1s pass the filter
    //    cout << "loop over sim" << endl;
    //    cout << "sim muon number:" << nSimMuon << endl;
    for (int i = 0; i < nSimMuon; i++) {
      //      cout << "sim muon pt eta phi" << (*simMuonPt).at(i) << " " << (*simMuonEta).at(i) << " " << (*simMuonPhi).at(i) << endl;
      if ((*simMuonPt).at(i) > 5 && (*simMuonEta).at(i) > -2.1 && (*simMuonEta).at(i) < 2.1) {
        level1LowPtDenom->Fill((*simMuonPt).at(i));
        level1PtDenom->Fill((*simMuonPt).at(i));
        level1EtaDenom->Fill((*simMuonEta).at(i));
        level1PhiDenom->Fill((*simMuonPhi).at(i));
	level1PtEtaDenom->Fill((*simMuonPt).at(i),(*simMuonEta).at(i));
	/*	for (map<int,int>::const_iterator fuckIt = (*simMuonNMuHits).begin(); fuckIt != (*simMuonNMuHits).end(); fuckIt++) {
		if (fuckIt->first == i) {
		if (fuckIt->second > 0) {
		simHitsPt->Fill((*simMuonPt).at(i),fuckIt->second);
		simHitsEta->Fill((*simMuonEta).at(i),fuckIt->second);
		}
		}
		}
	*/
      }
    }
    //    cout << "loop over L1" << endl;
    for (int i = 0; i < nL1; i++) {
      //      cout << "nL1 = " << nL1 << endl;
      level1Pt->Fill((*l1Pt).at(i));
      level1Eta->Fill((*l1Eta).at(i));
      level1Phi->Fill((*l1Phi).at(i));
      level1PtEta->Fill((*l1Pt).at(i),(*l1Eta).at(i));
      level1Quality->Fill((*l1Quality).at(i));
      level1EtaQuality->Fill((*l1Eta).at(i),(*l1Quality).at(i));
      level1PhiQuality->Fill((*l1Phi).at(i),(*l1Quality).at(i));
      if (nSimMuon > 0) {
	// make sure we have sim to associate to reco
        int simIndex = simIndexMatchingL1((*l1Eta).at(i),(*l1Phi).at(i),simMuonEta,simMuonPhi);
        if (simIndex > -1) { // make sure we did have association
	  if ((*simMuonPt).at(simIndex) < 100) { // remove for full range
	    level1ptRes->Fill(((*l1Pt).at(i) - (*simMuonPt).at(simIndex)) / (*simMuonPt).at(simIndex));
	    level1ptRes_vs_simpT->Fill((*simMuonPt).at(simIndex),((*l1Pt).at(i) - (*simMuonPt).at(simIndex)) / (*simMuonPt).at(simIndex));	    
	    level1ptRes_vs_simEta->Fill((*simMuonEta).at(simIndex),((*l1Pt).at(i) - (*simMuonPt).at(simIndex)) / (*simMuonPt).at(simIndex));
	  }
	  level1GenPt->Fill((*simMuonPt).at(simIndex),(*l1Pt).at(i));
	  if ((*l1Eta).at(i) > 2.1 || (*l1Eta).at(i) < -2.1) level1GenForwardPt->Fill((*simMuonPt).at(simIndex),(*l1Pt).at(i));
	  level1GenEta->Fill((*simMuonEta).at(simIndex),(*l1Eta).at(i));
	  level1GenPhi->Fill((*simMuonPhi).at(simIndex),(*l1Phi).at(i));

	  double delEta = (*l1Eta).at(i) - (*simMuonEta).at(simIndex);
	  double delPhi = (*l1Phi).at(i) - (*simMuonPhi).at(simIndex);
	  double delR = sqrt((delEta * delEta) + (delPhi * delPhi));
	  level1SimDeltaR->Fill(delR);

	  bool triggersL1 = false;
	  if ((*l1Pt).at(i) + 0.01 > 7 && (*l1Quality).at(i) >= 4 && (*l1SeedsL2).at(i)) { // l1 triggering and seeding
	    triggersL1 = true;
	    level1Quality_afterL1->Fill((*l1Quality).at(i));
	    level1GenPt_afterL1->Fill((*simMuonPt).at(simIndex),(*l1Pt).at(i));
	    //	cout << "l1 passes" << endl;
	    l1PassCount++;
	    int simIndex = simIndexMatchingL1((*l1Eta).at(i),(*l1Phi).at(i),simMuonEta,simMuonPhi);
	    //       	cout << "simIndex = " << simIndex << endl;
	    if ((*simMuonPt).at(simIndex) > 5 && (*simMuonEta).at(simIndex) > -2.1 && (*simMuonEta).at(simIndex) < 2.1) {	    
	      level1LowPtNum->Fill((*simMuonPt).at(simIndex));
	      level1PtNum->Fill((*simMuonPt).at(simIndex));
	      level1EtaNum->Fill((*simMuonEta).at(simIndex));
	      level1PhiNum->Fill((*simMuonPhi).at(simIndex));	
	      level1PtEtaNum->Fill((*simMuonPt).at(simIndex),(*simMuonEta).at(simIndex));
	    }
	  }
	  /*	  for (map<int,vector<int> >::const_iterator fuckIt2 = (*simMuonMuStationNumber).begin(); fuckIt2 != (*simMuonMuStationNumber).end(); fuckIt2++) {
	    //      cout << "here we go into the station loop." << endl;
	    bool firstRun = false;
	    if (i == 0) firstRun = true;
	    int nStations = 0;
	    bool hitFirst = false;
	    bool hitSecond = false;
	    bool hitThird = false;
	    bool hitFourth = false;
	    for (int j = 0; j < fuckIt2->second.size(); j++) { // over the hit station numbers
	      if (!hitFirst && fuckIt2->second.at(j) == 1) {
		nStations++;
		hitFirst = true;
	      }
	      if (!hitSecond && fuckIt2->second.at(j) == 2) {
		nStations++;
		hitSecond = true;
	      }
	      if (!hitThird && fuckIt2->second.at(j) == 3) {
		nStations++;
		hitThird = true;
	      }
	      if (!hitFourth && fuckIt2->second.at(j) == 4) {
		nStations++;
		hitFourth = true;
	      }
	    }
	    if (nStations > 0) {
	      if (firstRun) { 
		//		cout << "filling stations" << endl;
		simStationsBeforeL1->Fill(nStations); // Do it only once, since you were muppet enough to put a loop within a loop.
		//		cout << "fuckIt2.first, nStations = " << fuckIt2->first <<" " << nStations << endl;
		simStationsPt->Fill((*simMuonPt).at(fuckIt2->first),nStations);
		simStationsEta->Fill((*simMuonEta).at(fuckIt2->first),nStations);
		//		cout << "stations filled" << endl;
	      }
	      if (fuckIt2->first == simIndex && triggersL1) simStationsAfterL1->Fill(nStations);
	    } // nStations > 0
	  } // loop over simMuonStationNumbers
	  */
	  (*simMuonPt).erase((*simMuonPt).begin() + simIndex);
	  (*simMuonEta).erase((*simMuonEta).begin() + simIndex);
	  (*simMuonPhi).erase((*simMuonPhi).begin() + simIndex);
	}
      }
    }
  }
  divide_histos_and_errors(level1LowPtNum,level1LowPtDenom,level1LowPtEff);
  divide_histos_and_errors(level1PtNum,level1PtDenom,level1PtEff);
  divide_histos_and_errors(level1EtaNum,level1EtaDenom,level1EtaEff);
  divide_histos_and_errors(level1PhiNum,level1PhiDenom,level1PhiEff);
  divide_TH2_histos_and_errors(level1PtEtaNum,level1PtEtaDenom,level1PtEtaEff);
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

void divide_histos_and_errors(TH1F* HIS1, TH1F* HIS2, TH1F* HIS3) {
  // HIS1 = numerator
  // HIS2 = denominator
  // HIS3 = ratio
  Int_t nbins1 = HIS1->GetNbinsX(); // find out nr of bins to loop over.
  for(Int_t ibins1=1 ; ibins1 <= nbins1 ; ibins1 ++) {
    Float_t content_1 = HIS1->GetBinContent(ibins1);
    Float_t content_2 = HIS2->GetBinContent(ibins1);
    Float_t error_1 = HIS1->GetBinError(ibins1);
    Float_t error_2 = HIS2->GetBinError(ibins1);

    if (content_2 != 0) {
      HIS3->SetBinContent(ibins1, content_1/content_2);
      HIS3->SetBinError(ibins1, error_1/content_2);
    }
  }
}

void divide_TH2_histos_and_errors(TH2F* HIS1, TH2F* HIS2, TH2F* HIS3) {
  // HIS1 = numerator
  // HIS2 = denominator
  // HIS3 = ratio
  Int_t nbinsX = HIS1->GetNbinsX(); // find out nr of bins to loop over.
  Int_t nbinsY = HIS1->GetNbinsY();
  for(Int_t ibinsX=1; ibinsX <= nbinsX; ibinsX ++) {
    for (Int_t ibinsY=1; ibinsY <= nbinsY; ibinsY ++) {
      Float_t content_1 = HIS1->GetBinContent(ibinsX, ibinsY);
      Float_t content_2 = HIS2->GetBinContent(ibinsX, ibinsY);
      Float_t error_1 = HIS1->GetBinError(ibinsX, ibinsY);
      Float_t error_2 = HIS2->GetBinError(ibinsX, ibinsY);

      if (content_2 != 0) {
        HIS3->SetBinContent(ibinsX, ibinsY, content_1/content_2);
        HIS3->SetBinError(ibinsX, ibinsY, error_1/content_2);
      }
      else {
        HIS3->SetBinContent(ibinsX, ibinsY, 0);
        HIS3->SetBinError(ibinsX, ibinsY, 0);
      }
    }
  }
}
