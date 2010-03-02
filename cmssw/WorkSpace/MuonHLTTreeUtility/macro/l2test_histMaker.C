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
  TH1F *level2Pt =  new TH1F("level2Pt","level2Pt",150,0,750);
  TH1F *level2Eta = new TH1F("level2Eta","level2Eta",100,-2.5,2.5);
  TH1F *level2Phi = new TH1F("level2Phi","level2Phi",144,-3.18,3.1);

  TH1F *level2LowPtNum = new TH1F("level2LowPtNum","level2LowPtNum",nPtBins,pt_Edges);
  TH1F *level2PtNum = new TH1F("level2PtNum","level2PtNum",100,0,500);
  TH1F *level2EtaNum = new TH1F("level2EtaNum","level2EtaNum",36,-2.4,2.4);
  TH1F *level2PhiNum = new TH1F("level2PhiNum","level2PhiNum",32,-3.2,3.2);

  TH1F *level2OverSimLowPtDenom = new TH1F("level2OverSimLowPtDenom","level2OverSimLowPtDenom",nPtBins,pt_Edges);
  TH1F *level2OverSimPtDenom = new TH1F("level2OverSimPtDenom","level2OverSimPtDenom",100,0,500);
  TH1F *level2OverSimEtaDenom = new TH1F("level2OverSimEtaDenom","level2OverSimEtaDenom",36,-2.4,2.4);
  TH1F *level2OverSimPhiDenom = new TH1F("level2OverSimPhiDenom","level2OverSimPhiDenom",32,-3.2,3.2);

  TH1F *level2OverSimLowPtEff = new TH1F("level2OverSimLowPtEff","level2OverSimLowPtEff",nPtBins,pt_Edges);
  TH1F *level2OverSimPtEff = new TH1F("level2OverSimPtEff","level2OverSimPtEff",100,0,500);
  TH1F *level2OverSimEtaEff = new TH1F("level2OverSimEtaEff","level2OverSimEtaEff",36,-2.4,2.4);
  TH1F *level2OverSimPhiEff = new TH1F("level2OverSimPhiEff","level2OverSimPhiEff",32,-3.2,3.2);

  TH1F *level2OverL1LowPtDenom = new TH1F("level2OverL1LowPtDenom","level2OverL1LowPtDenom",nPtBins,pt_Edges);
  TH1F *level2OverL1PtDenom = new TH1F("level2OverL1PtDenom","level2OverL1PtDenom",100,0,500);
  TH1F *level2OverL1EtaDenom = new TH1F("level2OverL1EtaDenom","level2OverL1EtaDenom",36,-2.4,2.4);
  TH1F *level2OverL1PhiDenom = new TH1F("level2OverL1PhiDenom","level2OverL1PhiDenom",32,-3.2,3.2);

  TH1F *level2OverL1LowPtEff = new TH1F("level2OverL1LowPtEff","level2OverL1LowPtEff",nPtBins,pt_Edges);
  TH1F *level2OverL1PtEff = new TH1F("level2OverL1PtEff","level2OverL1PtEff",100,0,500);
  TH1F *level2OverL1EtaEff = new TH1F("level2OverL1EtaEff","level2OverL1EtaEff",36,-2.4,2.4);
  TH1F *level2OverL1PhiEff = new TH1F("level2OverL1PhiEff","level2OverL1PhiEff",32,-3.2,3.2);

  // efficiency maps
  TH2F *level2PtEtaNum = new TH2F("level2PtEtaNum","level2PtEtaNum",100,0,500,36,-2.4,2.4);

  TH2F *level2OverSimPtEtaDenom = new TH2F("level2OverSimPtEtaDenom","level2OverSimPtEtaDenom",100,0,500,36,-2.4,2.4);
  TH2F *level2OverSimPtEtaEff = new TH2F("level2OverSimPtEtaEff","level2OverSimPtEtaEff",100,0,500,36,-2.4,2.4);

  TH2F *level2OverL1PtEtaDenom = new TH2F("level2OverL1PtEtaDenom","level2OverL1PtEtaDenom",100,0,500,36,-2.4,2.4);
  TH2F *level2OverL1PtEtaEff = new TH2F("level2OverL1PtEtaEff","level2OverL1PtEtaEff",100,0,500,36,-2.4,2.4);

  // This turned out to work so well that I'm going to use this script for scatters and resolutions too
  TH2F *level2PtEta = new TH2F("level2PtEta","level2PtEta",150,0,750,100,-2.5,2.5);
  
  TH2F *level2GenPt = new TH2F("level2GenPt","level2GenPt",100,0,500,150,0,750); // Binned for single mu 0 - 500
  TH2F *level2GenPt_afterL2 = new TH2F("level2GenPt_afterL2","level2GenPt_afterL2",100,0,500,150,0,750); // Binned for single mu 0 - 500
  TH2F *level2GenEta = new TH2F("level2GenEta","level2GenEta",100,-2.5,2.5,100,-2.5,2.5);
  TH2F *level2GenPhi = new TH2F("level2GenPhi","level2GenPhi",144,-3.18,3.1,160,-3.2,3.2);

  TH1F *level2SimDeltaR = new TH1F("level2SimDeltaR","level2SimDeltaR",100,0,0.4);

  TH1F *level2ptRes = new TH1F("level2ptRes","level2ptRes",100,-1,1);  // binned exclusively for TeV.  100,-1,1 otherwise
  TH2F *level2ptRes_vs_simpT = new TH2F("level2ptRes_vs_simpT","level2ptRes_vs_simpT",100,0,500,100,-1,1);
  TH2F *level2ptRes_vs_simEta = new TH2F("level2ptRes_vs_simEta","level2ptRes_vs_simEta",100,-2.5,2.5,100,-1,1);

  // Stations before/after L2 Filter
  TH1F *l2StationsBeforeL2 = new TH1F("l2StationsBeforeL2","l2StationsBeforeL2",5,0,5);
  TH1F *l2StationsAfterL2 = new TH1F("l2StationsAfterL2","l2StationsAfterL2",5,0,5);

  TH2F *l2StationsPt = new TH2F("l2StationsPt","l2StationsPt",100,0,500,5,0,5);
  TH2F *l2StationsEta = new TH2F("l2StationsEta","l2StationsEta",48,-2.4,2.4,5,0,5);

  TH2F *simHitsPt = new TH2F("simHitsPt","simHitsPt",100,0,500,60,0,60);
  TH2F *simHitsEta = new TH2F("simHitsEta","simHitsEta",48,-2.4,2.4,60,0,60);

  int nEntries = tree->GetEntries();
  for( int iEntry = 0; iEntry < nEntries; iEntry++) {
    if (iEntry%1000 == 0) cout << "Event " << iEntry << endl;
    tree->GetEntry(iEntry);
    
    // First we fill in the sim information
    //    cout << "Loop over Sim" << endl;
    for (int i = 0; i < nSimMuon; i++) {
      //      cout << "sim muon pt eta phi" << (*simMuonPt).at(i) << " " << (*simMuonEta).at(i) << " " << (*simMuonPhi).at(i) << endl;
      if ((*simMuonPt).at(i) > 5 && (*simMuonEta).at(i) > -2.1 && (*simMuonEta).at(i) < 2.1) {
        level2OverSimLowPtDenom->Fill((*simMuonPt).at(i));
        level2OverSimPtDenom->Fill((*simMuonPt).at(i));
        level2OverSimEtaDenom->Fill((*simMuonEta).at(i));
        level2OverSimPhiDenom->Fill((*simMuonPhi).at(i));
	level2OverSimPtEtaDenom->Fill((*simMuonPt).at(i),(*simMuonEta).at(i));
	bool foundL2 = false;
	for (int iL2 = 0; iL2 < nL2; iL2++) {
	  if ((*simToL2Associated).at(i) && (*simToL2RecoIndex).at(i) == iL2 && (*simToL2AssociationVar).at(i) > -0.2) {
	    if (foundL2) cout << "Double-counting, you muppet" << endl;
	    foundL2 = true;
	    if ((*l1Pt).at((*indexL1SeedingL2).at(iL2)) + 0.01 > 7 && (*l1Quality).at((*indexL1SeedingL2).at(iL2)) >= 4) { 
	      level2OverL1LowPtDenom->Fill((*simMuonPt).at(i));
	      level2OverL1PtDenom->Fill((*simMuonPt).at(i));
	      level2OverL1EtaDenom->Fill((*simMuonEta).at(i));
	      level2OverL1PhiDenom->Fill((*simMuonPhi).at(i));
	      level2OverL1PtEtaDenom->Fill((*simMuonPt).at(i),(*simMuonEta).at(i));
	      if ((*l2Pt).at(iL2) > 7 && (*l2Eta).at(iL2) > -2.5 && (*l2Eta).at(iL2) < 2.5) { // does this L2 pass the filter requirement?		
		level2LowPtNum->Fill((*simMuonPt).at(i));
		level2PtNum->Fill((*simMuonPt).at(i));
		level2EtaNum->Fill((*simMuonEta).at(i));
		level2PhiNum->Fill((*simMuonPhi).at(i));
		level2PtEtaNum->Fill((*simMuonPt).at(i),(*simMuonEta).at(i));
	      } // l2 muon triggers
	      /*	      else { 
			      cout << "L1 pass, L2 fail: " << endl;
			      cout << "l1 quantities" << (*l1Pt).at((*indexL1SeedingL2).at(iL2)) << " " << (*l1Quality).at((*indexL1SeedingL2).at(iL2)) <<" "<< (*l1Eta).at((*indexL1SeedingL2).at(iL2)) << " " << (*l1Phi).at((*indexL1SeedingL2).at(iL2)) << endl;
			      cout << "L2 Quantities" << (*l2Pt).at(iL2) << " " << (*l2Eta).at(iL2) << " " << (*l2Phi).at(iL2) << endl;
	      
			      }
	      */
	    } // l1 seed triggers
	    //	    else cout << "sim pass, L1 fail" << endl;
	  } // sim to l2 associated
	} // loop over l2
	foundL2 = false;
      } // sim muon selection
    } // loop over sim
    
    for (int i = 0; i < nL2; i++) {
      level2Pt->Fill((*l2Pt).at(i));
      level2Eta->Fill((*l2Eta).at(i));
      level2Phi->Fill((*l2Phi).at(i));
      level2PtEta->Fill((*l2Pt).at(i),(*l2Eta).at(i));
      level2SimDeltaR->Fill(-1 * (*l2AssociationVar).at(i));
      //      cout << "Checking L2 association " << (*l2IsAssociated).at(i) <<" "<< (*l2AssociationVar).at(i) << endl;
      if ((*l2IsAssociated).at(i) && (*l2AssociationVar).at(i) > -0.2) {
	//	cout << "pt resolution" << endl;
	level2ptRes->Fill(((*l2Pt).at(i) - (*l2AssociatedSimMuonPt).at(i)) / (*l2AssociatedSimMuonPt).at(i));
	//	cout << "vs pT" << endl;
	level2ptRes_vs_simpT->Fill((*l2AssociatedSimMuonPt).at(i),((*l2Pt).at(i) - (*l2AssociatedSimMuonPt).at(i)) / (*l2AssociatedSimMuonPt).at(i));
	//	cout << "vs eta" << endl;
	level2ptRes_vs_simEta->Fill((*l2AssociatedSimMuonEta).at(i),((*l2Pt).at(i) - (*l2AssociatedSimMuonPt).at(i)) / (*l2AssociatedSimMuonPt).at(i));
	//	cout << "pT scatter" << endl;
	level2GenPt->Fill((*l2AssociatedSimMuonPt).at(i),(*l2Pt).at(i));
	//	cout << "eta scatter" << endl;
	level2GenEta->Fill((*l2AssociatedSimMuonEta).at(i),(*l2Eta).at(i));
	//	cout << "phi scatter" << endl;
	level2GenPhi->Fill((*l2AssociatedSimMuonPhi).at(i),(*l2Phi).at(i));	  
      } // l2 to sim associated
      //      cout << "Checking L2 trigger pass" << endl;
      if ((*l2Pt).at(i) > 7 && (*l2Eta).at(i) > -2.5 && (*l2Eta).at(i) < 2.5) { // does this L2 pass the filter requirement?
	// This is also the trigger efficiency per muon.  That means we also want to know if the L2's L1 seed passed.
	//	cout << "Checking L1 seed pass" << endl;
	if ((*l1Pt).at((*indexL1SeedingL2).at(i)) + 0.01 > 7 && (*l1Quality).at((*indexL1SeedingL2).at(i)) >= 4) {
	  if ((*l2IsAssociated).at(i) && (*l2AssociationVar).at(i) > -0.2) {
	    level2GenPt_afterL2->Fill((*l2AssociatedSimMuonPt).at(i),(*l2Pt).at(i));
	  }
	} // l1 fires trigger
      } // l2 fires trigger
    } // loop over L2
  } // loop over events
  divide_histos_and_errors(level2LowPtNum,level2OverSimLowPtDenom,level2OverSimLowPtEff);
  divide_histos_and_errors(level2PtNum,level2OverSimPtDenom,level2OverSimPtEff);
  divide_histos_and_errors(level2EtaNum,level2OverSimEtaDenom,level2OverSimEtaEff);
  divide_histos_and_errors(level2PhiNum,level2OverSimPhiDenom,level2OverSimPhiEff);
  divide_TH2_histos_and_errors(level2PtEtaNum,level2OverSimPtEtaDenom,level2OverSimPtEtaEff);

  divide_histos_and_errors(level2LowPtNum,level2OverL1LowPtDenom,level2OverL1LowPtEff);
  divide_histos_and_errors(level2PtNum,level2OverL1PtDenom,level2OverL1PtEff);
  divide_histos_and_errors(level2EtaNum,level2OverL1EtaDenom,level2OverL1EtaEff);
  divide_histos_and_errors(level2PhiNum,level2OverL1PhiDenom,level2OverL1PhiEff);
  divide_TH2_histos_and_errors(level2PtEtaNum,level2OverL1PtEtaDenom,level2OverL1PtEtaEff);

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
