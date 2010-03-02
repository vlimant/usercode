#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
//#include "Math/LorentzVector.h"
#include <vector>
#include "CMS1.h"
#include <TStopwatch.h>

//int ScanTree ( TTree* tree) {
int ScanTree ( TTree* tree, char *fileName) {

  // This reads in the tree.  As you might imagine.
  Init(tree);
  TFile *histFile = new TFile(fileName,"RECREATE");
  TDirectory *histDir = histFile->mkdir("eff_hist");
  TStopwatch timer;
  timer.Start();

  // Declare the usual constants
  double pi = 3.141592653;
  double ptCut = 5.0;
  vector<int> *l1UsedBySim = new vector<int>;
  vector<int> *l2UsedBySim = new vector<int>;
  vector<int> *l3UsedBySim = new vector<int>;

  histDir->cd();
  
  //  int nPtBins = 19;
  //  Double_t pt_Edges[20] = {0, 5, 7, 9, 11, 13, 15, 18, 21, 24, 27, 30, 34, 38, 42, 50, 60, 70, 80, 100};

  //  int nPtBins = 8;
  //  Double_t pt_Edges[9] = {0, 1, 2, 3, 4, 6, 8, 10, 20};

  int nPtBins = 11;
  Double_t pt_Edges[12] = {0, 3, 5, 7, 9, 11, 13, 15, 18, 21, 24, 30};

  TH1F *l3OverXEtaEffNum = new TH1F("l3OverXEtaEffNum","l3 associated sim #eta",24,-2.4,2.4);
  TH1F *l3OverXPtEffNum = new TH1F("l3OverXPtEffNum","l3 associated sim pt",nPtBins,pt_Edges);
  TH2F *l3OverXPtEtaEffNum = new TH2F("l3OverXPtEtaEffNum","blah",nPtBins,pt_Edges,24,-2.4,2.4);

  TH1F *l2OverXEtaEffNum = new TH1F("l2OverXEtaEffNum","l2 associated sim #eta",24,-2.4,2.4);
  TH1F *l2OverXPtEffNum = new TH1F("l2OverXPtEffNum","l2 associated sim pt",nPtBins,pt_Edges);
  TH2F *l2OverXPtEtaEffNum = new TH2F("l2OverXPtEtaEffNum","blah",nPtBins,pt_Edges,24,-2.4,2.4);

  TH1F *l1OverXEtaEffNum = new TH1F("l1OverXEtaEffNum","l1 associated sim #eta",24,-2.4,2.4);
  TH1F *l1OverXPtEffNum = new TH1F("l1OverXPtEffNum","l1 associated sim pt",nPtBins,pt_Edges);
  TH2F *l1OverXPtEtaEffNum = new TH2F("l1OverXPtEtaEffNum","blah",nPtBins,pt_Edges,24,-2.4,2.4);

  TH1F *lXOverSimEtaEffDenom = new TH1F("lXOverSimEtaEffDenom","sim #eta",24,-2.4,2.4);
  TH1F *lXOverSimPtEffDenom = new TH1F("lXOverSimPtEffDenom","lX associated sim pt",nPtBins,pt_Edges);
  TH2F *lXOverSimPtEtaEffDenom = new TH2F("lXOverSimPtEtaEffDenom","blah",nPtBins,pt_Edges,24,-2.4,2.4);
  TH1F *lXOverL1EtaEffDenom = new TH1F("lXOverL1EtaEffDenom","l1 associated sim #eta",24,-2.4,2.4);
  TH1F *lXOverL1PtEffDenom = new TH1F("lXOverL1PtEffDenom","lX associated sim pt",nPtBins,pt_Edges);
  TH2F *lXOverL1PtEtaEffDenom = new TH2F("lXOverL1PtEtaEffDenom","blah",nPtBins,pt_Edges,24,-2.4,2.4);
  TH1F *lXOverL2EtaEffDenom = new TH1F("lXOverL2EtaEffDenom","l2 associated sim #eta",24,-2.4,2.4);
  TH1F *lXOverL2PtEffDenom = new TH1F("lXOverL2PtEffDenom","lX associated sim pt",nPtBins,pt_Edges);
  TH2F *lXOverL2PtEtaEffDenom = new TH2F("lXOverL2PtEtaEffDenom","blah",nPtBins,pt_Edges,24,-2.4,2.4);

  TH1F *l3OverSimEtaEff = new TH1F("l3OverSimEtaEff","l3/sim eff vs. assoc. sim #eta",24,-2.4,2.4);
  TH1F *l3OverSimPtEff = new TH1F("l3OverSimPtEff","l3/sim eff vs. assoc. sim pt",nPtBins,pt_Edges);
  TH2F *l3OverSimPtEtaEff = new TH2F("l3OverSimPtEtaEff","l3/sim eff map: p_{T} vs #eta",nPtBins,pt_Edges,24,-2.4,2.4);

  TH1F *l3OverL1EtaEff = new TH1F("l3OverL1EtaEff","l3/l1 eff vs. assoc. sim #eta",24,-2.4,2.4);
  TH1F *l3OverL1PtEff = new TH1F("l3OverL1PtEff","l3/l1 eff vs. assoc. sim pt",nPtBins,pt_Edges);
  TH2F *l3OverL1PtEtaEff = new TH2F("l3OverL1PtEtaEff","l3/l1 eff map: p_{T} vs #eta",nPtBins,pt_Edges,24,-2.4,2.4);

  TH1F *l3OverL2EtaEff = new TH1F("l3OverL2EtaEff","l3/l2 eff vs. assoc. sim #eta",24,-2.4,2.4);  
  TH1F *l3OverL2PtEff = new TH1F("l3OverL2PtEff","l3/l2 eff vs. assoc. sim pt",nPtBins,pt_Edges);
  TH2F *l3OverL2PtEtaEff = new TH2F("l3OverL2PtEtaEff","l3/l2 eff map: p_{T} vs #eta",nPtBins,pt_Edges,24,-2.4,2.4);
  
  TH1F *l2OverSimEtaEff = new TH1F("l2OverSimEtaEff","l2/sim eff vs. assoc. sim #eta",24,-2.4,2.4);
  TH1F *l2OverSimPtEff = new TH1F("l2OverSimPtEff","l2/sim eff vs. assoc. sim pt",nPtBins,pt_Edges);
  TH2F *l2OverSimPtEtaEff = new TH2F("l2OverSimPtEtaEff","l2/sim eff map: p_{T} vs #eta",nPtBins,pt_Edges,24,-2.4,2.4);

  TH1F *l2OverL1EtaEff = new TH1F("l2OverL1EtaEff","l2/l1 eff vs. assoc. sim #eta",24,-2.4,2.4);
  TH1F *l2OverL1PtEff = new TH1F("l2OverL1PtEff","l2/l1 eff vs. assoc. sim pt",nPtBins,pt_Edges);
  TH2F *l2OverL1PtEtaEff = new TH2F("l2OverL1PtEtaEff","l2/l1 eff map: p_{T} vs #eta",nPtBins,pt_Edges,24,-2.4,2.4);
 
  TH1F *l1OverSimEtaEff = new TH1F("l1OverSimEtaEff","l1/sim eff vs. assoc. sim #eta",24,-2.4,2.4);
  TH1F *l1OverSimPtEff = new TH1F("l1OverSimPtEff","l1/sim eff vs. assoc. sim pt",nPtBins,pt_Edges);
  TH2F *l1OverSimPtEtaEff = new TH2F("l1OverSimPtEtaEff","l1/sim eff map: p_{T} vs #eta",nPtBins,pt_Edges,24,-2.4,2.4);


  int nEntries = tree->GetEntries();
    
  //Event Loop
  for( int iEntry = 0; iEntry < nEntries; iEntry++) {
    if (iEntry%1000 == 0) cout << "Event " << iEntry << " time " << timer.RealTime() << " cpu " << timer.CpuTime() << endl;
    timer.Continue();
    tree->GetEntry(iEntry);
    
    double max_simPt = -999;

    try {
      //      cout << "trying passL1" << endl;
      bool passL1 = false;
      for (int iL1 = 0; iL1 < nL1; iL1++) {
	if ((*l1Pt).at(iL1) + 0.01 > 7 && (*l1Quality).at(iL1) >= 4) {
	  passL1 = true;
	}
      }

      //      cout << "trying passL2" << endl;
      bool passL2 = false;
      for (int iL2 = 0; iL2 < nL2; iL2++) {
	if (passL1) {
	  if ((*l2Pt).at(iL2) > 7 && (*l2Eta).at(iL2) > -2.5 && (*l2Eta).at(iL2) < 2.5) {
	    passL2 = true;
	  }
	}
      }
      bool passL3 = false;
      for (int iL3 = 0; iL3 < nL3; iL3++) {
	if (passL2) {
	  if ((*l3Pt).at(iL3) > 9 && (*l3Eta).at(iL3) > -2.5 && (*l3Eta).at(iL3) < 2.5) {
	    passL3 = true;
	  }
	}
      }

      for (int i = 0; i < nSimMuon; i++) {
	if ((*simMuonPt).at(i) > ptCut) {
	
	  lXOverSimEtaEffDenom->Fill((*simMuonEta).at(i));
	  lXOverSimPtEffDenom->Fill((*simMuonPt).at(i));
	  lXOverSimPtEtaEffDenom->Fill((*simMuonPt).at(i),(*simMuonEta).at(i));

	  //	cout << "trying to associate sim to L1 by DR" << endl;
	  int l1Index = STR_Association_Index((*simMuonEta).at(i),(*simMuonPhi).at(i),l1Eta,l1Phi,1.0,l1UsedBySim);
	  int l2Index = STR_Association_Index((*simMuonEta).at(i),(*simMuonPhi).at(i),l2Eta,l2Phi,0.2,l2UsedBySim);
	  int l3Index = STR_Association_Index((*simMuonEta).at(i),(*simMuonPhi).at(i),l3Eta,l3Phi,0.1,l3UsedBySim);
	  if (l1Index >= 0) {
	    l1OverXEtaEffNum->Fill((*simMuonEta).at(i));
	    l1OverXPtEffNum->Fill((*simMuonPt).at(i));
	    l1OverXPtEtaEffNum->Fill((*simMuonPt).at(i),(*simMuonEta).at(i));
	    
	    lXOverL1EtaEffDenom->Fill((*simMuonEta).at(i));
	    lXOverL1PtEffDenom->Fill((*simMuonPt).at(i));
	    lXOverL1PtEtaEffDenom->Fill((*simMuonPt).at(i),(*simMuonEta).at(i));

	    if (l2Index >= 0) {
	      l2OverXEtaEffNum->Fill((*simMuonEta).at(i));
	      l2OverXPtEffNum->Fill((*simMuonPt).at(i));
	      l2OverXPtEtaEffNum->Fill((*simMuonPt).at(i),(*simMuonEta).at(i));	      

	      lXOverL2EtaEffDenom->Fill((*simMuonEta).at(i));
	      lXOverL2PtEffDenom->Fill((*simMuonPt).at(i));
	      lXOverL2PtEtaEffDenom->Fill((*simMuonPt).at(i),(*simMuonEta).at(i));

	      if (l3Index >= 0) {
		l3OverXEtaEffNum->Fill((*simMuonEta).at(i));
		l3OverXPtEffNum->Fill((*simMuonPt).at(i));
		l3OverXPtEtaEffNum->Fill((*simMuonPt).at(i),(*simMuonEta).at(i));	      
	      } // l3 associated
	    } // l2 associated
	  } // l1 associated
	} // cut on sim pT
      } // loop over sim
      (*l1UsedBySim).clear();
      (*l2UsedBySim).clear();
      (*l3UsedBySim).clear();
    } // try statement
    catch (...) {
      cout << "Fuck it.  Pressing on" << endl;
    }
  } // event loop
  
  divide_histos_and_errors(l3OverXEtaEffNum,lXOverSimEtaEffDenom,l3OverSimEtaEff);
  divide_histos_and_errors(l3OverXEtaEffNum,lXOverL1EtaEffDenom,l3OverL1EtaEff);
  divide_histos_and_errors(l3OverXEtaEffNum,lXOverL2EtaEffDenom,l3OverL2EtaEff);
  divide_histos_and_errors(l3OverXPtEffNum,lXOverSimPtEffDenom,l3OverSimPtEff);
  divide_histos_and_errors(l3OverXPtEffNum,lXOverL1PtEffDenom,l3OverL1PtEff);
  divide_histos_and_errors(l3OverXPtEffNum,lXOverL2PtEffDenom,l3OverL2PtEff);
  divide_TH2_histos_and_errors(l3OverXPtEtaEffNum,lXOverSimPtEtaEffDenom,l3OverSimPtEtaEff);
  divide_TH2_histos_and_errors(l3OverXPtEtaEffNum,lXOverL1PtEtaEffDenom,l3OverL1PtEtaEff);
  divide_TH2_histos_and_errors(l3OverXPtEtaEffNum,lXOverL2PtEtaEffDenom,l3OverL2PtEtaEff);

  divide_histos_and_errors(l2OverXEtaEffNum,lXOverSimEtaEffDenom,l2OverSimEtaEff);
  divide_histos_and_errors(l2OverXEtaEffNum,lXOverL1EtaEffDenom,l2OverL1EtaEff);
  divide_histos_and_errors(l2OverXPtEffNum,lXOverSimPtEffDenom,l2OverSimPtEff);
  divide_histos_and_errors(l2OverXPtEffNum,lXOverL1PtEffDenom,l2OverL1PtEff);  
  divide_TH2_histos_and_errors(l2OverXPtEtaEffNum,lXOverSimPtEtaEffDenom,l2OverSimPtEtaEff);
  divide_TH2_histos_and_errors(l2OverXPtEtaEffNum,lXOverL1PtEtaEffDenom,l2OverL1PtEtaEff);

  divide_histos_and_errors(l1OverXEtaEffNum,lXOverSimEtaEffDenom,l1OverSimEtaEff);
  divide_histos_and_errors(l1OverXPtEffNum,lXOverSimPtEffDenom,l1OverSimPtEff);
  divide_TH2_histos_and_errors(l1OverXPtEtaEffNum,lXOverSimPtEtaEffDenom,l1OverSimPtEtaEff);
 
  histDir->Write("",TObject::kOverwrite);
  
  return 0;
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

int STR_Association_Index (double simEta, double simPhi, vector<double> *recoEta, vector<double> *recoPhi, double DR_CUT, vector<int> *alreadyUsed)) {
  int index = -1;
  double DeltaR = 999.9;

  //  cout << "testing this method...looping over sim" << endl;
  //  cout << "simEta size" << (*simEta).size() << endl;

  for (int i = 0; i < (*recoEta).size(); i++) {
    //    cout << "eta phi at index" << (*simEta).at(i) <<" "<< (*simPhi).at(i) << endl;
    double temp_deleta = simEta - (*recoEta).at(i);
    double temp_delphi = simPhi - (*recoPhi).at(i);      
    if (temp_delphi > 3.141592653) {
      temp_delphi = (2 * 3.141592653) - temp_delphi;
    }
    double temp_DeltaR = sqrt((temp_deleta * temp_deleta) + (temp_delphi * temp_delphi));
    if (temp_DeltaR < DeltaR  && temp_DeltaR < DR_CUT) {
      // check to make sure you're not double counting
      bool used = false;
      for (int j = 0; j < (*alreadyUsed).size(); j++) {
	if (i == (*alreadyUsed).at(j)) used = true;
      }
      if (!used) {     
	index = i;
	DeltaR = temp_DeltaR;
      }
    }
  }
  //  cout << "DeltaR = " << DeltaR << endl;
  //  cout << "returning" << endl;
  (*alreadyUsed).push_back(index);
  return index;
}
