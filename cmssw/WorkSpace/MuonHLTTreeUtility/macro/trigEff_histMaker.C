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
  double ptCut = 0.0;
  double l1PtCut = 7.0;
  double l2PtCut = 7.0;
  double l3PtCut = 9.0;
  vector<int> *simUsedByL1 = new vector<int>;
  vector<int> *simUsedByL2 = new vector<int>;
  vector<int> *simUsedByL3 = new vector<int>;

  histDir->cd();

  int nPtBins = 19;
  Double_t pt_Edges[20] = {0, 5, 7, 9, 11, 13, 15, 18, 21, 24, 27, 30, 34, 38, 42, 50, 60, 70, 80, 100};

  //  int nPtBins = 10;
  //  Double_t pt_Edges[11] = {0, 25, 50, 75, 100, 150, 200, 250, 300, 400, 500};

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

    int nL1Pass = 0;
    int nL2Pass = 0;
    int nL3Pass = 0;

    try {
      //      cout << "trying passL1" << endl;
      for (int iL1 = 0; iL1 < nL1; iL1++) {
	if ((*l1Pt).at(iL1) + 0.01 > l1PtCut && (*l1Quality).at(iL1) >= 4) {
	  //	  if (iEntry ==290) cout << "we should have a pass here" << endl;
	  int simIndexL1 = RTS_Association_Index((*l1Eta).at(iL1),(*l1Phi).at(iL1),simMuonEta,simMuonPhi,simMuonPt,10.0,ptCut,simUsedByL1);
	  //	  if (iEntry ==290) cout << "simIndexL1 = " << simIndexL1 << endl;
	  //	  if (simIndexL1 < 0) cout << "we lose an L1 from association" << endl;
	  if (simIndexL1 >= 0 && (*simMuonPt).at(simIndexL1) > ptCut) {
	    nL1Pass++;
	    //	    cout <<  "pt, eta, phi, quality,: " << (*l1Pt).at(iL1) << " " << (*l1Eta).at(iL1) <<" " << (*l1Phi).at(iL1) << " " << (*l1Quality).at(iL1) << endl;
	    l1OverXEtaEffNum->Fill((*simMuonEta).at(simIndexL1));
	    l1OverXPtEffNum->Fill((*simMuonPt).at(simIndexL1));
	    l1OverXPtEtaEffNum->Fill((*simMuonPt).at(simIndexL1),(*simMuonEta).at(simIndexL1));
	    
	    lXOverL1EtaEffDenom->Fill((*simMuonEta).at(simIndexL1));	  
	    lXOverL1PtEffDenom->Fill((*simMuonPt).at(simIndexL1));	  
	    lXOverL1PtEtaEffDenom->Fill((*simMuonPt).at(simIndexL1),(*simMuonEta).at(simIndexL1));
	  }
	}
      }

      (*simUsedByL1).clear();

      //      cout << "loop over l2" << endl;

      for (int iL2 = 0; iL2 < nL2; iL2++) {
	if ((*l2Pt).at(iL2) > l2PtCut && (*l2Eta).at(iL2) > -2.5 && (*l2Eta).at(iL2) < 2.5) {
	  //  	  cout << "trying to find l1 seed" << endl;	  
	  int indexL1 = (*indexL1SeedingL2).at(iL2);
	  //	  if (iEntry == 290) cout << (*l1Pt).at(indexL1) << endl;
	  if ((*l1Pt).at(indexL1) + 0.01 > l1PtCut && (*l1Quality).at(indexL1) >= 4) {	    
	    int simIndexL2 = RTS_Association_Index((*l2Eta).at(iL2),(*l2Phi).at(iL2),simMuonEta,simMuonPhi,simMuonPt,10.0,ptCut,simUsedByL2);
	    if (simIndexL2 >= 0 && (*simMuonPt).at(simIndexL2) > ptCut) {	
	      nL2Pass++;
	      l2OverXEtaEffNum->Fill((*simMuonEta).at(simIndexL2));
	      l2OverXPtEffNum->Fill((*simMuonPt).at(simIndexL2));
	      l2OverXPtEtaEffNum->Fill((*simMuonPt).at(simIndexL2),(*simMuonEta).at(simIndexL2));

	      lXOverL2EtaEffDenom->Fill((*simMuonEta).at(simIndexL2));
	      lXOverL2PtEffDenom->Fill((*simMuonPt).at(simIndexL2));	  
	      lXOverL2PtEtaEffDenom->Fill((*simMuonPt).at(simIndexL2),(*simMuonEta).at(simIndexL2));
	    }
	  }
	}
      }

      (*simUsedByL2).clear();

      //            cout << "loop over l3" << endl;
      for (int iL3 = 0; iL3 < nL3; iL3++) {
	if ((*l3Pt).at(iL3) > l3PtCut && (*l3Eta).at(iL3) > -2.5 && (*l3Eta).at(iL3)< 2.5) {
	  //	  cout << "trying to get l2 seed for l3" << endl;
	  int indexL2 = (*indexL2SeedingL3).at(iL3);
	  //	  cout << "got it" << endl;
	  if ((*l2Pt).at(indexL2) > l2PtCut && (*l2Eta).at(indexL2) > -2.5 && (*l2Eta).at(indexL2) < 2.5) {
	    //	    cout << "trying to get l1 seed for l2 seeding l3" << endl;
	    int indexL1 = (*indexL1SeedingL2).at(indexL2);
	    //	    cout << "got it" << endl;
	    if ((*l1Pt).at(indexL1) + 0.01 > l1PtCut && (*l1Quality).at(indexL1) >= 4) {
	      int simIndexL3 = RTS_Association_Index((*l3Eta).at(iL3),(*l3Phi).at(iL3),simMuonEta,simMuonPhi,simMuonPt,10.0,ptCut,simUsedByL3);
	      if (simIndexL3 >= 0 && (*simMuonPt).at(simIndexL3) > ptCut) {
		nL3Pass++;
		l3OverXEtaEffNum->Fill((*simMuonEta).at(simIndexL3));
		l3OverXPtEffNum->Fill((*simMuonPt).at(simIndexL3));
		l3OverXPtEtaEffNum->Fill((*simMuonPt).at(simIndexL3),(*simMuonEta).at(simIndexL3));
	      }
	    }
	  }
	}
      }

      (*simUsedByL3).clear();

      for (int i = 0; i < nSimMuon; i++) {
	if ((*simMuonPt).at(i) > ptCut && (*simMuonEta).at(i) < 2.1 && (*simMuonEta).at(i) > -2.1) {
	  lXOverSimEtaEffDenom->Fill((*simMuonEta).at(i));
	  lXOverSimPtEffDenom->Fill((*simMuonPt).at(i));	  
	  lXOverSimPtEtaEffDenom->Fill((*simMuonPt).at(i),(*simMuonEta).at(i));	
	}
      }
    }

    catch (...) {
      cout << "Fuck it.  Pressing on" << endl;
    }
  }
  
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

int STR_Association_Index (double simEta, double simPhi, vector<double> *recoEta, vector<double> *recoPhi, double DR_CUT, vector<int> *alreadyUsed) {
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

int RTS_Association_Index (double recoEta, double recoPhi, vector<double> *simEta, vector<double> *simPhi, vector<double> *simPt, double DR_CUT, double PT_CUT, vector<int> *alreadyUsed) {
  int index = -1;
  double DeltaR = 999.9;

  for (int i = 0; i < (*simEta).size(); i++) {
    if ((*simPt).at(i) > PT_CUT && (*simEta).at(i) > -2.1 &&
	(*simEta).at(i) < 2.1) {
      double temp_deleta = recoEta - (*simEta).at(i);
      double temp_delphi = fabs(recoPhi - (*simPhi).at(i));
      if (temp_delphi > 3.141592653) {
	temp_delphi = (2 * 3.141592653) - temp_delphi;
      }
      double temp_DeltaR = sqrt((temp_deleta * temp_deleta) + (temp_delphi * temp_delphi));
      if (temp_DeltaR < DeltaR  && temp_DeltaR < DR_CUT) {
	// check to make sure you're not double counting
	bool used = false;
	for (int j = 0; j < (*alreadyUsed).size(); j++) {
	  if (i == (*alreadyUsed).at(j)) {
	    used = true;
	  }
	}
	if (!used) {
	  index = i;
	  DeltaR = temp_DeltaR;
	}
      }
    }
  }
  //  cout << "DeltaR = " << DeltaR << endl;
  //  cout << "returning " << index << endl;
  (*alreadyUsed).push_back(index);
  return index;
}

