#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
//#include "Math/LorentzVector.h"
#include <vector>
#include "CMS1.h"

//int ScanTree ( TTree* tree) {
int ScanTree ( TTree* tree, char *fileName) {

  // This reads in the tree.  As you might imagine.
  Init(tree);
  TFile *histFile = new TFile(fileName,"RECREATE");
  TDirectory *histDir = histFile->mkdir("rate_hist");

  histDir->cd();

  // Time for some histograms and stacks

  double ptUpperLimit = 15;
  TH1F *l2PtRate_motherBin_1 = new TH1F("l2PtRate_motherBin_1"," #pi^{+/-}",60,0,ptUpperLimit);
  TH1F *l2PtRate_motherBin_2 = new TH1F("l2PtRate_motherBin_2"," K",60,0,ptUpperLimit);
  TH1F *l2PtRate_motherBin_3 = new TH1F("l2PtRate_motherBin_3"," D",60,0,ptUpperLimit);
  TH1F *l2PtRate_motherBin_4 = new TH1F("l2PtRate_motherBin_4"," B",60,0,ptUpperLimit);
  TH1F *l2PtRate_motherBin_5 = new TH1F("l2PtRate_motherBin_5"," #Lambda_{b}",60,0,ptUpperLimit);
  TH1F *l2PtRate_motherBin_6 = new TH1F("l2PtRate_motherBin_6"," J/#Psi",60,0,ptUpperLimit);
  TH1F *l2PtRate_motherBin_10 = new TH1F("l2PtRate_motherBin_10"," #tau^{+/-}",60,0,ptUpperLimit);
  TH1F *l2PtRate_motherBin_12 = new TH1F("l2PtRate_motherBin_12"," other",60,0,ptUpperLimit);
  TH1F *l2PtRate_motherBin_13 = new TH1F("l2PtRate_motherBin_13"," non-associated",60,0,ptUpperLimit);

  TH1F *l2PtRate_motherBin_1_cd = new TH1F("l2PtRate_motherBin_1_cd"," #pi^{+/-}",60,0,ptUpperLimit);
  TH1F *l2PtRate_motherBin_2_cd = new TH1F("l2PtRate_motherBin_2_cd"," K",60,0,ptUpperLimit);
  TH1F *l2PtRate_motherBin_3_cd = new TH1F("l2PtRate_motherBin_3_cd"," D",60,0,ptUpperLimit);
  TH1F *l2PtRate_motherBin_4_cd = new TH1F("l2PtRate_motherBin_4_cd"," B",60,0,ptUpperLimit);
  TH1F *l2PtRate_motherBin_5_cd = new TH1F("l2PtRate_motherBin_5_cd"," #Lambda_{b}",60,0,ptUpperLimit);
  TH1F *l2PtRate_motherBin_6_cd = new TH1F("l2PtRate_motherBin_6_cd"," J/#Psi",60,0,ptUpperLimit);
  TH1F *l2PtRate_motherBin_10_cd = new TH1F("l2PtRate_motherBin_10_cd"," #tau^{+/-}",60,0,ptUpperLimit);
  TH1F *l2PtRate_motherBin_12_cd = new TH1F("l2PtRate_motherBin_12_cd"," other",60,0,ptUpperLimit);
  TH1F *l2PtRate_motherBin_13_cd = new TH1F("l2PtRate_motherBin_13_cd"," non-associated",60,0,ptUpperLimit);

  TH1F *l2IsoPtRate_motherBin_1 = new TH1F("l2IsoPtRate_motherBin_1"," #pi^{+/-}",60,0,ptUpperLimit);
  TH1F *l2IsoPtRate_motherBin_2 = new TH1F("l2IsoPtRate_motherBin_2"," K",60,0,ptUpperLimit);
  TH1F *l2IsoPtRate_motherBin_3 = new TH1F("l2IsoPtRate_motherBin_3"," D",60,0,ptUpperLimit);
  TH1F *l2IsoPtRate_motherBin_4 = new TH1F("l2IsoPtRate_motherBin_4"," B",60,0,ptUpperLimit);
  TH1F *l2IsoPtRate_motherBin_5 = new TH1F("l2IsoPtRate_motherBin_5"," #Lambda_{b}",60,0,ptUpperLimit);
  TH1F *l2IsoPtRate_motherBin_6 = new TH1F("l2IsoPtRate_motherBin_6"," J/#Psi",60,0,ptUpperLimit);
  TH1F *l2IsoPtRate_motherBin_10 = new TH1F("l2IsoPtRate_motherBin_10"," #tau^{+/-}",60,0,ptUpperLimit);
  TH1F *l2IsoPtRate_motherBin_12 = new TH1F("l2IsoPtRate_motherBin_12"," other",60,0,ptUpperLimit);
  TH1F *l2IsoPtRate_motherBin_13 = new TH1F("l2IsoPtRate_motherBin_13"," non-associated",60,0,ptUpperLimit);

  TH1F *l2IsoPtRate_motherBin_1_cd = new TH1F("l2IsoPtRate_motherBin_1_cd"," #pi^{+/-}",60,0,ptUpperLimit);
  TH1F *l2IsoPtRate_motherBin_2_cd = new TH1F("l2IsoPtRate_motherBin_2_cd"," K",60,0,ptUpperLimit);
  TH1F *l2IsoPtRate_motherBin_3_cd = new TH1F("l2IsoPtRate_motherBin_3_cd"," D",60,0,ptUpperLimit);
  TH1F *l2IsoPtRate_motherBin_4_cd = new TH1F("l2IsoPtRate_motherBin_4_cd"," B",60,0,ptUpperLimit);
  TH1F *l2IsoPtRate_motherBin_5_cd = new TH1F("l2IsoPtRate_motherBin_5_cd"," #Lambda_{b}",60,0,ptUpperLimit);
  TH1F *l2IsoPtRate_motherBin_6_cd = new TH1F("l2IsoPtRate_motherBin_6_cd"," J/#Psi",60,0,ptUpperLimit);
  TH1F *l2IsoPtRate_motherBin_10_cd = new TH1F("l2IsoPtRate_motherBin_10_cd"," #tau^{+/-}",60,0,ptUpperLimit);
  TH1F *l2IsoPtRate_motherBin_12_cd = new TH1F("l2IsoPtRate_motherBin_12_cd"," other",60,0,ptUpperLimit);
  TH1F *l2IsoPtRate_motherBin_13_cd = new TH1F("l2IsoPtRate_motherBin_13_cd"," non-associated",60,0,ptUpperLimit);

  THStack *l2PtRate = new THStack("l2PtRate","L2 rate as f(p_{T,L2})");
  THStack *l2PtRate_cd = new THStack("l2PtRate_cd","L2 rate as f(p_{T,L2})");
  THStack *l2IsoPtRate = new THStack("l2IsoPtRate","L2Iso rate as f(p_{T,L2})");
  THStack *l2IsoPtRate_cd = new THStack("l2IsoPtRate_cd","L2Iso rate as f(p_{T,L2})");

  TH1F *l3PtRate_motherBin_1 = new TH1F("l3PtRate_motherBin_1"," #pi^{+/-}",60,0,ptUpperLimit);
  TH1F *l3PtRate_motherBin_2 = new TH1F("l3PtRate_motherBin_2"," K",60,0,ptUpperLimit);
  TH1F *l3PtRate_motherBin_3 = new TH1F("l3PtRate_motherBin_3"," D",60,0,ptUpperLimit);
  TH1F *l3PtRate_motherBin_4 = new TH1F("l3PtRate_motherBin_4"," B",60,0,ptUpperLimit);
  TH1F *l3PtRate_motherBin_5 = new TH1F("l3PtRate_motherBin_5"," #Lambda_{b}",60,0,ptUpperLimit);
  TH1F *l3PtRate_motherBin_6 = new TH1F("l3PtRate_motherBin_6"," J/#Psi",60,0,ptUpperLimit);
  TH1F *l3PtRate_motherBin_10 = new TH1F("l3PtRate_motherBin_10"," #tau^{+/-}",60,0,ptUpperLimit);
  TH1F *l3PtRate_motherBin_12 = new TH1F("l3PtRate_motherBin_12"," other",60,0,ptUpperLimit);
  TH1F *l3PtRate_motherBin_13 = new TH1F("l3PtRate_motherBin_13"," non-associated",60,0,ptUpperLimit);

  TH1F *l3PtRate_motherBin_1_cd = new TH1F("l3PtRate_motherBin_1_cd"," #pi^{+/-}",60,0,ptUpperLimit);
  TH1F *l3PtRate_motherBin_2_cd = new TH1F("l3PtRate_motherBin_2_cd"," K",60,0,ptUpperLimit);
  TH1F *l3PtRate_motherBin_3_cd = new TH1F("l3PtRate_motherBin_3_cd"," D",60,0,ptUpperLimit);
  TH1F *l3PtRate_motherBin_4_cd = new TH1F("l3PtRate_motherBin_4_cd"," B",60,0,ptUpperLimit);
  TH1F *l3PtRate_motherBin_5_cd = new TH1F("l3PtRate_motherBin_5_cd"," #Lambda_{b}",60,0,ptUpperLimit);
  TH1F *l3PtRate_motherBin_6_cd = new TH1F("l3PtRate_motherBin_6_cd"," J/#Psi",60,0,ptUpperLimit);
  TH1F *l3PtRate_motherBin_10_cd = new TH1F("l3PtRate_motherBin_10_cd"," #tau^{+/-}",60,0,ptUpperLimit);
  TH1F *l3PtRate_motherBin_12_cd = new TH1F("l3PtRate_motherBin_12_cd"," other",60,0,ptUpperLimit);
  TH1F *l3PtRate_motherBin_13_cd = new TH1F("l3PtRate_motherBin_13_cd"," non-associated",60,0,ptUpperLimit);

  TH1F *l3PreIsoPtRate_motherBin_1 = new TH1F("l3PreIsoPtRate_motherBin_1"," #pi^{+/-}",60,0,ptUpperLimit);
  TH1F *l3PreIsoPtRate_motherBin_2 = new TH1F("l3PreIsoPtRate_motherBin_2"," K",60,0,ptUpperLimit);
  TH1F *l3PreIsoPtRate_motherBin_3 = new TH1F("l3PreIsoPtRate_motherBin_3"," D",60,0,ptUpperLimit);
  TH1F *l3PreIsoPtRate_motherBin_4 = new TH1F("l3PreIsoPtRate_motherBin_4"," B",60,0,ptUpperLimit);
  TH1F *l3PreIsoPtRate_motherBin_5 = new TH1F("l3PreIsoPtRate_motherBin_5"," #Lambda_{b}",60,0,ptUpperLimit);
  TH1F *l3PreIsoPtRate_motherBin_6 = new TH1F("l3PreIsoPtRate_motherBin_6"," J/#Psi",60,0,ptUpperLimit);
  TH1F *l3PreIsoPtRate_motherBin_10 = new TH1F("l3PreIsoPtRate_motherBin_10"," #tau^{+/-}",60,0,ptUpperLimit);
  TH1F *l3PreIsoPtRate_motherBin_12 = new TH1F("l3PreIsoPtRate_motherBin_12"," other",60,0,ptUpperLimit);
  TH1F *l3PreIsoPtRate_motherBin_13 = new TH1F("l3PreIsoPtRate_motherBin_13"," non-associated",60,0,ptUpperLimit);

  TH1F *l3PreIsoPtRate_motherBin_1_cd = new TH1F("l3PreIsoPtRate_motherBin_1_cd"," #pi^{+/-}",60,0,ptUpperLimit);
  TH1F *l3PreIsoPtRate_motherBin_2_cd = new TH1F("l3PreIsoPtRate_motherBin_2_cd"," K",60,0,ptUpperLimit);
  TH1F *l3PreIsoPtRate_motherBin_3_cd = new TH1F("l3PreIsoPtRate_motherBin_3_cd"," D",60,0,ptUpperLimit);
  TH1F *l3PreIsoPtRate_motherBin_4_cd = new TH1F("l3PreIsoPtRate_motherBin_4_cd"," B",60,0,ptUpperLimit);
  TH1F *l3PreIsoPtRate_motherBin_5_cd = new TH1F("l3PreIsoPtRate_motherBin_5_cd"," #Lambda_{b}",60,0,ptUpperLimit);
  TH1F *l3PreIsoPtRate_motherBin_6_cd = new TH1F("l3PreIsoPtRate_motherBin_6_cd"," J/#Psi",60,0,ptUpperLimit);
  TH1F *l3PreIsoPtRate_motherBin_10_cd = new TH1F("l3PreIsoPtRate_motherBin_10_cd"," #tau^{+/-}",60,0,ptUpperLimit);
  TH1F *l3PreIsoPtRate_motherBin_12_cd = new TH1F("l3PreIsoPtRate_motherBin_12_cd"," other",60,0,ptUpperLimit);
  TH1F *l3PreIsoPtRate_motherBin_13_cd = new TH1F("l3PreIsoPtRate_motherBin_13_cd"," non-associated",60,0,ptUpperLimit);

  TH1F *l3IsoPtRate_motherBin_1 = new TH1F("l3IsoPtRate_motherBin_1"," #pi^{+/-}",60,0,ptUpperLimit);
  TH1F *l3IsoPtRate_motherBin_2 = new TH1F("l3IsoPtRate_motherBin_2"," K",60,0,ptUpperLimit);
  TH1F *l3IsoPtRate_motherBin_3 = new TH1F("l3IsoPtRate_motherBin_3"," D",60,0,ptUpperLimit);
  TH1F *l3IsoPtRate_motherBin_4 = new TH1F("l3IsoPtRate_motherBin_4"," B",60,0,ptUpperLimit);
  TH1F *l3IsoPtRate_motherBin_5 = new TH1F("l3IsoPtRate_motherBin_5"," #Lambda_{b}",60,0,ptUpperLimit);
  TH1F *l3IsoPtRate_motherBin_6 = new TH1F("l3IsoPtRate_motherBin_6"," J/#Psi",60,0,ptUpperLimit);
  TH1F *l3IsoPtRate_motherBin_10 = new TH1F("l3IsoPtRate_motherBin_10"," #tau^{+/-}",60,0,ptUpperLimit);
  TH1F *l3IsoPtRate_motherBin_12 = new TH1F("l3IsoPtRate_motherBin_12"," other",60,0,ptUpperLimit);
  TH1F *l3IsoPtRate_motherBin_13 = new TH1F("l3IsoPtRate_motherBin_13"," non-associated",60,0,ptUpperLimit);

  TH1F *l3IsoPtRate_motherBin_1_cd = new TH1F("l3IsoPtRate_motherBin_1_cd"," #pi^{+/-}",60,0,ptUpperLimit);
  TH1F *l3IsoPtRate_motherBin_2_cd = new TH1F("l3IsoPtRate_motherBin_2_cd"," K",60,0,ptUpperLimit);
  TH1F *l3IsoPtRate_motherBin_3_cd = new TH1F("l3IsoPtRate_motherBin_3_cd"," D",60,0,ptUpperLimit);
  TH1F *l3IsoPtRate_motherBin_4_cd = new TH1F("l3IsoPtRate_motherBin_4_cd"," B",60,0,ptUpperLimit);
  TH1F *l3IsoPtRate_motherBin_5_cd = new TH1F("l3IsoPtRate_motherBin_5_cd"," #Lambda_{b}",60,0,ptUpperLimit);
  TH1F *l3IsoPtRate_motherBin_6_cd = new TH1F("l3IsoPtRate_motherBin_6_cd"," J/#Psi",60,0,ptUpperLimit);
  TH1F *l3IsoPtRate_motherBin_10_cd = new TH1F("l3IsoPtRate_motherBin_10_cd"," #tau^{+/-}",60,0,ptUpperLimit);
  TH1F *l3IsoPtRate_motherBin_12_cd = new TH1F("l3IsoPtRate_motherBin_12_cd"," other",60,0,ptUpperLimit);
  TH1F *l3IsoPtRate_motherBin_13_cd = new TH1F("l3IsoPtRate_motherBin_13_cd"," non-associated",60,0,ptUpperLimit);

  THStack *l3PtRate = new THStack("l3PtRate","L3 rate as f(p_{T,L3})");
  THStack *l3PtRate_cd = new THStack("l3PtRate_cd","L3 rate as f(p_{T,L3})");
  THStack *l3PreIsoPtRate = new THStack("l3PreIsoPtRate","L3PreIso rate as f(p_{T,L3})");
  THStack *l3PreIsoPtRate_cd = new THStack("l3PreIsoPtRate_cd","L3PreIso rate as f(p_{T,L3})");
  THStack *l3IsoPtRate = new THStack("l3IsoPtRate","L3Iso rate as f(p_{T,L3})");
  THStack *l3IsoPtRate_cd = new THStack("l3IsoPtRate_cd","L3Iso rate as f(p_{T,L3})");

  // There are many comments to be made about ROOT.  Most of them Rated-R or Rated-X.
  // Here, I'll limit myself to saying this is necessary for the TStacks to go into the TDirectory.
  gDirectory->Append(l2PtRate);
  gDirectory->Append(l3PtRate);
  gDirectory->Append(l2IsoPtRate);
  gDirectory->Append(l3PreIsoPtRate);
  gDirectory->Append(l3IsoPtRate);
  gDirectory->Append(l2PtRate_cd);
  gDirectory->Append(l3PtRate_cd);
  gDirectory->Append(l2IsoPtRate_cd);
  gDirectory->Append(l3PreIsoPtRate_cd);
  gDirectory->Append(l3IsoPtRate_cd);

  int l2OverFlow[13] = {0,0,0,0,0,0,0,0,0,0,0,0,0};
  int l3OverFlow[13] = {0,0,0,0,0,0,0,0,0,0,0,0,0};
  int l2IsoOverFlow[13] = {0,0,0,0,0,0,0,0,0,0,0,0,0};
  int l3PreIsoOverFlow[13] = {0,0,0,0,0,0,0,0,0,0,0,0,0};
  int l3IsoOverFlow[13] = {0,0,0,0,0,0,0,0,0,0,0,0,0};

  int nEntries = tree->GetEntries();
    
  //Event Loop
  for( int iEntry = 0; iEntry < nEntries; iEntry++) {
    if (iEntry%1000 == 0) cout << "Event " << iEntry << endl;
    tree->GetEntry(iEntry);
    
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
        //      if ((*l2Pt).at(i) > 7 && abs((*l2Eta).at(i)) < 2.5) {
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
    bool passL3iso = false;
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
	  if ((*l3TrackIsoDeposit).at(i) < etaTrackIsoDepositCut((*l2Eta).at(i))) {
	    passL3iso = true;
	    continue;
	  }
        }
      }
    }

    if (passL2) {
      double pt_L2 = findMaxPt(l2Pt);
      int l2_index = findIndexOfMaxPt(l2Pt);
      if ((*l2IsAssociated).at(l2_index) == 1 && (*l2AssociationVar).at(l2_index) > -0.2) {
	switch ( (*l2MotherBinNumber).at(l2_index) ) {
	case 1 :
	  l2PtRate_motherBin_1->Fill(pt_L2);
	  if (pt_L2 > ptUpperLimit) l2OverFlow[0] ++;
	  break;
	case 2 :
	  l2PtRate_motherBin_2->Fill(pt_L2);
	  if (pt_L2 > ptUpperLimit) l2OverFlow[1] ++;
	  break;
	case 3 :
	  l2PtRate_motherBin_3->Fill(pt_L2);
	  if (pt_L2 > ptUpperLimit) l2OverFlow[2] ++;
	  break;
	case 4 :
	  l2PtRate_motherBin_4->Fill(pt_L2);
	  if (pt_L2 > ptUpperLimit) l2OverFlow[3] ++;
	  break;
	case 5 :
	  l2PtRate_motherBin_5->Fill(pt_L2);
	  if (pt_L2 > ptUpperLimit) l2OverFlow[4] ++;
	  break;
	case 6 :
	  l2PtRate_motherBin_6->Fill(pt_L2);
	  if (pt_L2 > ptUpperLimit) l2OverFlow[5] ++;
	  break;
	case 10 :
	  l2PtRate_motherBin_10->Fill(pt_L2);
	  if (pt_L2 > ptUpperLimit) l2OverFlow[9] ++;
	  break;
	case 12 :
	  l2PtRate_motherBin_12->Fill(pt_L2);
	  if (pt_L2 > ptUpperLimit) l2OverFlow[11] ++;
	  break;
	default :
	  //      cout << "this is interesting" << endl;
	  l2PtRate_motherBin_12->Fill(pt_L2);
	  if (pt_L2 > ptUpperLimit) l2OverFlow[11] ++;
	}
      }
      else {
	l2PtRate_motherBin_13->Fill(pt_L2);
	if (pt_L2 > ptUpperLimit) l2OverFlow[12] ++;
      }
    }

    if (passL2Iso) {
      double pt_L2 = findMaxPt(l2Pt);
      int l2_index = findIndexOfMaxPt(l2Pt);
      if ((*l2IsAssociated).at(l2_index) == 1 && (*l2AssociationVar).at(l2_index) > -0.2) {
	switch ( (*l2MotherBinNumber).at(l2_index) ) {
	case 1 :
	  l2IsoPtRate_motherBin_1->Fill(pt_L2);
	  if (pt_L2 > ptUpperLimit) l2IsoOverFlow[0] ++;
	  break;
	case 2 :
	  l2IsoPtRate_motherBin_2->Fill(pt_L2);
	  if (pt_L2 > ptUpperLimit) l2IsoOverFlow[1] ++;
	  break;
	case 3 :
	  l2IsoPtRate_motherBin_3->Fill(pt_L2);
	  if (pt_L2 > ptUpperLimit) l2IsoOverFlow[2] ++;
	  break;
	case 4 :
	  l2IsoPtRate_motherBin_4->Fill(pt_L2);
	  if (pt_L2 > ptUpperLimit) l2IsoOverFlow[3] ++;
	  break;
	case 5 :
	  l2IsoPtRate_motherBin_5->Fill(pt_L2);
	  if (pt_L2 > ptUpperLimit) l2IsoOverFlow[4] ++;
	  break;
	case 6 :
	  l2IsoPtRate_motherBin_6->Fill(pt_L2);
	  if (pt_L2 > ptUpperLimit) l2IsoOverFlow[5] ++;
	  break;
	case 10 :
	  l2IsoPtRate_motherBin_10->Fill(pt_L2);
	  if (pt_L2 > ptUpperLimit) l2IsoOverFlow[9] ++;
	  break;
	case 12 :
	  l2IsoPtRate_motherBin_12->Fill(pt_L2);
	  if (pt_L2 > ptUpperLimit) l2IsoOverFlow[11] ++;
	  break;
	default :
	  //           cout << "this is interesting" << endl;
	  l2IsoPtRate_motherBin_12->Fill(pt_L2);
	  if (pt_L2 > ptUpperLimit) l2IsoOverFlow[11] ++;
	}
      }
      else {
	l2IsoPtRate_motherBin_13->Fill(pt_L2);
	if (pt_L2 > ptUpperLimit) l2IsoOverFlow[12] ++;
      }
    }

    if (passL3) {
      double pt_L3 = findMaxPt(l3Pt,l3Eta,l3D0);
      int l3_index = findIndexOfMaxPt(l3Pt,l3Eta,l3D0);
      if ((*l3IsAssociated).at(l3_index) == 1 && (*l3AssociationVar).at(l3_index) > -0.1) {
	switch ( (*l3MotherBinNumber).at(l3_index) ) {
	case 1 :
	  l3PtRate_motherBin_1->Fill(pt_L3);
	  if (pt_L3 > ptUpperLimit) l3OverFlow[0] ++;
	  break;
	case 2 :
	  l3PtRate_motherBin_2->Fill(pt_L3);
	  if (pt_L3 > ptUpperLimit) l3OverFlow[1] ++;
	  break;
	case 3 :
	  l3PtRate_motherBin_3->Fill(pt_L3);
	  if (pt_L3 > ptUpperLimit) l3OverFlow[2] ++;
	  break;
	case 4 :
	  l3PtRate_motherBin_4->Fill(pt_L3);
	  if (pt_L3 > ptUpperLimit) l3OverFlow[3] ++;
	  break;
	case 5 :
	  l3PtRate_motherBin_5->Fill(pt_L3);
	  if (pt_L3 > ptUpperLimit) l3OverFlow[4] ++;
	  break;
	case 6 :
	  l3PtRate_motherBin_6->Fill(pt_L3);
	  if (pt_L3 > ptUpperLimit) l3OverFlow[5] ++;
	  break;
	case 10 :
	  l3PtRate_motherBin_10->Fill(pt_L3);
	  if (pt_L3 > ptUpperLimit) l3OverFlow[9] ++;
	  break;
	case 12 :
	  l3PtRate_motherBin_12->Fill(pt_L3);
	  if (pt_L3 > ptUpperLimit) l3OverFlow[11] ++;
	  break;
	default :
	  //           cout << "oh, cock, this is interesting" << endl;
	  l3PtRate_motherBin_12->Fill(pt_L3);
	  if (pt_L3 > ptUpperLimit) l3OverFlow[11] ++;
	}
      }
      else {
	l3PtRate_motherBin_13->Fill(pt_L3);
	if (pt_L3 > ptUpperLimit) l3OverFlow[12] ++;
      }
    }
    
    if (passL3pre) {
      double pt_L3 = findMaxPt(l3Pt,l3Eta,l3D0);
      int l3_index = findIndexOfMaxPt(l3Pt,l3Eta,l3D0);
      if ((*l3IsAssociated).at(l3_index) == 1 && (*l3AssociationVar).at(l3_index) > -0.1) {
	switch ( (*l3MotherBinNumber).at(l3_index) ) {
	case 1 :
	  l3PreIsoPtRate_motherBin_1->Fill(pt_L3);
	  if (pt_L3 > ptUpperLimit) l3PreIsoOverFlow[0] ++;
	  break;
	case 2 :
	  l3PreIsoPtRate_motherBin_2->Fill(pt_L3);
	  if (pt_L3 > ptUpperLimit) l3PreIsoOverFlow[1] ++;
	  break;
	case 3 :
	  l3PreIsoPtRate_motherBin_3->Fill(pt_L3);
	  if (pt_L3 > ptUpperLimit) l3PreIsoOverFlow[2] ++;
	  break;
	case 4 :
	  l3PreIsoPtRate_motherBin_4->Fill(pt_L3);
	  if (pt_L3 > ptUpperLimit) l3PreIsoOverFlow[3] ++;
	  break;
	case 5 :
	  l3PreIsoPtRate_motherBin_5->Fill(pt_L3);
	  if (pt_L3 > ptUpperLimit) l3PreIsoOverFlow[4] ++;
	  break;
	case 6 :
	  l3PreIsoPtRate_motherBin_6->Fill(pt_L3);
	  if (pt_L3 > ptUpperLimit) l3PreIsoOverFlow[5] ++;
	  break;
	case 10 :
	  l3PreIsoPtRate_motherBin_10->Fill(pt_L3);
	  if (pt_L3 > ptUpperLimit) l3PreIsoOverFlow[9] ++;
	  break;
	case 12 :
	  l3PreIsoPtRate_motherBin_12->Fill(pt_L3);
	  if (pt_L3 > ptUpperLimit) l3PreIsoOverFlow[11] ++;
	  break;
	default :
	  //           cout << "oh, cock, this is interesting" << endl;
	  l3PreIsoPtRate_motherBin_12->Fill(pt_L3);
	  if (pt_L3 > ptUpperLimit) l3PreIsoOverFlow[11] ++;
	}
      }
      else {
	l3PreIsoPtRate_motherBin_13->Fill(pt_L3);
	if (pt_L3 > ptUpperLimit) l3PreIsoOverFlow[12] ++;
      }
    }

    if (passL3iso) {
      double pt_L3 = findMaxPt(l3Pt,l3Eta,l3D0);
      int l3_index = findIndexOfMaxPt(l3Pt,l3Eta,l3D0);
      if ((*l3IsAssociated).at(l3_index) == 1 && (*l3AssociationVar).at(l3_index) > -0.1) {
	switch ( (*l3MotherBinNumber).at(l3_index) ) {
	case 1 :
	  l3IsoPtRate_motherBin_1->Fill(pt_L3);
	  if (pt_L3 > ptUpperLimit) l3IsoOverFlow[0] ++;
	  break;
	case 2 :
	  l3IsoPtRate_motherBin_2->Fill(pt_L3);
	  if (pt_L3 > ptUpperLimit) l3IsoOverFlow[1] ++;
	  break;
	case 3 :
	  l3IsoPtRate_motherBin_3->Fill(pt_L3);
	  if (pt_L3 > ptUpperLimit) l3IsoOverFlow[2] ++;
	  break;
	case 4 :
	  l3IsoPtRate_motherBin_4->Fill(pt_L3);
	  if (pt_L3 > ptUpperLimit) l3IsoOverFlow[3] ++;
	  break;
	case 5 :
	  l3IsoPtRate_motherBin_5->Fill(pt_L3);
	  if (pt_L3 > ptUpperLimit) l3IsoOverFlow[4] ++;
	  break;
	case 6 :
	  l3IsoPtRate_motherBin_6->Fill(pt_L3);
	  if (pt_L3 > ptUpperLimit) l3IsoOverFlow[5] ++;
	  break;
	case 10 :
	  l3IsoPtRate_motherBin_10->Fill(pt_L3);
	  if (pt_L3 > ptUpperLimit) l3IsoOverFlow[9] ++;
	  break;
	case 12 :
	  l3IsoPtRate_motherBin_12->Fill(pt_L3);
	  if (pt_L3 > ptUpperLimit) l3IsoOverFlow[11] ++;
	  break;
	default :
	  //           cout << "oh, cock, this is interesting" << endl;
	  l3IsoPtRate_motherBin_12->Fill(pt_L3);
	  if (pt_L3 > ptUpperLimit) l3IsoOverFlow[11] ++;
	}
      }
      else {
	l3IsoPtRate_motherBin_13->Fill(pt_L3);
	if (pt_L3 > ptUpperLimit) l3IsoOverFlow[12] ++;
      }
    }
  }

  fill_overflow(l2PtRate_motherBin_1,l2OverFlow[0]);
  fill_overflow(l2PtRate_motherBin_2,l2OverFlow[1]);
  fill_overflow(l2PtRate_motherBin_3,l2OverFlow[2]);
  fill_overflow(l2PtRate_motherBin_4,l2OverFlow[3]);
  fill_overflow(l2PtRate_motherBin_5,l2OverFlow[4]);
  fill_overflow(l2PtRate_motherBin_6,l2OverFlow[5]);
  fill_overflow(l2PtRate_motherBin_10,l2OverFlow[9]);
  fill_overflow(l2PtRate_motherBin_12,l2OverFlow[11]);
  fill_overflow(l2PtRate_motherBin_13,l2OverFlow[12]);

  fill_overflow(l3PtRate_motherBin_1,l3OverFlow[0]);
  fill_overflow(l3PtRate_motherBin_2,l3OverFlow[1]);
  fill_overflow(l3PtRate_motherBin_3,l3OverFlow[2]);
  fill_overflow(l3PtRate_motherBin_4,l3OverFlow[3]);
  fill_overflow(l3PtRate_motherBin_5,l3OverFlow[4]);
  fill_overflow(l3PtRate_motherBin_6,l3OverFlow[5]);
  fill_overflow(l3PtRate_motherBin_10,l3OverFlow[9]);
  fill_overflow(l3PtRate_motherBin_12,l3OverFlow[11]);
  fill_overflow(l3PtRate_motherBin_13,l3OverFlow[12]);

  fill_overflow(l2IsoPtRate_motherBin_1,l2IsoOverFlow[0]);
  fill_overflow(l2IsoPtRate_motherBin_2,l2IsoOverFlow[1]);
  fill_overflow(l2IsoPtRate_motherBin_3,l2IsoOverFlow[2]);
  fill_overflow(l2IsoPtRate_motherBin_4,l2IsoOverFlow[3]);
  fill_overflow(l2IsoPtRate_motherBin_5,l2IsoOverFlow[4]);
  fill_overflow(l2IsoPtRate_motherBin_6,l2IsoOverFlow[5]);
  fill_overflow(l2IsoPtRate_motherBin_10,l2IsoOverFlow[9]);
  fill_overflow(l2IsoPtRate_motherBin_12,l2IsoOverFlow[11]);
  fill_overflow(l2IsoPtRate_motherBin_13,l2IsoOverFlow[12]);

  fill_overflow(l3PreIsoPtRate_motherBin_1,l3PreIsoOverFlow[0]);
  fill_overflow(l3PreIsoPtRate_motherBin_2,l3PreIsoOverFlow[1]);
  fill_overflow(l3PreIsoPtRate_motherBin_3,l3PreIsoOverFlow[2]);
  fill_overflow(l3PreIsoPtRate_motherBin_4,l3PreIsoOverFlow[3]);
  fill_overflow(l3PreIsoPtRate_motherBin_5,l3PreIsoOverFlow[4]);
  fill_overflow(l3PreIsoPtRate_motherBin_6,l3PreIsoOverFlow[5]);
  fill_overflow(l3PreIsoPtRate_motherBin_10,l3PreIsoOverFlow[9]);
  fill_overflow(l3PreIsoPtRate_motherBin_12,l3PreIsoOverFlow[11]);
  fill_overflow(l3PreIsoPtRate_motherBin_13,l3PreIsoOverFlow[12]);

  fill_overflow(l3IsoPtRate_motherBin_1,l3IsoOverFlow[0]);
  fill_overflow(l3IsoPtRate_motherBin_2,l3IsoOverFlow[1]);
  fill_overflow(l3IsoPtRate_motherBin_3,l3IsoOverFlow[2]);
  fill_overflow(l3IsoPtRate_motherBin_4,l3IsoOverFlow[3]);
  fill_overflow(l3IsoPtRate_motherBin_5,l3IsoOverFlow[4]);
  fill_overflow(l3IsoPtRate_motherBin_6,l3IsoOverFlow[5]);
  fill_overflow(l3IsoPtRate_motherBin_10,l3IsoOverFlow[9]);
  fill_overflow(l3IsoPtRate_motherBin_12,l3IsoOverFlow[11]);
  fill_overflow(l3IsoPtRate_motherBin_13,l3IsoOverFlow[12]);

  make_cd_from_histo(l2PtRate_motherBin_1,l2PtRate_motherBin_1_cd);
  make_cd_from_histo(l2PtRate_motherBin_2,l2PtRate_motherBin_2_cd);
  make_cd_from_histo(l2PtRate_motherBin_3,l2PtRate_motherBin_3_cd);
  make_cd_from_histo(l2PtRate_motherBin_4,l2PtRate_motherBin_4_cd);
  make_cd_from_histo(l2PtRate_motherBin_5,l2PtRate_motherBin_5_cd);
  make_cd_from_histo(l2PtRate_motherBin_6,l2PtRate_motherBin_6_cd);
  make_cd_from_histo(l2PtRate_motherBin_10,l2PtRate_motherBin_10_cd);
  make_cd_from_histo(l2PtRate_motherBin_12,l2PtRate_motherBin_12_cd);
  make_cd_from_histo(l2PtRate_motherBin_13,l2PtRate_motherBin_13_cd);

  make_cd_from_histo(l3PtRate_motherBin_1,l3PtRate_motherBin_1_cd);
  make_cd_from_histo(l3PtRate_motherBin_2,l3PtRate_motherBin_2_cd);
  make_cd_from_histo(l3PtRate_motherBin_3,l3PtRate_motherBin_3_cd);
  make_cd_from_histo(l3PtRate_motherBin_4,l3PtRate_motherBin_4_cd);
  make_cd_from_histo(l3PtRate_motherBin_5,l3PtRate_motherBin_5_cd);
  make_cd_from_histo(l3PtRate_motherBin_6,l3PtRate_motherBin_6_cd);
  make_cd_from_histo(l3PtRate_motherBin_10,l3PtRate_motherBin_10_cd);
  make_cd_from_histo(l3PtRate_motherBin_12,l3PtRate_motherBin_12_cd);
  make_cd_from_histo(l3PtRate_motherBin_13,l3PtRate_motherBin_13_cd);

  make_cd_from_histo(l2IsoPtRate_motherBin_1,l2IsoPtRate_motherBin_1_cd);
  make_cd_from_histo(l2IsoPtRate_motherBin_2,l2IsoPtRate_motherBin_2_cd);
  make_cd_from_histo(l2IsoPtRate_motherBin_3,l2IsoPtRate_motherBin_3_cd);
  make_cd_from_histo(l2IsoPtRate_motherBin_4,l2IsoPtRate_motherBin_4_cd);
  make_cd_from_histo(l2IsoPtRate_motherBin_5,l2IsoPtRate_motherBin_5_cd);
  make_cd_from_histo(l2IsoPtRate_motherBin_6,l2IsoPtRate_motherBin_6_cd);
  make_cd_from_histo(l2IsoPtRate_motherBin_10,l2IsoPtRate_motherBin_10_cd);
  make_cd_from_histo(l2IsoPtRate_motherBin_12,l2IsoPtRate_motherBin_12_cd);
  make_cd_from_histo(l2IsoPtRate_motherBin_13,l2IsoPtRate_motherBin_13_cd);

  make_cd_from_histo(l3PreIsoPtRate_motherBin_1,l3PreIsoPtRate_motherBin_1_cd);
  make_cd_from_histo(l3PreIsoPtRate_motherBin_2,l3PreIsoPtRate_motherBin_2_cd);
  make_cd_from_histo(l3PreIsoPtRate_motherBin_3,l3PreIsoPtRate_motherBin_3_cd);
  make_cd_from_histo(l3PreIsoPtRate_motherBin_4,l3PreIsoPtRate_motherBin_4_cd);
  make_cd_from_histo(l3PreIsoPtRate_motherBin_5,l3PreIsoPtRate_motherBin_5_cd);
  make_cd_from_histo(l3PreIsoPtRate_motherBin_6,l3PreIsoPtRate_motherBin_6_cd);
  make_cd_from_histo(l3PreIsoPtRate_motherBin_10,l3PreIsoPtRate_motherBin_10_cd);
  make_cd_from_histo(l3PreIsoPtRate_motherBin_12,l3PreIsoPtRate_motherBin_12_cd);
  make_cd_from_histo(l3PreIsoPtRate_motherBin_13,l3PreIsoPtRate_motherBin_13_cd);

  make_cd_from_histo(l3IsoPtRate_motherBin_1,l3IsoPtRate_motherBin_1_cd);
  make_cd_from_histo(l3IsoPtRate_motherBin_2,l3IsoPtRate_motherBin_2_cd);
  make_cd_from_histo(l3IsoPtRate_motherBin_3,l3IsoPtRate_motherBin_3_cd);
  make_cd_from_histo(l3IsoPtRate_motherBin_4,l3IsoPtRate_motherBin_4_cd);
  make_cd_from_histo(l3IsoPtRate_motherBin_5,l3IsoPtRate_motherBin_5_cd);
  make_cd_from_histo(l3IsoPtRate_motherBin_6,l3IsoPtRate_motherBin_6_cd);
  make_cd_from_histo(l3IsoPtRate_motherBin_10,l3IsoPtRate_motherBin_10_cd);
  make_cd_from_histo(l3IsoPtRate_motherBin_12,l3IsoPtRate_motherBin_12_cd);
  make_cd_from_histo(l3IsoPtRate_motherBin_13,l3IsoPtRate_motherBin_13_cd);

  l2PtRate->Add(l2PtRate_motherBin_1);
  l2PtRate->Add(l2PtRate_motherBin_2);
  l2PtRate->Add(l2PtRate_motherBin_3);
  l2PtRate->Add(l2PtRate_motherBin_4);
  l2PtRate->Add(l2PtRate_motherBin_5);
  l2PtRate->Add(l2PtRate_motherBin_6);
  l2PtRate->Add(l2PtRate_motherBin_10);
  l2PtRate->Add(l2PtRate_motherBin_12);
  l2PtRate->Add(l2PtRate_motherBin_13);

  l2PtRate_cd->Add(l2PtRate_motherBin_1_cd);
  l2PtRate_cd->Add(l2PtRate_motherBin_2_cd);
  l2PtRate_cd->Add(l2PtRate_motherBin_3_cd);
  l2PtRate_cd->Add(l2PtRate_motherBin_4_cd);
  l2PtRate_cd->Add(l2PtRate_motherBin_5_cd);
  l2PtRate_cd->Add(l2PtRate_motherBin_6_cd);
  l2PtRate_cd->Add(l2PtRate_motherBin_10_cd);
  l2PtRate_cd->Add(l2PtRate_motherBin_12_cd);
  l2PtRate_cd->Add(l2PtRate_motherBin_13_cd);

  l3PtRate->Add(l3PtRate_motherBin_1);
  l3PtRate->Add(l3PtRate_motherBin_2);
  l3PtRate->Add(l3PtRate_motherBin_3);
  l3PtRate->Add(l3PtRate_motherBin_4);
  l3PtRate->Add(l3PtRate_motherBin_5);
  l3PtRate->Add(l3PtRate_motherBin_6);
  l3PtRate->Add(l3PtRate_motherBin_10);
  l3PtRate->Add(l3PtRate_motherBin_12);
  l3PtRate->Add(l3PtRate_motherBin_13);

  l3PtRate_cd->Add(l3PtRate_motherBin_1_cd);
  l3PtRate_cd->Add(l3PtRate_motherBin_2_cd);
  l3PtRate_cd->Add(l3PtRate_motherBin_3_cd);
  l3PtRate_cd->Add(l3PtRate_motherBin_4_cd);
  l3PtRate_cd->Add(l3PtRate_motherBin_5_cd);
  l3PtRate_cd->Add(l3PtRate_motherBin_6_cd);
  l3PtRate_cd->Add(l3PtRate_motherBin_10_cd);
  l3PtRate_cd->Add(l3PtRate_motherBin_12_cd);
  l3PtRate_cd->Add(l3PtRate_motherBin_13_cd);

  l2IsoPtRate->Add(l2IsoPtRate_motherBin_1);
  l2IsoPtRate->Add(l2IsoPtRate_motherBin_2);
  l2IsoPtRate->Add(l2IsoPtRate_motherBin_3);
  l2IsoPtRate->Add(l2IsoPtRate_motherBin_4);
  l2IsoPtRate->Add(l2IsoPtRate_motherBin_5);
  l2IsoPtRate->Add(l2IsoPtRate_motherBin_6);
  l2IsoPtRate->Add(l2IsoPtRate_motherBin_10);
  l2IsoPtRate->Add(l2IsoPtRate_motherBin_12);
  l2IsoPtRate->Add(l2IsoPtRate_motherBin_13);

  l2IsoPtRate_cd->Add(l2IsoPtRate_motherBin_1_cd);
  l2IsoPtRate_cd->Add(l2IsoPtRate_motherBin_2_cd);
  l2IsoPtRate_cd->Add(l2IsoPtRate_motherBin_3_cd);
  l2IsoPtRate_cd->Add(l2IsoPtRate_motherBin_4_cd);
  l2IsoPtRate_cd->Add(l2IsoPtRate_motherBin_5_cd);
  l2IsoPtRate_cd->Add(l2IsoPtRate_motherBin_6_cd);
  l2IsoPtRate_cd->Add(l2IsoPtRate_motherBin_10_cd);
  l2IsoPtRate_cd->Add(l2IsoPtRate_motherBin_12_cd);
  l2IsoPtRate_cd->Add(l2IsoPtRate_motherBin_13_cd);

  l3PreIsoPtRate->Add(l3PreIsoPtRate_motherBin_1);
  l3PreIsoPtRate->Add(l3PreIsoPtRate_motherBin_2);
  l3PreIsoPtRate->Add(l3PreIsoPtRate_motherBin_3);
  l3PreIsoPtRate->Add(l3PreIsoPtRate_motherBin_4);
  l3PreIsoPtRate->Add(l3PreIsoPtRate_motherBin_5);
  l3PreIsoPtRate->Add(l3PreIsoPtRate_motherBin_6);
  l3PreIsoPtRate->Add(l3PreIsoPtRate_motherBin_10);
  l3PreIsoPtRate->Add(l3PreIsoPtRate_motherBin_12);
  l3PreIsoPtRate->Add(l3PreIsoPtRate_motherBin_13);

  l3PreIsoPtRate_cd->Add(l3PreIsoPtRate_motherBin_1_cd);
  l3PreIsoPtRate_cd->Add(l3PreIsoPtRate_motherBin_2_cd);
  l3PreIsoPtRate_cd->Add(l3PreIsoPtRate_motherBin_3_cd);
  l3PreIsoPtRate_cd->Add(l3PreIsoPtRate_motherBin_4_cd);
  l3PreIsoPtRate_cd->Add(l3PreIsoPtRate_motherBin_5_cd);
  l3PreIsoPtRate_cd->Add(l3PreIsoPtRate_motherBin_6_cd);
  l3PreIsoPtRate_cd->Add(l3PreIsoPtRate_motherBin_10_cd);
  l3PreIsoPtRate_cd->Add(l3PreIsoPtRate_motherBin_12_cd);
  l3PreIsoPtRate_cd->Add(l3PreIsoPtRate_motherBin_13_cd);

  l3IsoPtRate->Add(l3IsoPtRate_motherBin_1);
  l3IsoPtRate->Add(l3IsoPtRate_motherBin_2);
  l3IsoPtRate->Add(l3IsoPtRate_motherBin_3);
  l3IsoPtRate->Add(l3IsoPtRate_motherBin_4);
  l3IsoPtRate->Add(l3IsoPtRate_motherBin_5);
  l3IsoPtRate->Add(l3IsoPtRate_motherBin_6);
  l3IsoPtRate->Add(l3IsoPtRate_motherBin_10);
  l3IsoPtRate->Add(l3IsoPtRate_motherBin_12);
  l3IsoPtRate->Add(l3IsoPtRate_motherBin_13);

  l3IsoPtRate_cd->Add(l3IsoPtRate_motherBin_1_cd);
  l3IsoPtRate_cd->Add(l3IsoPtRate_motherBin_2_cd);
  l3IsoPtRate_cd->Add(l3IsoPtRate_motherBin_3_cd);
  l3IsoPtRate_cd->Add(l3IsoPtRate_motherBin_4_cd);
  l3IsoPtRate_cd->Add(l3IsoPtRate_motherBin_5_cd);
  l3IsoPtRate_cd->Add(l3IsoPtRate_motherBin_6_cd);
  l3IsoPtRate_cd->Add(l3IsoPtRate_motherBin_10_cd);
  l3IsoPtRate_cd->Add(l3IsoPtRate_motherBin_12_cd);
  l3IsoPtRate_cd->Add(l3IsoPtRate_motherBin_13_cd);

  //  double crossSection_mb = 75.28; // MinBias
  //  double crossSection_mb = 51.6; // ppMuX
  double crossSection_mb = 48.44; // ppMuX @ 7 TeV
  //  double crossSection_mb = 0.000317 ; // ttbar
  //  double filterEff = 1.0; //MinBias
  //  double filterEff = 0.00289; // ppMuX
  double filterEff = 0.00176; // ppMuX @ 7 TeV
  //   double filterEff = 0.33;  // ttbar

  double eventNumberUnit= 1000000;
  //   double filterEff = inclusiveppMuX_filterEff;
  double numberOfEvents = nEntries /eventNumberUnit /filterEff ;
  cout << "numberOfEvents = " << numberOfEvents << endl;

  double crossSection=crossSection_mb ;
  double L_E32= 100000; //10^32 cm-2.s-1 = 10^5 mb-1.s-1
  double L_E31= 10000; //10^31 cm-2.s-1 = 10^4 mb-1.s-1
  double L_E30= 1000; //10^30 cm-2.s-1 = 10^3 mb-1.s-1
  double L=L_E30;

  double mbInvToHz=L/eventNumberUnit;
  double rateFactor=crossSection / numberOfEvents * mbInvToHz ;
  cout << "rateFactor = " << rateFactor << endl;

  histDir->cd();
  scaleToRate("*PtRate*", rateFactor);


  histDir->Write("",TObject::kOverwrite);
  
  return 0;
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

double etaTrackIsoDepositCut (double eta) {
  if (abs(eta) < 0.435) return 1.1;
  else if (abs(eta) < 0.1305) return 1.1;
  else if (abs(eta) < 0.2175) return 1.1;
  else if (abs(eta) < 0.3045) return 1.1;
  else if (abs(eta) < 0.3915) return 1.2;
  else if (abs(eta) < 0.4785) return 1.1;
  else if (abs(eta) < 0.5655) return 1.2;
  else if (abs(eta) < 0.6525) return 1.1;
  else if (abs(eta) < 0.7395) return 1.2;
  else if (abs(eta) < 0.8265) return 1.0;
  else if (abs(eta) < 0.9135) return 1.1;
  else if (abs(eta) < 1.0005) return 1.0;
  else if (abs(eta) < 1.0875) return 1.0;
  else if (abs(eta) < 1.1745) return 1.1;
  else if (abs(eta) < 1.2615) return 1.0;
  else if (abs(eta) < 1.3485) return 1.0;
  else if (abs(eta) < 1.4355) return 1.1;
  else if (abs(eta) < 1.5225) return 0.9;
  else if (abs(eta) < 1.6095) return 1.1;
  else if (abs(eta) < 1.6965) return 0.9;
  else if (abs(eta) < 1.785) return 1.1;
  else if (abs(eta) < 1.88) return 1.0;
  else if (abs(eta) < 1.9865) return 1.0;
  else if (abs(eta) < 2.1075) return 0.9;
  else if (abs(eta) < 2.247) return 0.8;
  else if (abs(eta) < 2.411) return 0.1;
  else return 0.1;
}

double findMaxPt(vector<double> *pts) { // just a simple "find max"
  double maxPt = -999;
  for (int i = 0; i < pts->size(); i++) {
    if (pts->at(i) > maxPt) maxPt = pts->at(i);
  }
  return maxPt;
}

double findMaxPt(vector<double> *pts, vector<double> *etas) { // find max for l2
  double maxPt = -999;
  double etaCut = 2.5;
  for (int i = 0; i < pts->size(); i++) {
    if (pts->at(i) > maxPt && fabs(etas->at(i)) < etaCut) maxPt = pts->at(i);
  }
  return maxPt;
}

double findMaxPt(vector<double> *pts, vector<double> *etas, vector<double> *d0s) { // find max for l3
  double maxPt = -999;
  double etaCut = 1000000000; // junk, just to make sure I get no fucking seg faults
  double d0Cut = 2;
  for (int i = 0; i < pts->size(); i++) {
    if (pts->at(i) > maxPt && fabs(etas->at(i)) < etaCut && fabs(d0s->at(i)) < d0Cut) maxPt = pts->at(i);
  }
  return maxPt;
}

int findIndexOfMaxPt(vector<double> *pts) { // just a simple "find max"
  double maxPt = -999;
  int index = -999;
  for (int i = 0; i < pts->size(); i++) {
    if (pts->at(i) > maxPt) {
      index = i;
      maxPt = pts->at(i);
    }
  }
  return index;
}

int findIndexOfMaxPt(vector<double> *pts, vector<double> *etas) { // find max for l2
  double maxPt = -999;
  double etaCut = 2.5;
  int index = -999;
  for (int i = 0; i < pts->size(); i++) {
    if (pts->at(i) > maxPt && fabs(etas->at(i)) < etaCut) {
      index = i;
      maxPt = pts->at(i);
    }
  }
  return index;
}

int findIndexOfMaxPt(vector<double> *pts, vector<double> *etas, vector<double> *d0s) { // find max for l3
  double maxPt = -999;
  double etaCut = 999;
  double d0Cut = 2;
  int index = -999;
  for (int i = 0; i < pts->size(); i++) {
    if (pts->at(i) > maxPt && fabs(etas->at(i)) < etaCut && fabs(d0s->at(i)) < d0Cut) {
      index = i;
      maxPt = pts->at(i);
    }
  }
  return index;
}

void fill_overflow(TH1F* hist, double overflow) {
  int n_max = hist->GetNbinsX();
  hist->SetBinContent(n_max + 1, overflow);
}

// This one takes a hist and makes a cd.  Doesn't have to be same number of bins, but the number of bins should be the same
void make_cd_from_histo(TH1F* hist, TH1F* cd) {

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

  /*  double norm_factor = hist->GetEntries();
      
  if (norm_factor != 0) {
  cd->Scale(1/norm_factor);
  hist->Scale(1/norm_factor);
  }
  */
  return;
}

void scaleToRate(const char* collectionName, double rateFactor) {
  TRegexp reg(collectionName, kTRUE);
  
  //    gDirectory->ls();
  TList* list = gDirectory->GetList() ;
  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;
  
  while (obj = iter->Next()) {
    if (! obj->InheritsFrom(TH1::Class())) {
      //      cout << "bugger" << endl;
      continue;
    }
    
    
    TString name = obj->GetName();
    cout << "Testing name: " << name << " against " << collectionName << endl;
    
    if (TString(collectionName).MaybeRegexp()) {
      cout << "we have a possible match" << endl;
      cout << "Trying to match to " << TString(obj->GetName()) << endl;
      if (TString(obj->GetName()).Index(reg) < 0 ) {
	cout << "failure here.  Argument returns " << TString(obj->GetName()).Index(reg) << endl;
	continue;
      }
    }
    else if (! name.BeginsWith(collectionName)) continue;
    
    cout << "We're trying to scale" << name << endl;
    ((TH1*)obj)->Scale(rateFactor);
    
  }
  
}
