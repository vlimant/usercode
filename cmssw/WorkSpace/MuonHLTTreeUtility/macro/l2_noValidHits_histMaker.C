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
  
  TH1F *DR_toL2WithValidHits = new TH1F("DR_toL2WithValidHits","#Delta R to nearest L2",150,0,6);
  TH1F *l2_seeding_success = new TH1F("l2_seeding_success","check if l2 seeds l3",2,-0.5,1.5);
  TH1F *pTNearestL2 = new TH1F("pTNearestL2","log10(p_{T}) of nearest L2",100,-3,6);
  TH1F *pTProblemL2 = new TH1F("pTProblemL2","p_{T} of L2 w/o valid hits",100,-3,6);
  TH1F *fractionInvalidHits = new TH1F("fractionInvalidHits","n_invalid/n_hits",101,0,1.01);
  TH1F *nHits_all = new TH1F("nHits_all","n_hits for all L2",55,-0.5,54.5);
  TH1F *nHits_noValid = new TH1F("nHits_noValid","n_hits L2 with no valid hits",13,-0.5,12.5);
  TH1F *calIsoDeposit_ = new TH1F("calIsoDeposit_all","l2 cal dep, all",100,0,10);

  TH1F *l1SeedQuality = new TH1F("l1SeedQuality","quality of L1 seed",8,-0.5,7.5);
  TH1F *l1SeedEta = new TH1F("l1SeedEta","#eta of L1 seed",50,-2.5,2.5);
  TH1F *l1SeedPt = new TH1F("l1SeedPt","p_{T} of L1 seed",150,0,150);

  TH1F *l2IsAssoc = new TH1F("l2IsAssoc","l2->sim assoc",2,-0.5,1.5);
  TH1F *l2Parent = new TH1F("l2Parent","mother bin, assoc w/ no valid hits",13,-0.5,12.5);
  TH1F *l2ToSimDR = new TH1F("l2ToSimDR","#Delta R (sim, L2), no valid hits",100,0,1);
  TH1F *l2ParentPt = new TH1F("l2ParentPt","pT of assoc sim #mu",100,0,100);

  TH1F *l3SeededPt = new TH1F("l3SeededPt","log10(L3 p_{T}) from N.V.H. L2",100,-3,6);
  TH1F *l3SeededEta = new TH1F("l3SeededEta","L3 #eta from N.V.H. L2",50,-2.5,2.5);
  TH1F *l3SeededD0 = new TH1F("l3SeededD0","L3 d_{0} from N.V.H. L2",100,-0.06,0.06);

  TH2F *DR_toL2WithValidHits_vsEta = new TH2F("DR_toL2WithValidHits_vsEta","#Delta R to nearest L2 vs #eta",100,-2.5,2.5,150,0,6);
  TH2F *pTProblemL2_vsEta = new TH2F("pTProblemL2_vsEta","log10(p_{T}) of L2 w/o valid hits vs #eta",100,-2.5,2.5,100,-3,6);
  TH2F *fractionInvalidHits_vsEta = new TH2F("fractionInvalidHits_vsEta","n_invalid/n_hits vs #eta",100,-2.5,2.5,101,0,1.01);
  TH2F *nHits_all_vsEta = new TH2F("nHits_all_vsEta","n_hits for all L2 vs #eta",100,-2.5,2.5,55,-0.5,54.5);
  TH2F *nHits_noValid_vsEta = new TH2F("nHits_noValid_vsEta","n_hits for all L2 vs #eta",100,-2.5,2.5,13,-0.5,12.5);

  TH2F *l1SeedQuality_vsEta = new TH2F("l1SeedQuality_vsEta","L1 quality vs #eta",50,-2.5,2.5,8,-0.5,7.5);
  TH2F *l1SeedPt_vsEta = new TH2F("l1SeedPt_vsEta","L1 pt vs #eta",50,-2.5,2.5,150,0,150);

  int nEntries = tree->GetEntries();
				     
  for( int iEntry = 0; iEntry < nEntries; iEntry++) {
    if (iEntry%1000 == 0) cout << "Event " << iEntry << endl;
    tree->GetEntry(iEntry);
    
    // First we fill in the sim information
    //    cout << "Loop over Sim" << endl;
    map<int,int>::iterator fuckIt;
    map<int,int>::iterator fuckIt2;  
    map<int,int>::iterator fuckIt3;

    for (fuckIt = (*l2NMuHits).begin(); fuckIt != (*l2NMuHits).end(); fuckIt++) {
      int index = fuckIt->first;     
      int indexL1 = (*indexL1SeedingL2).at(index);
      fractionInvalidHits->Fill(1 - ((double)(fuckIt->second)/(double)(*l2NHits).at(index)));
      fractionInvalidHits_vsEta->Fill((*l2Eta).at(index),1 - ((double)(fuckIt->second)/(double)(*l2NHits).at(index)));
      nHits_all->Fill((*l2NHits).at(index));
      nHits_all_vsEta->Fill((*l2Eta).at(index),(*l2NHits).at(index));

      if (fuckIt->second == 0) {
	l2IsAssoc->Fill((*l2IsAssociated).at(index));
	
	if ((*l2IsAssociated).at(index)) {
	  l2Parent->Fill((*l2MotherBinNumber).at(index));
	  l2ToSimDR->Fill(-1 * (*l2AssociationVar).at(index));
	  l2ParentPt->Fill((*l2AssociatedSimMuonPt).at(index));
	}
	
	double DR = dist_to_nearestL2((*l2Eta).at(index),(*l2Phi).at(index), index, l2Eta, l2Phi);
	l2_seeding_success->Fill((*l2SeedsL3).at(index));

	for (fuckIt2 = (*l3NMuHits).begin(); fuckIt2 != (*l3NMuHits).end(); fuckIt2++) {
	  int indexL3 = fuckIt2->first;
	  int indexL2 = (*indexL2SeedingL3).at(indexL3);
	  if (indexL2 == index) {
	    l3SeededPt->Fill(log10((*l3Pt).at(indexL3)));
	    l3SeededEta->Fill((*l3Eta).at(indexL3));
	    l3SeededD0->Fill((*l3D0).at(indexL3));
	  }
	}
	
	pTProblemL2->Fill(log10((*l2Pt).at(index)));
	pTProblemL2_vsEta->Fill((*l2Eta).at(index),log10((*l2Pt).at(index)));
	nHits_noValid->Fill((*l2NHits).at(index));
	nHits_noValid_vsEta->Fill((*l2Eta).at(index),(*l2NHits).at(index));
	if (DR < 1.0) {
	  int indexNearest = index_of_nearestL2((*l2Eta).at(index),(*l2Phi).at(index), index, l2Eta, l2Phi);
	  pTNearestL2->Fill(log10((*l2Pt).at(indexNearest)));
	}

	l1SeedPt->Fill((*l1Pt).at(indexL1));
	l1SeedEta->Fill((*l1Eta).at(indexL1));
	l1SeedQuality->Fill((*l1Quality).at(indexL1));

	l1SeedQuality_vsEta->Fill((*l1Eta).at(indexL1),(*l1Quality).at(indexL1));
	l1SeedPt_vsEta->Fill((*l1Eta).at(indexL1),(*l1Pt).at(indexL1));

	DR_toL2WithValidHits->Fill(DR);
	DR_toL2WithValidHits_vsEta->Fill((*l2Eta).at(index),DR);	
	//	cout << "DR, index of nearest = " << DR <<" " << indexNearest << endl;
      }
    }

  } // loop over events

  histFile->Write("",TObject::kOverwrite);
  return 0;
}

double dist_to_nearestL2(double thisEta, double thisPhi, int thisIndex, vector<double> *allEta, vector<double> *allPhi) {

  double DeltaR = 999.9;
  
  for (int i = 0; i < (*allEta).size(); i++) {
    if (i != thisIndex) {
      double temp_deleta = thisEta - (*allEta).at(i);
      double temp_delphi = thisPhi - (*allPhi).at(i);      
      if (temp_delphi > 3.141592653) {
	temp_delphi = (2 * 3.141592653) - temp_delphi;
      }
      double temp_DeltaR = sqrt((temp_deleta * temp_deleta) + (temp_delphi * temp_delphi));
      
      if (temp_DeltaR < DeltaR) {
	DeltaR = temp_DeltaR;
      }
    }
  }
  //  cout << "DeltaR = " << DeltaR << endl;
  //  cout << "returning" << endl;
  return DeltaR;
}

double index_of_nearestL2(double thisEta, double thisPhi, int thisIndex, vector<double> *allEta, vector<double> *allPhi) {

  double DeltaR = 999.9;
  int index = -1;

  for (int i = 0; i < (*allEta).size(); i++) {
    if (i != thisIndex) {
      double temp_deleta = thisEta - (*allEta).at(i);
      double temp_delphi = thisPhi - (*allPhi).at(i);      
      if (temp_delphi > 3.141592653) {
	temp_delphi = (2 * 3.141592653) - temp_delphi;
      }
      double temp_DeltaR = sqrt((temp_deleta * temp_deleta) + (temp_delphi * temp_delphi));
      
      if (temp_DeltaR < DeltaR) {
	DeltaR = temp_DeltaR;
	index = i;
      }
    }
  }
  //  cout << "DeltaR = " << DeltaR << endl;
  //  cout << "returning" << endl;
  return index;
}
