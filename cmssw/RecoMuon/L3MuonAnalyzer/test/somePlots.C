void doEta(){

  TCut where="RoadCandidateArray.best.nRoadHit<=5";

  TH1F * h1 = new TH1F("h1_eta","Simulated muon |#eta|",30,0,2.5);
  h1->SetXTitle("muon |#eta_{sim}|");

  TH2F * h2 = new TH2F("h2_eta","Simulated muon |#eta| versus reconstructed muon |#eta|",30,0,2.5,30,0,2.5);
  h2->SetXTitle("muon |#eta_{sim}|");
  h2->SetYTitle("muon |#eta_{best}|");

  data->Draw("RoadCandidateArray.sim.eta(0)>>"+h1->GetName(),where && "hasBest","");
  data->Draw("RoadCandidateArray.best.updated.eta(0):RoadCandidateArray.sim.eta(0)>>"+h2->GetName(),where && "hasBest","");

  h1->Scale(1/(float)h1->GetMaximum());

  h2->Draw("colz");
  h1->Draw("same");
}


TCanvas * seperate(TString tag , TString var , TString xlabel , int Nb, double min, double max, TCut separation="RoadCandidateArray.best.nRoadHit>4", TCut selection ="hasBest"){
  gStyle->SetOptStat(0);
  TCanvas * c = new TCanvas ("c_"+tag,TString("comparison with seperation: ")+separation.GetTitle(),400,400);
  c->cd()->SetGrid();

  TH1F * h_below = new TH1F("h_below_"+tag,"distribution of "+xlabel+" not verfying "+separation.GetTitle(),Nb,min,max);
  h_below->SetXTitle(xlabel);
  h_below->SetLineColor(4);
  h_below->SetLineWidth(2);

  TH1F * h_above = new TH1F("h_above_"+tag,"distribution of "+xlabel+" verfying "+separation.GetTitle(),Nb,min,max);
  h_above->SetXTitle(xlabel);
  h_above->SetLineColor(2);
  h_above->SetLineWidth(2);

  data->Draw(var+">>"+h_below->GetName(), !separation && selection);
  data->Draw(var+">>"+h_above->GetName(), separation && selection);

  h_below->Scale(1/h_below->Integral());
  h_above->Scale(1/h_above->Integral());


  TLegend * leg = new TLegend(0.7,0.8,0.99,0.90);
  leg->AddEntry(h_below,Form("(%s)=False",separation.GetTitle()),"l");
  leg->AddEntry(h_above,Form("(%s)=True",separation.GetTitle()),"l");

  h_below->Draw();
  h_above->Draw("same");
  leg->Draw();

  return c;
}

void doSeparation()
{
  TCut separation = "RoadCandidateArray.best.nRoadHit>4";
  
  seperate("seta","RoadCandidateArray.sim.eta(0)","|#eta_{sim}|", 30,0,3,separation,"hasBest && hasSim");
  seperate("sphi","RoadCandidateArray.sim.phi()","#phi_{sim}", 30,-TMath::Pi(),TMath::Pi(),separation,"hasBest && hasSim");

  seperate("ratio","RoadCandidateArray.nRoadHit/RoadCandidateArray.best.nRoadHit","#frac{Nb.road}{Nb.best}", 30,0,6,separation,"hasBest");

  seperate("nmuRH","RoadCandidateArray.Nrechits","Nb. muon RH",71,0.5,70.5,separation,"hasBest"); 

  seperate("ntkRH","RoadCandidateArray.best.nRoadHit","Nb. Tk RH",20,0.5,19.5,separation,"hasBest");
  seperate("nroadRH","RoadCandidateArray.nRoadHit","Nb. road RH",20,0.5,19.5,separation,"hasBest");

  seperate("pull_sta_x","(RoadCandidateArray.Cpca.x-RoadCandidateArray.sim.x)/sqrt(RoadCandidateArray.Cpca.c11)","#frac{#Delta_{x}}{#sigma_{x}}",20,-10,10,separation,"hasBest && hasSim && hasCpca");
  seperate("pull_sta_y","(RoadCandidateArray.Cpca.y-RoadCandidateArray.sim.y)/sqrt(RoadCandidateArray.Cpca.c22)","#frac{#Delta_{y}}{#sigma_{y}}",20,-10,10,separation,"hasBest && hasSim && hasCpca");
  seperate("pull_sta_z","(RoadCandidateArray.Cpca.z-RoadCandidateArray.sim.z)/sqrt(RoadCandidateArray.Cpca.c33)","#frac{#Delta_{z}}{#sigma_{z}}",20,-10,10,separation,"hasBest && hasSim && hasCpca");

  seperate("pull_sta_px","(RoadCandidateArray.Cpca.px-RoadCandidateArray.sim.px)/sqrt(RoadCandidateArray.Cpca.c44)","#frac{#Delta_{px}}{#sigma_{px}}",30,-15,15,separation,"hasBest && hasSim && hasCpca");
  seperate("pull_sta_py","(RoadCandidateArray.Cpca.py-RoadCandidateArray.sim.py)/sqrt(RoadCandidateArray.Cpca.c55)","#frac{#Delta_{py}}{#sigma_{py}}",30,-15,15,separation,"hasBest && hasSim && hasCpca");
  seperate("pull_sta_pz","(RoadCandidateArray.Cpca.pz-RoadCandidateArray.sim.pz)/sqrt(RoadCandidateArray.Cpca.c66)","#frac{#Delta_{pz}}{#sigma_{pz}}",30,-15,15,separation,"hasBest && hasSim && hasCpca");

  seperate("pull_x","(RoadCandidateArray.best.updated.x-RoadCandidateArray.sim.x)/sqrt(RoadCandidateArray.best.updated.c11)","#frac{#Delta_{x}}{#sigma_{x}}",20,-10,10,separation,"hasBest && hasSim");
  seperate("pull_y","(RoadCandidateArray.best.updated.y-RoadCandidateArray.sim.y)/sqrt(RoadCandidateArray.best.updated.c22)","#frac{#Delta_{y}}{#sigma_{y}}",20,-10,10,separation,"hasBest && hasSim");
  seperate("pull_z","(RoadCandidateArray.best.updated.z-RoadCandidateArray.sim.z)/sqrt(RoadCandidateArray.best.updated.c33)","#frac{#Delta_{z}}{#sigma_{z}}",20,-10,10,separation,"hasBest && hasSim");

  seperate("pull_px","(RoadCandidateArray.best.updated.px-RoadCandidateArray.sim.px)/sqrt(RoadCandidateArray.best.updated.c44)","#frac{#Delta_{px}}{#sigma_{px}}",30,-15,15,separation,"hasBest && hasSim");
  seperate("pull_py","(RoadCandidateArray.best.updated.py-RoadCandidateArray.sim.py)/sqrt(RoadCandidateArray.best.updated.c55)","#frac{#Delta_{py}}{#sigma_{py}}",30,-15,15,separation,"hasBest && hasSim");
  seperate("pull_pz","(RoadCandidateArray.best.updated.pz-RoadCandidateArray.sim.pz)/sqrt(RoadCandidateArray.best.updated.c66)","#frac{#Delta_{pz}}{#sigma_{pz}}",30,-15,15,separation,"hasBest && hasSim");

}
