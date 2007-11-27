#include "TH1F.h"
#include "TAxis.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include <vector>
using namespace std;
#include "TString.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGraph.h"
#include <iostream>

TVirtualPad * JoinPad(TCanvas * c , int ip1, int ip2){
  double xl,yl,xh,yh;
  xl=c->GetPad(ip1)->GetXlowNDC();
  yl=c->GetPad(ip1)->GetYlowNDC();
  xh=c->GetPad(ip2)->GetXlowNDC() + c->GetPad(ip2)->GetWNDC();
  yh=c->GetPad(ip2)->GetYlowNDC() + c->GetPad(ip2)->GetHNDC();
  c->GetPad(ip2)->Delete();
  c->GetPad(ip1)->SetPad(xl,yl,xh,yh); 
  return c->cd(ip1);}

double GetMean(TH1F* h, double xmin){
  int Nmin = h->FindBin(xmin);
  int Nmax = h->GetNbinsX();
  double sum=0;
  double entries=0;
  for (int ib=Nmin;ib<=Nmax;++ib){
    entries+= h->GetBinContent(ib);
    sum+=h->GetBinContent(ib)*h->GetBinCenter(ib);}
  double mean=0;
  if (entries!=0){
    mean = sum/entries;}
  else {mean=0;}
  return mean;}




TPaveText * pave(TH1F * h,double min){
  TPaveText * p= new TPaveText(0.6,0.6,0.90,0.90,"NDC");
  p->AddText("Module:");
  p->AddText(h->GetName());
  p->AddText("Global mean: ");
  p->AddText(Form("%6.2f ms",h->GetMean()));
  if (min!=0){  
    p->AddText(Form(">%3.f ms mean: ",min));
    p->AddText(Form("%6.2f ms",GetMean(h,min)));}
  return p;}

class module {
 public:
  module(TString n,double mi,double ma=-1,bool l=false):name(n),min(mi),max(ma),log(l){;}
  TString name;
  double min;//used for the mean
  double max;//used to resize the histogram
  bool log;
};

vector<module> ckfPIXEL;
vector<module> muRS;
vector<module> ckfRS;
vector<module> ckfRSOI;

vector<module> ckfPIXELSigma1;
vector<module> ckfPIXELSigma3;
vector<module> ckfPIXELSigma6;
vector<module> ckfPIXELSigma10;

int bsize=300;
int Nbin=100;
TFile * masterFile=NULL;
TH1F * get(TString hname){return (TH1F*) masterFile->Get(hname);}
void SetX(TH1F* h , float xmax){
  int imax=h->GetXaxis()->FindBin(xmax);
  if (imax>Nbin){
    //need to rebin
    h->Rebin(imax/Nbin);
    imax=h->GetXaxis()->FindBin(xmax);
    h->GetXaxis()->SetRange(0,imax);}
  else{
    h->GetXaxis()->SetRange(0,imax);}}

TGraph * truncatedMean(TH1F* h,double max){
  TGraph * gr = new TGraph();
  gr->SetName(Form("cumul_mean_%s",h->GetName()));
  double past=-1;
  int g=0;
  const int iP=Nbin;
  int imax=h->FindBin(max);
  int iS = std::max(1,imax/iP);
  cout <<"calculating truncated mean for module "<<h->GetName()<<"..."<<endl;
  for (int ib=1;ib<=imax;ib+=iS){
    double xmin=h->GetBinCenter(ib);
    //    cout <<ib<<" "<<xmin<<"..."<<endl;
    double mean=GetMean(h,xmin);
    //    cout<<mean<<endl;
    //    if (mean!=past)
      {gr->SetPoint(g++,xmin,mean);}
    past=mean;
  }
  //  gr->SetMarkerStyle(4);
  cout<<"done"<<endl;
  return gr;}

TList * aPath(TString tag, TString label,vector<module> modules){
  TList * result=new TList();

  bool canProceed=false;
  int Ndiv=0;
  for (int i=0;i!=(int)modules.size();++i){if(get(modules[i].name)) {canProceed=true;Ndiv++;}}
  if (!canProceed) return result;

  TCanvas * c = new TCanvas("c_"+tag,"timing for "+label+"path",modules.size()*bsize,bsize);
  c->Divide(modules.size(),1);
  //  c->Divide(Ndiv,1);
  TCanvas * d = new TCanvas("d_"+tag,"truncated mean for "+label+"path",modules.size()*bsize,bsize);
  d->Divide(modules.size(),1);
  //  d->Divide(Ndiv,1);

  result->Add(c);
  result->Add(d);
  
  for (int i=0;i!=(int)modules.size();++i){
    TH1F * h = get(modules[i].name);
    if (!h) {
      //      i--;
      continue;}
    c->cd(i+1)->SetGrid();
    if (modules[i].log) 
      c->cd(i+1)->SetLogy();

    if (modules[i].max!=-1){
      //resize the axis
      SetX(h,modules[i].max);}
    h->SetXTitle("time [ms]");
    h->Draw();
    //plot the means
    pave(h,modules[i].min)->Draw();
    d->cd(i+1)->SetGrid();
    TGraph * tr=truncatedMean(h,modules[i].max);
    tr->Draw("apl");
    tr->GetXaxis()->SetTitle("min time [ms]");
    tr->GetYaxis()->SetTitle("mean time starting from min [ms]");
    tr->GetXaxis()->SetRange(0,tr->GetXaxis()->FindBin(modules[i].max));
  }
  return result;}

void Init(char * filename){
  //some style
  gStyle->SetOptStat(false);
  gStyle->SetOptTitle(false);
  gStyle->SetLineWidth(2);

  //the file name
  masterFile=TFile::Open(filename);

  float seedmin=0;
  float seedmax=100;

  float trackcandmin=0;
  float trackcandmax=3000;  

  float trackpromin=0;
  float trackpromax=100;

  float assocmin=0;
  float assocmax=10;

  //ckf path
  ckfPIXEL.push_back(module("siPixelClusters",0,2500));
  ckfPIXEL.push_back(module("siPixelRecHits",0,100));
  //  ckfPIXEL.push_back(module("muPIXELseed",seedmin,seedmax,true));
  ckfPIXEL.push_back(module("muPIXELseed",seedmin,6000,true));
  ckfPIXEL.push_back(module("ckfPIXELtrackcandidate",trackcandmin,trackcandmax,true));
  ckfPIXEL.push_back(module("ckfPIXELtrackproducer",trackpromin,trackpromax));

  //muRS path
  muRS.push_back(module("muRSseed",seedmin,seedmax));
  muRS.push_back(module("muRStrackcandidate",trackcandmin,trackcandmax,true));
  muRS.push_back(module("muRStrackproducer",trackpromin,trackpromax));
  
  //ckfRS path
  ckfRS.push_back(module("muRSseed",seedmin,seedmax));
  ckfRS.push_back(module("ckfRStrackcandidate",trackcandmin,trackcandmax,true));
  ckfRS.push_back(module("ckfRStrackproducer",trackpromin,trackpromax));

  //ckfRS-OI path
  ckfRSOI.push_back(module("muRSseedOI",seedmin,seedmax));
  ckfRSOI.push_back(module("ckfRStrackcandidateOI",trackcandmin,trackcandmax,true));
  ckfRSOI.push_back(module("ckfRStrackproducerOI",trackpromin,trackpromax));

  //ckf path sigma 1
  ckfPIXELSigma1.push_back(module("siPixelClusters",0,2500));
  ckfPIXELSigma1.push_back(module("siPixelRecHits",0,100));
  ckfPIXELSigma1.push_back(module("muPIXELseedSigma1",seedmin,seedmax,true));
  ckfPIXELSigma1.push_back(module("ckfPIXELtrackcandidateSigma1",trackcandmin,trackcandmax,true));
  ckfPIXELSigma1.push_back(module("ckfPIXELtrackproducerSigma1",trackpromin,trackpromax));
  ckfPIXELSigma1.push_back(module("ckfPIXELassocSigma1",assocmin,assocmax));

  //ckf path sigma 3
  ckfPIXELSigma3.push_back(module("siPixelClusters",0,2500));
  ckfPIXELSigma3.push_back(module("siPixelRecHits",0,100));
  ckfPIXELSigma3.push_back(module("muPIXELseedSigma3",seedmin,seedmax,true));
  ckfPIXELSigma3.push_back(module("ckfPIXELtrackcandidateSigma3",trackcandmin,trackcandmax,true));
  ckfPIXELSigma3.push_back(module("ckfPIXELtrackproducerSigma3",trackpromin,trackpromax));
  ckfPIXELSigma3.push_back(module("ckfPIXELassocSigma3",assocmin,assocmax));

  //ckf path sigma 6
  ckfPIXELSigma6.push_back(module("siPixelClusters",0,2500));
  ckfPIXELSigma6.push_back(module("siPixelRecHits",0,100));
  ckfPIXELSigma6.push_back(module("muPIXELseedSigma6",seedmin,seedmax,true));
  ckfPIXELSigma6.push_back(module("ckfPIXELtrackcandidateSigma6",trackcandmin,trackcandmax,true));
  ckfPIXELSigma6.push_back(module("ckfPIXELtrackproducerSigma6",trackpromin,trackpromax));
  ckfPIXELSigma6.push_back(module("ckfPIXELassocSigma6",assocmin,assocmax));

  //ckf path sigma 10
  ckfPIXELSigma10.push_back(module("siPixelClusters",0,2500));
  ckfPIXELSigma10.push_back(module("siPixelRecHits",0,100));
  ckfPIXELSigma10.push_back(module("muPIXELseedSigma10",seedmin,seedmax,true));
  ckfPIXELSigma10.push_back(module("ckfPIXELtrackcandidateSigma10",trackcandmin,trackcandmax,true));
  ckfPIXELSigma10.push_back(module("ckfPIXELtrackproducerSigma10",trackpromin,trackpromax));
  ckfPIXELSigma10.push_back(module("ckfPIXELassocSigma10",assocmin,assocmax));
}

void ShowMe(char * fname){
  Init(fname);
  
  aPath("ckfPIXEL","PIXEL seed + Ckf TC",ckfPIXEL);
  aPath("muRS","muon road seed + muon road TC",muRS);
  aPath("ckfRS","muon road seed + Ckf TC",ckfRS);
  aPath("ckfRSOI","muon road seed OUT + Ckf TC",ckfRSOI);

  aPath("ckfPIXELSigma1","PIXEL seed (1#sigma) + Ckf TC",ckfPIXELSigma1);
  aPath("ckfPIXELSigma3","PIXEL seed (3#sigma) + Ckf TC",ckfPIXELSigma3);
  aPath("ckfPIXELSigma6","PIXEL seed (6#sigma) + Ckf TC",ckfPIXELSigma6);
  aPath("ckfPIXELSigma10","PIXEL seed (10#sigma) + Ckf TC",ckfPIXELSigma10);


  return;
}
