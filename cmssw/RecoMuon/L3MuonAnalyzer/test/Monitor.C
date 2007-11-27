#include "TString.h"
#include "TTree.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TCut.h"
#include "TF1.h"
#include "TList.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TTreeFormula.h"
#include "TDirectory.h"
#include "TObject.h"
#include "TFile.h"
#include "TChain.h"
#include "TProfile3D.h"
#include "TSystem.h"
#include "TKey.h"
#include "TVirtualPad.h"

#include <iostream>
#include <fstream>
using namespace std;


const int xblock=300;
const int yblock=300;

TTree * data=NULL;
TCut globalCut="(1==1)";


void applyCut(TTree * t , TString cut="old")
{
  /*  static TEventList * oldList = t->GetEventList();
      if (cut=="old")
      t->SetEventList(oldList);
  */

  cout<<"applying cut "<<cut<<endl;
  /*
    TEventList * l=new TEventList("applyCutList","list of events");
    t->Draw(TString(">>")+l->GetName(),cut);
    cout<<"number of events selected: "<<l->GetN();

    t->SetEventList(l);
  */

  globalCut=cut;
}

TTree * openMany(char * filelist,int N=-1,TString withCut="")
{
  TDirectory *oldDir = gDirectory;
  ifstream  file(filelist);
  TChain * chain = new TChain("dump","");
  string filename;
  if (N==-1) N=100000;
  int i=0;
  while (!file.eof() && i!=N)
    {
      file >> filename;
      if (filename.size()!=0)
	{
	  cout<<"adding "<<filename<<" to the chain"<<endl;
	  chain->Add(filename.c_str());
	}
      i++;
    }
  data = (TTree*) chain;
  if (withCut!="")
    applyCut(data,withCut);

  oldDir->cd();
  return data;
}

TTree * open(char * filename,TString withCut="")
{ cout <<" opening :"<<filename<<"..."; 
  TDirectory *oldDir = gDirectory; 
  data=(TTree*)( TFile::Open(filename)->Get("dump"));
  oldDir->cd();
  cout<<"...done"<<endl;
  if (withCut!="")
    applyCut(data,withCut);
  return data;}


//transform a profile to a TH1D of the RMS in each bin of the profile
//very handy for resolution plots
TH1D * RMS(const TProfile * pf)
{
  TH1D * h= new TH1D(pf->GetName()+TString("_RMS"),
		     pf->GetTitle()+TString(" RMS"),
		     pf->GetNbinsX(),
		     pf->GetXaxis()->GetXmin(),
		     pf->GetXaxis()->GetXmax());
		     
  for (int ib=0;ib!=pf->GetNbinsX()+1;ib++)
    {
      double Nb = pf->GetBinEntries(ib);
      double rms = pf->GetBinError(ib);
      if (Nb!=0)
	h->SetBinError(ib,rms/(2*sqrt((double)Nb)));
      else
	h->SetBinError(ib,0);

      h->SetBinContent(ib,rms);
    }

  gROOT->cd();
  return h;
}


TH1D * CoreWidth(const TH2D * h2,double limit=0.3)
{
  TH1D * h= new TH1D(h2->GetName()+TString("_CoreW"),
		     h2->GetTitle()+TString(" core width"),
		     h2->GetNbinsX(),
		     h2->GetXaxis()->GetXmin(),
		     h2->GetXaxis()->GetXmax());
  TF1 f("f","gaus",-10,10);
  for (int ib=1;ib!=h2->GetNbinsX();ib++)
    {
      TH1D * ph=h2->ProjectionY("_py",ib,ib);
      //      int status=
      ph->Fit(&f,"","",-limit,limit);
      /* while(status!=4)
	 {
	 double scale=2;
	 status=h2->ProjectionY("_py",ib,ib)->Fit(&f,"","",-scale*limit,scale*limit);
	 scale++;
	 }*/
      double val=f.GetParameter(2);
      double e=f.GetParError(2);
      if (val>ph->GetRMS() && ph->GetEntries()!=0)
	{ val=ph->GetRMS();
	  e=ph->GetRMS()/sqrt(4*ph->GetEntries());}

      h->SetBinContent(ib,val);
      h->SetBinError(ib,e);
    }
  return h;
}


/*
  enum layer {TIB_1,TIB_2,TIB_3,TIB_4,
  TOB_1,TOB_2,TOB_3,TOB_4,TOB_5,TOB_6,
  TID_1,TID_2,TID_3,
  TEC_1,TEC_2,TEC_3,TEC_4,TEC_5,TEC_6,TEC_7,TEC_8,TEC_9};
  double radiuses[10]={24,32,40,48,62,68,78,85,95,106};
  double Zes[12]={77.5, 90 , 105,125 , 138 , 151 , 168 , 185 , 205, 225 , 248 , 266};
*/
TString lnames[10+12]={
  "TIB_1","TIB_2","TIB_3","TIB_4",
  "TOB_1","TOB_2","TOB_3","TOB_4","TOB_5","TOB_6",
  "TID_1","TID_2","TID_3",
  "TEC_1","TEC_2","TEC_3","TEC_4","TEC_5","TEC_6","TEC_7","TEC_8","TEC_9"};


TCut layercut(int l)
{
  TString layer ="RoadCandidateArray[].RoadSimHitArray[].simhit.layer(";
  layer+=l;
  layer+=")";
  return TCut(layer);
  /*
    double radius=radiuses[min((int)l,9)];
    double z=Zes[max(0,(int)l-10)];
    TString layer;
    if (l<10)
    {
    layer="(abs(RoadCandidateArray[].RoadSimHitArray[].rechit.r()-";
    layer+=radius;
    layer+=")<4)";
    layer+="&&(abs(RoadCandidateArray[].RoadSimHitArray[].rechit.z)<";
    if (l<4)
    layer+=Zes[0]-8; //firs TID z
    else
    layer+=Zes[3]-8; //first TEC z
    layer+=")";
    }
    else
    {
    layer="(abs(RoadCandidateArray[].RoadSimHitArray[].rechit.z-";
    layer+=z;
    layer+=")<8)";
    if (l<13)
    {
    layer+="&&(RoadCandidateArray[].RoadSimHitArray[].rechit.r()<";
    layer+=radiuses[4]-5;//first TOB layer
    layer+=")";
    }
    }
    return TCut(layer);*/
}



TCanvas * doRoadSize(TCut where, TString tag,int Nbinx, double xmin,double xmax,int Nbiny,double ymin,double ymax)
{
  cout<<"doRoadSize: "<<where<<" "<<tag<<endl;
  gStyle->SetOptStat(true);
  gStyle->SetOptFit(true);
  
  TH1F * road_x = new TH1F(tag+"_roadwidth_x","road x width at "+tag,Nbinx,xmin,xmax);
  road_x->SetXTitle("#sigma^{road}_{x^{local}} [cm]");
  TH1F * road_y = new TH1F(tag+"_roadwidth_y","road y width at "+tag,Nbiny,ymin,ymax);
  road_y->SetXTitle("#sigma^{road}_{y^{local}} [cm]");
  
  TCanvas * c = new TCanvas(tag+"_roadidual","road width at "+tag+"",2*xblock,1*yblock);
  data->Draw(TString("sqrt(RoadCandidateArray[].RoadSimHitArray[].road.v11_l)>>")+road_x->GetName(),where && globalCut);
  data->Draw(TString("sqrt(RoadCandidateArray[].RoadSimHitArray[].road.v22_l)>>")+road_y->GetName(),where && globalCut);
  
  c->Clear();
  c->Divide(2,1);
  c->cd(1)->SetGrid();
  road_x->Draw();

  c->cd(2)->SetGrid();
  road_y->Draw();

 return c;
}

TCanvas * doResidual(TCut where, TString tag,int Nbinx, double xmin,double xmax,int Nbiny,double ymin,double ymax)
{
  cout<<"doResidual: "<<where<<" "<<tag<<endl;
  gStyle->SetOptStat(true);
  gStyle->SetOptFit(true);
  
  TH1F * res_x = new TH1F(tag+"_residual_x","local x residuals at "+tag,Nbinx,xmin,xmax);
  res_x->SetXTitle("x^{local}_{road}-x^{local}_{hit} [cm]");
  TH1F * res_y = new TH1F(tag+"_residual_y","local y residuals at "+tag,Nbiny,ymin,ymax);
  res_y->SetXTitle("y^{local}_{road}-y^{local}_{hit} [cm]");
  
  TCanvas * c = new TCanvas(tag+"_residual","residuals at "+tag+"",2*xblock,1*yblock);
  data->Draw(TString("RoadCandidateArray[].RoadSimHitArray[].road.x_l-RoadCandidateArray[].RoadSimHitArray[].simhit.x_l>>")+res_x->GetName(),where && globalCut);
  data->Draw(TString("RoadCandidateArray[].RoadSimHitArray[].road.y_l-RoadCandidateArray[].RoadSimHitArray[].simhit.y_l>>")+res_y->GetName(),where && globalCut);
  
  c->Clear();
  c->Divide(2,1);
  c->cd(1)->SetGrid();
  res_x->Draw();
  res_x->Fit("gaus","q+");

  c->cd(2)->SetGrid();
  res_y->Draw();
  res_y->Fit("gaus","q+");

 return c;
}

TCanvas * doPull(TCut where, TString tag,int Nbinx, double xmin,double xmax,int Nbiny,double ymin,double ymax,int Nbinc,double cmin,double cmax)
{
  cout<<"doPull: "<<where<<" "<<tag<<endl;
  gStyle->SetOptStat(true);
  gStyle->SetOptFit(true);
  
  TH1F * res_x = new TH1F(tag+"_pull_x","local x pulls at "+tag,Nbinx,xmin,xmax);
  res_x->SetXTitle("(x^{local}_{road}-x^{local}_{hit})/#sigma^{road}_{x^{local}} [100%]");
  TH1F * res_y = new TH1F(tag+"_pull_y","local y pulls at "+tag,Nbiny,ymin,ymax);
  res_y->SetXTitle("(y^{local}_{road}-y^{local}_{hit})/#sigma^{road}_{y^{local}} [100%]");
  TH1F * chi2 = new TH1F(tag+"_chi2_road","road chi2 at "+tag,Nbinc,cmin,cmax);
  chi2->SetXTitle("#chi^{2}_{road}");

  TCanvas * c = new TCanvas(tag+"_pull","pulls at "+tag+"",3*xblock,1*yblock);

  data->Draw(TString("(RoadCandidateArray[].RoadSimHitArray[].road.x_l-RoadCandidateArray[].RoadSimHitArray[].simhit.x_l)/sqrt(RoadCandidateArray[].RoadSimHitArray[].road.v11_l)>>")+res_x->GetName(),where && globalCut);
  data->Draw(TString("(RoadCandidateArray[].RoadSimHitArray[].road.y_l-RoadCandidateArray[].RoadSimHitArray[].simhit.y_l)/sqrt(RoadCandidateArray[].RoadSimHitArray[].road.v22_l)>>")+res_y->GetName(),where && globalCut);

  data->Draw(TString("RoadCandidateArray[].RoadSimHitArray[].rechit.chi2_road>>")+chi2->GetName(),where && globalCut);

  c->Clear();
  c->Divide(3,1);
  c->cd(1)->SetGrid();
  res_x->Draw();
  res_x->Fit("gaus","q+");

  c->cd(2)->SetGrid();
  res_y->Draw();
  res_y->Fit("gaus","q+");

  c->cd(3)->SetGrid();
  c->cd(3)->SetLogy();
  chi2->Draw();
  //  chi2->Fit("gaus","q+");
  
  return c;
}


TList * doRsizeResiduPull(TCut where, TString tag,int Nbinx, double xmin,double xmax,int Nbiny,double ymin,double ymax,int Nbinc,double cmin,double cmax)
{ TList * canvas = new TList();
  canvas->AddLast(doRoadSize(where,tag,Nbinx,xmin,xmax,Nbiny,ymin,ymax));
  canvas->AddLast(doResidual(where,tag,Nbinx,xmin,xmax,Nbiny,ymin,ymax));
  canvas->AddLast(doPull(where,tag,Nbinx,xmin,xmax,Nbiny,ymin,ymax,Nbinc,cmin,cmax));
  return canvas;}

TList * doAllLayer_RsizeResiduPull(int Nbinx, double xmin,double xmax,int Nbiny,double ymin,double ymax,int Nbinc,double cmin,double cmax,int start=0,int end=22,TCut additionalcut="hasSim")
{ TList * canvas = new TList();
  TCanvas * m =  new TCanvas("geometry","geometry of tracker",(int)3*xblock,(int)1.2*xblock);
  TList * mapgen = new TList();
  TH2F * genmap = new TH2F("genmap","genmap",300,0,300,120,0,120);
  TLegend * mleg = new TLegend(0.85,0.2,1,1);
  genmap->SetStats(false);
  m->cd();
  data->Draw("RoadCandidateArray[].RoadSimHitArray[].rechit.r():abs(RoadCandidateArray[].RoadSimHitArray[].rechit.z)>>genmap",globalCut);  

  if (end<start) end=start+1;
  for (int il=start;il!=end;il++)
    {
      TCut where = layercut(il);

      TString tag=lnames[il];
      canvas->AddAll(doRsizeResiduPull(where&& additionalcut,tag, Nbinx,xmin,xmax, Nbiny,ymin,ymax, Nbinc,cmin,cmax));

      m->cd();
      TH2F * h2=new TH2F(tag+"_map",tag+" map",90,0,300,36,0,120);
      mleg->AddEntry(h2,tag,"f");
      h2->SetMarkerColor(il%8+2);
      h2->SetFillColor(il%8+2);
      h2->SetStats(false);
      data->Draw(TString("RoadCandidateArray[].RoadSimHitArray[].rechit.r():abs(RoadCandidateArray[].RoadSimHitArray[].rechit.z)>>")+h2->GetName(),where && additionalcut && globalCut,"");
      mapgen->AddLast(h2);
    }
  m->cd();
  m->Clear();
  TListIter it(mapgen);
  TH2F * h2;
  genmap->Draw();
  while ((bool)(h2=(TH2F*)it()))
    h2->Draw("same box");
  genmap->Draw("same");
  canvas->AddLast(m);    
  mleg->Draw();
  return canvas;
}


//parametristion of the seed error matrix as a function of pT, phi, eta of the seed (can be changed to sim)
TFile * _errorParametrisationFile= NULL;
TFile * openParametrisationFile(char * filename)
{  TDirectory * oldDir = gDirectory; _errorParametrisationFile= TFile::Open(filename,"RECREATE"); oldDir->cd();return _errorParametrisationFile;}


TCanvas * doOneParametrisation(TString tag, TString var, TString label, int Nb, double min, double max)
{
  cout<<"doOneParametrisation: "<<tag<<" "<<var<<" "<<label<<endl;
  TString opt="";

  TList createdhere;

  TCanvas * c = new TCanvas (tag+"_parametrization","parametrisation of "+label,2*xblock,3*yblock);
  TString who="seed";

  int NPt=10;
  double minPt=1;
  double maxPt=200;
  
  TH2D * h2_Pt = new TH2D("h2_"+tag+"_Pt",label+" as a function of muon p_{T}",NPt,minPt,maxPt,Nb,min,max);
  h2_Pt->SetXTitle("moun p^{sta}_{T} [GeV]");
  h2_Pt->SetYTitle(label);
  createdhere.AddLast(h2_Pt);
  
  int NEta=20;
  double minEta=0;
  double maxEta=2.5;
  TH2D * h2_Eta = new TH2D("h2_"+tag+"_Eta",label+" as a function of muon |#eta|",NEta,minEta,maxEta,Nb,min,max);
  h2_Eta->SetXTitle("moun |#eta^{sta}|");
  h2_Eta->SetYTitle(label);
  createdhere.AddLast(h2_Eta);
  
  int NPhi=20;
  double minPhi=-TMath::Pi();
  double maxPhi=TMath::Pi();
  TH2D * h2_Phi = new TH2D("h2_"+tag+"_Phi",label+" as a function of muon #varphi",NPhi,minPhi,maxPhi,Nb,min,max);
  h2_Phi->SetXTitle("moun #varphi^{sta}");
  h2_Phi->SetYTitle(label);
  createdhere.AddLast(h2_Phi);
  
  TProfile * pf_Pt = h2_Pt->ProfileX("pf_"+tag+"_Pt",-1,-1,opt);
  pf_Pt->SetTitle(label+" as a function of muon p_{T}");
  pf_Pt->SetXTitle("moun p_{T} [GeV]");
  pf_Pt->SetYTitle(label);
  pf_Pt->SetMarkerStyle(7);
  pf_Pt->SetMarkerColor(2);
  createdhere.AddLast(pf_Pt);

  TProfile * pf_Eta = h2_Eta->ProfileX("pf_"+tag+"_Eta",-1,-1,opt);
  pf_Eta->SetTitle(label+" as a function of muon |#eta|");
  pf_Eta->SetXTitle("moun |#eta|");
  pf_Eta->SetYTitle(label);
  pf_Eta->SetMarkerStyle(7);
  pf_Eta->SetMarkerColor(2);
  createdhere.AddLast(pf_Eta);

  TProfile * pf_Phi = h2_Phi->ProfileX("pf_"+tag+"_Phi",-1,-1,opt);
  pf_Phi->SetTitle(label+" as a function of muon #varphi");
  pf_Phi->SetXTitle("moun #varphi");
  pf_Phi->SetYTitle(label);
  pf_Phi->SetMarkerStyle(7);
  pf_Phi->SetMarkerColor(2);
  createdhere.AddLast(pf_Phi);

  TProfile3D * pf3 = new TProfile3D("pf3_"+tag,label+" in triple parametrisation",NPt,minPt,maxPt,NEta,minEta,maxEta,NPhi,minPhi,maxPhi,"S");
  pf3->SetXTitle("moun p_{T} [GeV]");
  pf3->SetYTitle("moun |#eta|");
  pf3->SetZTitle("moun #varphi");
  createdhere.AddLast(pf3);

  
  int N = data->GetEntries();
  data->GetEntry(0);

  TTreeFormula * cutFormula=NULL; 
  TTreeFormula *  ncand=NULL;
  TTreeFormula * phiFormula=NULL;
  TTreeFormula *etaFormula=NULL;
  TTreeFormula *pTFormula=NULL;
  TTreeFormula *varFormula=NULL;

  int fT=-999;
  /*
    TEventList * l = data->GetEventList();
    
    int Nl=l->GetN();
    for (int il=0;il!=Nl,++il)
    {
    int ie=l->GetEntry(il);*/

  TCut selection=globalCut;
  for (int ie=0;ie!=N;ie++)
    {

      if (ie%100==0) cout <<"doing entry :"<<ie<<"...."<<endl;
      data->GetEntry(ie);
      if (ie%100==0) cout<<"get entry done"<<endl;

      if (fT!=data->GetTreeNumber())
	{
	  if (ncand){
	    cout<<"new tree... deleting the formula"<<endl;
	    delete ncand;
	    delete phiFormula;
	    delete etaFormula;
	    delete pTFormula;
	    delete varFormula;
	    delete cutFormula;
	  }
	  cout<<"new tree... new formula"<<endl;
	  cutFormula = new TTreeFormula("cutFormula",selection.GetTitle(),data);
	  ncand= new TTreeFormula("ncand","Nseeds",data);
	  phiFormula= new TTreeFormula("phiFormula","RoadCandidateArray."+who+".phi()",data);
	  etaFormula= new TTreeFormula("etaFormula","RoadCandidateArray."+who+".eta(0)",data);
	  pTFormula= new TTreeFormula("pTFormula","RoadCandidateArray."+who+".pt()",data);
	  varFormula= new TTreeFormula("varFormula",var,data);
	  fT=data->GetTreeNumber();
	}

      int Ns =(int) ncand->EvalInstance();
      if (ie%100==0) cout<<"Ns done"<<endl;

      if (Ns<=0) continue;
      if (ie%100==0) cout <<"Ns="<<Ns<<endl;

      for (int is=0;is!=Ns;is++)
	{
	  bool sel = cutFormula->EvalInstance(is);
	  if (!sel) continue;

	  double eta = etaFormula->EvalInstance(is);
	  double phi =  phiFormula->EvalInstance(is);
	  double pT = pTFormula->EvalInstance(is);
	  double value = varFormula->EvalInstance(is);
	  pf3->Fill(pT,eta,phi,value,1);
	}
      if (ie%100==0) cout <<"done with entry :"<<ie<<"."<<endl;
    }
  

  TH1D* h_pT = (TH1D*) pf3->Project3D("xe");
  TH1D* h_eta = (TH1D*) pf3->Project3D("ye");
  TH1D* h_phi = (TH1D*) pf3->Project3D("ze");  
  createdhere.AddLast(h_pT);
  createdhere.AddLast(h_eta);
  createdhere.AddLast(h_phi);

  TH2D* h_pT_eta = (TH2D*) pf3->Project3D("xye");
  createdhere.AddLast(h_pT_eta);
    
  c->Clear();
  c->Divide(2,3);
  c->cd(1)->SetGrid();
  h2_Pt->Draw();
  pf_Pt->Draw("same");
  c->cd(2)->SetGrid();
  h2_Eta->Draw();
  pf_Eta->Draw("same");
  c->cd(3)->SetGrid();
  h2_Phi->Draw();
  pf_Phi->Draw("same");
  c->cd(4)->SetGrid();
  h_pT->Draw();
  c->cd(5)->SetGrid();
  h_eta->Draw();
  c->cd(6)->SetGrid();
  h_phi->Draw();
  if (_errorParametrisationFile)
    { 
      TDirectory * oldDir = gDirectory;
      _errorParametrisationFile->cd();
      TListIter it(&createdhere);
      TObject * o;
      while ((bool)(o=it()))
	o->Write();
      oldDir->cd();
    }
	
  return c;
}

TList * doParametrisation(bool doError = true,bool doCov =false)
{ TList * canvas =new TList();
  TString vars[5]={"#frac{q}{|p|}","#lambda","#varphi_{0}","X_{T}","Y_{T}"};
  int NbCov =20;
  double minCov=-1;
  double maxCov=1;
  int Nb[5]={100,100,100,100,100};
  double min[5]={0,0,0,0,0};
  double max[5]={0.02,0.008,0.01,0.2,0.5};

  int Nmax=5;
  for (int i=0;i!=Nmax;i++)
    {
      TString v11="RoadCandidateArray[].seed.v";
      v11+=(i+1);v11+=(i+1);

      for (int j=i;j!=5;j++)
	{
	  TString v22="RoadCandidateArray[].seed.v";
	  v22+=(j+1);v22+=(j+1);

     	  TString v12="RoadCandidateArray[].seed.v";
	  v12+=(i+1);
	  v12+=(j+1);

	  TString tag="V";
	  tag+=(i+1);
	  tag+=(j+1);
	  //	  cout <<" doing "<<tag<<endl;
	  TString label;
	  TString var;
	  if (i==j)
	    {
	      if (doError){
		//error plots
		label="#sigma^{2}_{"+vars[i]+"}";
		var="sqrt("+v11+")";
		canvas->AddLast(doOneParametrisation(tag,var,label,Nb[i],min[i],max[i]));
	      }
	     }
	  else
	    {
	      if (doCov){
		//covariance plots
		label="Cov("+vars[i]+","+vars[j]+")";
		var=v12+"/sqrt("+v11+"*"+v22+")";
		canvas->AddLast(doOneParametrisation(tag,var,label,NbCov,minCov,maxCov));
	      }
	    }
	}
    }

  if (_errorParametrisationFile)
    _errorParametrisationFile->Close();
  return canvas;
}


//plot the content of the designatted object






//plot the content of the designated object
//as a function of the required variable
//need primary(),secondary(),noise(),other(),nRoadHit to be available
TCanvas * doContent(TString which,TString tag,TString var,TString xlabel, int Nbin, double min,double max,TCut valid="hasBest")
{
  cout<<"doContent: "<<which<<" "<<tag<<" "<<var<<" "<<xlabel<<" "<<valid.GetTitle()<<endl;
  //road content as a function of the input
  gStyle->SetOptStat(false);
  gStyle->SetOptFit(true);
  
  TCanvas * c= new TCanvas(tag+"_Content","content versus "+xlabel,2*xblock,3*yblock);
  
  TProfile * h_Nhots_mc = new TProfile(tag+"h_Nhots_mc","# hits (MC)",Nbin,min,max,"S");
  TH2F * h2_Nhots_mc = new TH2F(tag+"h2_Nhots_mc","# hits (MC)",Nbin,min,max,100,-0.5,20.5);
  h2_Nhots_mc->SetXTitle(xlabel);
  h_Nhots_mc->SetXTitle(xlabel);
  h_Nhots_mc->SetMarkerStyle(7);
  h_Nhots_mc->SetMarkerColor(2);
  h_Nhots_mc->SetLineColor(4);

  data->Draw("RoadCandidateArray[].nRoadSimHit:"+var+">>"+h_Nhots_mc->GetName(),valid && globalCut);
  data->Draw("RoadCandidateArray[].nRoadSimHit:"+var+">>"+h2_Nhots_mc->GetName(),valid && globalCut);


  TProfile * h_Nhots = new TProfile(tag+"h_Nhots","# hits found",Nbin,min,max,"S");
  TH2F * h2_Nhots = new TH2F(tag+"h2_Nhots","# hits found",Nbin,min,max,100,-0.5,20.5);
  h2_Nhots->SetXTitle(xlabel);
  h_Nhots->SetXTitle(xlabel);
  h_Nhots->SetMarkerStyle(7);
  h_Nhots->SetMarkerColor(2);
  h_Nhots->SetLineColor(4);

  data->Draw(which+".nRoadHit:"+var+">>"+h_Nhots->GetName(),valid && globalCut);
  data->Draw(which+".nRoadHit:"+var+">>"+h2_Nhots->GetName(),valid && globalCut);

  TProfile * h_Ntrue = new TProfile(tag+"h_Ntrue","# true hits found",Nbin,min,max,"S");
  TH2F * h2_Ntrue = new TH2F(tag+"h2_Ntrue","# true hits found",Nbin,min,max,100,-0.5,20.5);
  h2_Ntrue->SetXTitle(xlabel);
  h_Ntrue->SetXTitle(xlabel);
  h_Ntrue->SetMarkerStyle(7);
  h_Ntrue->SetMarkerColor(2);
  h_Ntrue->SetLineColor(4);

  data->Draw(which+".primary():"+var+">>"+h_Ntrue->GetName(),valid && globalCut);
  data->Draw(which+".primary():"+var+">>"+h2_Ntrue->GetName(),valid && globalCut);

  TProfile * h_Nsecondary = new TProfile(tag+"h_Nsecondary","# secondary hits found",Nbin,min,max,"S");
  TH2F * h2_Nsecondary = new TH2F(tag+"h2_Nsecondary","# secondary hits found",Nbin,min,max,100,-0.5,5.5);
  h2_Nsecondary->SetXTitle(xlabel);
  h_Nsecondary->SetXTitle(xlabel);
  h_Nsecondary->SetMarkerStyle(7);
  h_Nsecondary->SetMarkerColor(2);
  h_Nsecondary->SetLineColor(4);
 
  data->Draw(which+".secondary():"+var+">>"+h_Nsecondary->GetName(),valid && globalCut);
  data->Draw(which+".secondary():"+var+">>"+h2_Nsecondary->GetName(),valid && globalCut);
    

  TProfile * h_Nothertrack = new TProfile(tag+"h_Nothertrack","# othertrack hits found",Nbin,min,max,"S");
  TH2F * h2_Nothertrack = new TH2F(tag+"h2_Nothertrack","# othertrack hits found",Nbin,min,max,100,-0.5,5.5);
  h2_Nothertrack->SetXTitle(xlabel);
  h_Nothertrack->SetXTitle(xlabel);
  h_Nothertrack->SetMarkerStyle(7);
  h_Nothertrack->SetMarkerColor(2);
  h_Nothertrack->SetLineColor(4);
 
  data->Draw(which+".other():"+var+">>"+h_Nothertrack->GetName(),valid && globalCut);
  data->Draw(which+".other():"+var+">>"+h2_Nothertrack->GetName(),valid && globalCut);
  
  TProfile * h_Nnoise = new TProfile(tag+"h_Nnoise","# noise hits found",Nbin,min,max);
  TH2F * h2_Nnoise = new TH2F(tag+"h2_Nnoise","# noise hits found",Nbin,min,max,100,-0.5,5.5);
  h2_Nnoise->SetXTitle(xlabel);
  h_Nnoise->SetXTitle(xlabel);
  h_Nnoise->SetMarkerStyle(7);
  h_Nnoise->SetMarkerColor(2);
  h_Nnoise->SetLineColor(4);

  data->Draw(which+".noise():"+var+">>"+h_Nnoise->GetName(),valid && globalCut);
  data->Draw(which+".noise():"+var+">>"+h2_Nnoise->GetName(),valid && globalCut);
  
  c->Clear();
  c->Divide(2,3);
  c->cd(1)->SetGrid();
  h2_Nhots->Draw();
  h_Nhots->Draw("p same");

  c->cd(2)->SetGrid();
  h2_Nhots_mc->Draw();
  h_Nhots_mc->Draw("p same");

  c->cd(3)->SetGrid();
  h2_Ntrue->Draw();
  h_Ntrue->Draw("p same");

  c->cd(4)->SetGrid();
  h2_Nsecondary->Draw();
  h_Nsecondary->Draw("p same");
  
  c->cd(5)->SetGrid();
  h2_Nothertrack->Draw();
  h_Nothertrack->Draw("p same");

  c->cd(6)->SetGrid();
  h2_Nnoise->Draw();
  h_Nnoise->Draw("p same");
  return c;
}
TCanvas * doRoadContent(TString tag,TString var,TString xlabel, int Nbin, double min,double max)
{  return doContent("RoadCandidateArray[]",tag+"_Road",var,xlabel,Nbin,min,max,"hasSim");}
TCanvas * doBestContent(TString tag,TString var,TString xlabel, int Nbin, double min,double max)
{  return doContent("RoadCandidateArray[].best",tag+"_BestTrack",var,xlabel,Nbin,min,max,"hasBest");}
TCanvas * doOtherContent(TString tag,TString var,TString xlabel, int Nbin, double min,double max)
{  return doContent("RoadCandidateArray[].CandidateArray[]",tag+"_OtherTrack",var,xlabel,Nbin,min,max,"hasSim");}

//plot the efficiency /purity of the designated object
//as a function of the required variable
//need primary(),nRoadHit to be available
TCanvas * doEffPurity(TString which,TString tag,TString var,TString xlabel, int Nbin, double min,double max,TString numerator=".primary()")
{
  TCut valid = "RoadCandidateArray[].hasBest==1";
  cout<<"doEffPurity: "<<which<<" "<<tag<<" "<<var<<" "<<xlabel<<" "<<numerator<<endl;
  
  gStyle->SetOptStat(false);
  gStyle->SetOptFit(true);

  double effMax=1.1;
  
  TF1 * road_p0 = new TF1("raod_p0","pol0");
  
  TCanvas * c= new TCanvas(tag+"PurityEfficiency",tag+" efficiency and purity versus "+xlabel,3*xblock,1*yblock);
  
  TProfile * hf_eff = new TProfile(tag+"hf_eff","Fraction of actual hits found",Nbin,min,max);
  hf_eff->SetXTitle(xlabel);
  hf_eff->SetYTitle("efficiency=#frac{# true hits found }{# hits (MC)}");
  hf_eff->SetLineColor(2);
  TH2F * h2_eff = new TH2F(tag+"h2_eff","Fraction of actual hits found",Nbin,min,max,100,0,effMax);
  h2_eff->SetXTitle(xlabel);
  h2_eff->SetYTitle("efficiency=#frac{# true hits found }{# hits (MC)}");
  h2_eff->SetMarkerStyle(6);

  TProfile * hf_true = new TProfile(tag+"hf_true","Fraction of hits found that are true",Nbin,min,max);
  hf_true->SetXTitle(xlabel);
  hf_true->SetYTitle("purity=#frac{# true hits found }{# hits found}");
  hf_true->SetLineColor(2);
  TH2F * h2_true = new TH2F(tag+"h2_true","Fraction of actual hits found",Nbin,min,max,100,0,effMax);
  h2_true->SetXTitle(xlabel);
  h2_true->SetYTitle("purity=#frac{# true hits found }{# hits (MC)}");
  h2_true->SetMarkerStyle(6);
 
  TH1F * h_ = new TH1F(tag+"Disribution",xlabel,Nbin,min,max);
  h_->SetXTitle(xlabel);


  data->Draw(which+numerator+"/RoadCandidateArray[].nRoadSimHit:"+var+">>"+hf_eff->GetName(),valid && globalCut);
  data->Draw(which+numerator+"/"+which+".nRoadHit:"+var+">>"+hf_true->GetName(),valid && globalCut);
  data->Draw(which+numerator+"/RoadCandidateArray[].nRoadSimHit:"+var+">>"+h2_eff->GetName(),valid && globalCut);
  data->Draw(which+numerator+"/"+which+".nRoadHit:"+var+">>"+h2_true->GetName(),valid && globalCut);
  data->Draw(var+">>"+h_->GetName(),valid && globalCut);

  c->Clear();
  c->Divide(3,1);
  
  c->cd(1)->SetGrid();
  c->cd(1)->SetLeftMargin(0.2);
  hf_eff->SetMinimum(0);  
  hf_eff->SetMaximum(effMax);  
  hf_eff->Draw("e1");
  hf_eff->GetYaxis()->SetTitleOffset(2);
  road_p0->SetParName(0,"efficiency");
  hf_eff->Fit(road_p0,"q","e1");
  h2_eff->Draw("same");

  c->cd(2)->SetGrid();
  c->cd(2)->SetLeftMargin(0.2);
  hf_true->SetMinimum(0);
  hf_true->SetMaximum(effMax);
  hf_true->Draw("e1");
  hf_true->GetYaxis()->SetTitleOffset(2);
  road_p0->SetParName(0,"purity");
  hf_true->Fit(road_p0,"q","e1");
  h2_true->Draw("same");  

  c->cd(3);
  h_->Draw("he");
  return c;
}

TCanvas * doRoadOwnEffPurity(TString tag,TString var,TString xlabel, int Nbin, double min,double max)
{ return doEffPurity("RoadCandidateArray[]",tag+"_own__Road",var,xlabel,Nbin,min,max,".own()");}
TCanvas * doBestOwnEffPurity(TString tag,TString var,TString xlabel, int Nbin, double min,double max)
{ return doEffPurity("RoadCandidateArray[].best",tag+"_own_BestTrack",var,xlabel,Nbin,min,max,".own()");}
TCanvas * doOtherOwnEffPurity(TString tag,TString var,TString xlabel, int Nbin, double min,double max)
{ return doEffPurity("RoadCandidateArray[].CandidateArray[]",tag+"_own_OtherTrack",var,xlabel,Nbin,min,max,".own()");}


TCanvas * doRoadEffPurity(TString tag,TString var,TString xlabel, int Nbin, double min,double max)
{ return doEffPurity("RoadCandidateArray[]",tag+"_Road",var,xlabel,Nbin,min,max);}
TCanvas * doBestEffPurity(TString tag,TString var,TString xlabel, int Nbin, double min,double max)
{ return doEffPurity("RoadCandidateArray[].best",tag+"_BestTrack",var,xlabel,Nbin,min,max);}
TCanvas * doOtherEffPurity(TString tag,TString var,TString xlabel, int Nbin, double min,double max)
{ return doEffPurity("RoadCandidateArray[].CandidateArray[]",tag+"_OtherTrack",var,xlabel,Nbin,min,max);}

TList * doRoadContentAndEfficiency(TString tag,TString var,TString xlabel, int Nbin, double min,double max)
{  TList * canvas =new TList();
  canvas->AddLast(doRoadContent(tag,var,xlabel,Nbin,min,max));
  canvas->AddLast(doRoadEffPurity(tag,var,xlabel,Nbin,min,max));
  canvas->AddLast(doRoadOwnEffPurity(tag,var,xlabel,Nbin,min,max));
  return canvas;}  
TList * doBestContentAndEfficiency(TString tag,TString var,TString xlabel, int Nbin, double min,double max)
{TList * canvas = new TList();
  canvas->AddLast(doBestContent(tag,var,xlabel,Nbin,min,max));
  canvas->AddLast(doBestEffPurity(tag,var,xlabel,Nbin,min,max));
  canvas->AddLast(doBestOwnEffPurity(tag,var,xlabel,Nbin,min,max));
  return canvas;}
TList * doOtherContentAndEfficiency(TString tag,TString var,TString xlabel, int Nbin, double min,double max)
{ TList * canvas = new TList();
  canvas->AddLast(doOtherContent(tag,var,xlabel,Nbin,min,max));
  canvas->AddLast(doOtherEffPurity(tag,var,xlabel,Nbin,min,max));
  canvas->AddLast(doOtherOwnEffPurity(tag,var,xlabel,Nbin,min,max));
  return canvas;}  

TList * doOneVar_ContentAndEfficiency(TString tag,TString var,TString xlabel, int Nbin, double min,double max)
{  TList * canvas = new TList();
  canvas->AddAll(doRoadContentAndEfficiency(tag,var,xlabel,Nbin,min,max));
  canvas->AddAll(doBestContentAndEfficiency(tag,var,xlabel,Nbin,min,max));
  canvas->AddAll(doOtherContentAndEfficiency(tag,var,xlabel,Nbin,min,max));
  return canvas;}


//do the resolution between the two objects
//as a function of the required variable
//pt() required for both
TCanvas * doRes(TString which,TString w_label,TString which_ref,TString w_r_label,TString res_var, TString r_v_label,TString unit,TString tag,TString var,TString xlabel, int Nbin, double min,double max,TCut selection="1")
{
  cout<<"doRes: "<<which<<" "<<w_label<<" "<<which_ref<<" "<<w_r_label<<" "<<res_var<<" "<<r_v_label<<" "<<unit<<" "<<tag<<" "<<var<<" "<<xlabel<<" "<<selection.GetTitle()<<endl;
  
  gStyle->SetOptStat(false);
  gStyle->SetOptFit(false);
  
  TString pf_opt="";

  TCanvas * c= new TCanvas(tag+"_resolution",r_v_label+" resolution versus "+xlabel,2*xblock,2*yblock);
  
  TH2D * h2_rel = new TH2D(tag+"_h2_resrel",r_v_label+" rel.diff. versus "+xlabel,Nbin,min,max,300,-1,3);
  h2_rel->SetXTitle(xlabel);
  h2_rel->SetYTitle("("+w_label+"-"+w_r_label+")/"+ w_r_label);

  TH2D * h2 = new TH2D(tag+"_h2_res",r_v_label+" diff. versus "+xlabel,Nbin,min,max,300,-500,500);
  h2->SetXTitle(xlabel);
  h2->SetYTitle("("+ w_label+"-"+ w_r_label+") "+unit);
  
  TProfile * pf_rel = new TProfile(tag+"_pf_rel",r_v_label+" rel.diff. versus "+xlabel,Nbin,min,max,-1,3,pf_opt);
  pf_rel->SetXTitle(xlabel);
  pf_rel->SetYTitle("("+w_label+"-"+w_r_label+")/"+ w_r_label+" [100%]");

  TProfile * pf = new TProfile(tag+"_pf",r_v_label+" diff. versus "+xlabel,Nbin,min,max,-500,500,pf_opt);
  pf->SetXTitle(xlabel);
  pf->SetYTitle("("+ w_label+"-"+ w_r_label+") [GeV]");


  pf->SetMarkerStyle(6);
  pf->SetMarkerColor(2);
  pf->SetLineColor(4);

  pf_rel->SetMarkerStyle(6);
  pf_rel->SetMarkerColor(2);
  pf_rel->SetLineColor(4);
  

  data->Draw("("+which+"."+res_var+"-"+which_ref+"."+res_var+"):"+var+">>"+pf->GetName(),selection && globalCut,"prof"+pf_opt);
  data->Draw("("+which+"."+res_var+"-"+which_ref+"."+res_var+"):"+var+">>"+h2->GetName(),selection && globalCut);

  data->Draw("("+which+"."+res_var+"-"+which_ref+"."+res_var+")/"+which_ref+"."+res_var+":"+var+">>"+pf_rel->GetName(),selection && globalCut,"prof"+pf_opt);
  data->Draw("("+which+"."+res_var+"-"+which_ref+"."+res_var+")/"+which_ref+"."+res_var+":"+var+">>"+h2_rel->GetName(),selection && globalCut);

  
  //  TH1D * pf_rel_rms = RMS(pf_rel);
  TH1D * pf_rel_rms = CoreWidth(h2_rel);
  pf_rel_rms->SetTitle(r_v_label+" resolution versus "+xlabel);
  pf_rel_rms->SetXTitle(xlabel);
  pf_rel_rms->SetYTitle("#sigma_{("+w_label+"-"+w_r_label+")/"+ w_r_label+"} [100%]");
  pf_rel_rms->SetMarkerStyle(6);

  pf_rel_rms->SetMinimum(0.00001);
  pf_rel_rms->SetMaximum(1);


  //  TH1D * pf_rms = RMS(pf);
  TH1D * pf_rms = CoreWidth(h2,100);
  pf_rms->SetTitle(r_v_label+" abs. resolution versus "+xlabel);
  pf_rms->SetXTitle(xlabel);
  pf_rms->SetYTitle("#sigma_{("+ w_label+"-"+ w_r_label+")} "+unit);
  pf_rms->SetMarkerStyle(6);

  pf_rms->SetMinimum(0);
  pf_rms->SetMaximum(100);


  c->Clear();
  c->Divide(2,2);
  c->cd(1)->SetGrid();
  h2->Draw();
  pf->Draw("same");

  c->cd(2)->SetGrid();
  pf_rms->Draw();


  c->cd(3)->SetGrid();
  h2_rel->Draw();
  pf_rel->Draw("same");

  c->cd(4)->SetGrid();
  pf_rel_rms->Draw();

  return c;
}

TCanvas * doPtRes(TString which,TString w_label,TString which_ref,TString w_r_label,TString tag,TString var,TString xlabel, int Nbin, double min,double max,TCut selection="1")
{  return doRes(which,w_label,which_ref,w_r_label,"pt()","p_{T}","[GeV]",tag+"_pT",var,xlabel,Nbin,min,max,selection);}
TCanvas * doPRes(TString which,TString w_label,TString which_ref,TString w_r_label,TString tag,TString var,TString xlabel, int Nbin, double min,double max,TCut selection="1")
{  return doRes(which,w_label,which_ref,w_r_label,"p()","momentum","[GeV]",tag+"_p",var,xlabel,Nbin,min,max,selection);}
TCanvas * doEtaRes(TString which,TString w_label,TString which_ref,TString w_r_label,TString tag,TString var,TString xlabel, int Nbin, double min,double max,TCut selection="1")
{  return doRes(which,w_label,which_ref,w_r_label,"eta()","eta","",tag+"_eta",var,xlabel,Nbin,min,max,selection);}


TCanvas * doPtResBest(TString tag,TString var,TString xlabel, int Nbin, double min,double max)
{ return doPtRes("RoadCandidateArray[].best.updated","p_{T}^{best}","RoadCandidateArray[].sim","p_{T}^{true}",tag+"_best_true",var,xlabel,Nbin,min,max,"hasBest");}
TCanvas * doPtResSTA(TString tag,TString var,TString xlabel, int Nbin, double min,double max)
{ return doPtRes("RoadCandidateArray[].seed","p_{T}^{seed}","RoadCandidateArray[].sim","p_{T}^{true}",tag+"_sta_true",var,xlabel,Nbin,min,max);}
TCanvas * doPtRespca(TString tag,TString var,TString xlabel, int Nbin, double min,double max)
{ return doPtRes("RoadCandidateArray[].pca","p_{T}^{PCA}","RoadCandidateArray[].sim","p_{T}^{true}",tag+"_cpa_true",var,xlabel,Nbin,min,max,"hasPca");}
TCanvas * doPtResCpca(TString tag,TString var,TString xlabel, int Nbin, double min,double max)
{ return doPtRes("RoadCandidateArray[].Cpca","p_{T}^{cPCA}","RoadCandidateArray[].sim","p_{T}^{true}",tag+"_Ccpa_true",var,xlabel,Nbin,min,max,"hasCpca");}
TList * doOneVar_PtRes(TString tag,TString var,TString xlabel, int Nbin, double min,double max)
{ TList * canvas = new TList();
  canvas->AddLast(doPtResSTA(tag,var,xlabel,Nbin,min,max));
  canvas->AddLast(doPtRespca(tag,var,xlabel,Nbin,min,max));
  canvas->AddLast(doPtResCpca(tag,var,xlabel,Nbin,min,max));
  canvas->AddLast(doPtResBest(tag,var,xlabel,Nbin,min,max));
  return canvas;}

TCanvas * doPResBest(TString tag,TString var,TString xlabel, int Nbin, double min,double max)
{ return doPRes("RoadCandidateArray[].best.updated","p^{best}","RoadCandidateArray[].sim","p^{true}",tag+"_best_true",var,xlabel,Nbin,min,max,"hasBest");}
TCanvas * doPResSTA(TString tag,TString var,TString xlabel, int Nbin, double min,double max)
{ return doPRes("RoadCandidateArray[].seed","p^{seed}","RoadCandidateArray[].sim","p^{true}",tag+"_sta_true",var,xlabel,Nbin,min,max);}
TCanvas * doPRespca(TString tag,TString var,TString xlabel, int Nbin, double min,double max)
{ return doPRes("RoadCandidateArray[].pca","p^{PCA}","RoadCandidateArray[].sim","p^{true}",tag+"_cpa_true",var,xlabel,Nbin,min,max,"hasPca");}
TCanvas * doPResCpca(TString tag,TString var,TString xlabel, int Nbin, double min,double max)
{ return doPRes("RoadCandidateArray[].Cpca","p^{Cpca}","RoadCandidateArray[].sim","p^{true}",tag+"_Ccpa_true",var,xlabel,Nbin,min,max,"hasCpca");}
TList * doOneVar_PRes(TString tag,TString var,TString xlabel, int Nbin, double min,double max)
{ TList * canvas = new TList();
  canvas->AddLast(doPResSTA(tag,var,xlabel,Nbin,min,max));
  canvas->AddLast(doPRespca(tag,var,xlabel,Nbin,min,max));
  canvas->AddLast(doPResCpca(tag,var,xlabel,Nbin,min,max));
  canvas->AddLast(doPResBest(tag,var,xlabel,Nbin,min,max));
  return canvas;}

TCanvas * doEtaResBest(TString tag,TString var,TString xlabel, int Nbin, double min,double max)
{ return doEtaRes("RoadCandidateArray[].best.updated","#eta^{best}","RoadCandidateArray[].sim","#eta^{true}",tag+"_best_true",var,xlabel,Nbin,min,max,"hasBest");}
TCanvas * doEtaResSTA(TString tag,TString var,TString xlabel, int Nbin, double min,double max)
{ return doEtaRes("RoadCandidateArray[].seed","#eta^{seed}","RoadCandidateArray[].sim","#eta^{true}",tag+"_sta_true",var,xlabel,Nbin,min,max);}
TCanvas * doEtaRespca(TString tag,TString var,TString xlabel, int Nbin, double min,double max)
{ return doEtaRes("RoadCandidateArray[].pca","#eta^{PCA}","RoadCandidateArray[].sim","#eta^{true}",tag+"_cpa_true",var,xlabel,Nbin,min,max,"hasPca");}
TCanvas * doEtaResCpca(TString tag,TString var,TString xlabel, int Nbin, double min,double max)
{ return doEtaRes("RoadCandidateArray[].Cpca","#eta^{cPCA}","RoadCandidateArray[].sim","#eta^{true}",tag+"_Ccpa_true",var,xlabel,Nbin,min,max,"hasCpca");}
TList *doOneVar_EtaRes(TString tag,TString var,TString xlabel, int Nbin, double min,double max)
{ TList * canvas = new TList();
  canvas->AddLast(doEtaResSTA(tag,var,xlabel,Nbin,min,max));
  canvas->AddLast(doEtaRespca(tag,var,xlabel,Nbin,min,max));
  canvas->AddLast(doEtaResCpca(tag,var,xlabel,Nbin,min,max));
  canvas->AddLast(doEtaResBest(tag,var,xlabel,Nbin,min,max));
  return canvas;}

TList * doOneVar_Res(TString tag,TString var,TString xlabel, int Nbin, double min,double max)
{ TList * canvas = new TList();
  canvas->AddAll(doOneVar_PtRes(tag,var,xlabel,Nbin,min,max));
  //  canvas->AddAll(doOneVar_PRes(tag,var,xlabel,Nbin,min,max));
  //  canvas->AddAll(doOneVar_EtaRes(tag,var,xlabel,Nbin,min,max));
  return canvas;}


//make a comparison of the given variable
TCanvas * doComp(TString tag, TString var , TString xlabel, int Nbin,double min,double max)
{
  cout<<"doComp: "<<tag<<" "<<var<<" "<<xlabel<<endl;

  gStyle->SetOptStat(false);
  gStyle->SetOptFit(true);
  TCanvas * c = new TCanvas(tag+"_comparison","comparison of "+xlabel,xblock,yblock);

  TH1F * h_gen = new TH1F(tag+"_h_gen",xlabel+" for generated muons",Nbin,min,max);
  h_gen->SetXTitle(xlabel);
  h_gen->SetLineWidth(2);
  h_gen->SetLineColor(1);
  TH1F * h_sta = new TH1F(tag+"_h_sta",xlabel+" for STA muons",Nbin,min,max);
  h_sta->SetXTitle(xlabel);
  h_sta->SetLineWidth(2);
  h_sta->SetLineColor(2);
  TH1F * h_pca = new TH1F(tag+"_h_pca",xlabel+" for STA muons propagated to PCA",Nbin,min,max);
  h_pca->SetXTitle(xlabel);
  h_pca->SetLineWidth(2);
  h_pca->SetLineColor(3);
  TH1F * h_cpca = new TH1F(tag+"_h_cpca",xlabel+" for STA muons combined with IP",Nbin,min,max);
  h_cpca->SetXTitle(xlabel);  
  h_cpca->SetLineWidth(2);
  h_cpca->SetLineColor(4);
  TH1F * h_best = new TH1F(tag+"_h_best",xlabel+" for STA muons combined with IP",Nbin,min,max);
  h_best->SetXTitle(xlabel);  
  h_best->SetLineWidth(2);
  h_best->SetLineColor(7);

  cout<<"draw sim"<<endl;
  data->GetEntry(0);
  //  cout << data<<endl;
  data->Draw("RoadCandidateArray[].sim."+var+">>"+h_gen->GetName(),"hasSim" && globalCut);
  //  cout<<"draw sta"<<endl;
  data->Draw("RoadCandidateArray[].seed."+var+">>"+h_sta->GetName(),"hasSeed" && globalCut);
  //  cout<<"draw spca"<<endl;
  data->Draw("RoadCandidateArray[].pca."+var+">>"+h_pca->GetName(),"hasPca" && globalCut);
  //  cout<<"draw cpca"<<endl;
  data->Draw("RoadCandidateArray[].Cpca."+var+">>"+h_cpca->GetName(),"hasCpca" && globalCut);
  //  cout<<"draw best"<<endl;
  data->Draw("RoadCandidateArray[].best.updated."+var+">>"+h_best->GetName(),"hasBest" && globalCut);
  
  cout<<"clear"<<endl;
  c->Clear();
  TLegend * leg = new TLegend(0.7,0.7,0.95,0.95);
  leg->AddEntry(h_gen,"generated","l");
  leg->AddEntry(h_sta,"seed","l");
  leg->AddEntry(h_pca,"at PCA","l");
  leg->AddEntry(h_cpca,"at PCA(+IP)","l");
  leg->AddEntry(h_best,"best track","l");
  
  cout<<"redraw"<<endl;
  h_sta->Draw();
  h_pca->Draw("same");
  h_cpca->Draw("same");
  h_gen->Draw("same");
  h_best->Draw("same");
  leg->Draw();
  return c;
}



//algorithm step efficiency as a function the given variable
TCanvas * doAlgoStepEff(TCut num,TString num_label,TCut denom,TString denom_label,TString tag,TString var, TString xlabel,int Nbin,double min,double max)
{
  cout<<"doAlgoStepEff: "<<num.GetTitle()<<" "<<num_label<<" "<<denom.GetTitle()<<" "<<denom_label<<" "<<tag<<" "<<var<<" "<<xlabel<<endl;

  gStyle->SetOptStat(false);
  gStyle->SetOptFit(true);

  TH1D * h_denom = new TH1D(tag+"_hdenom",TString("seed ")+xlabel+" distribution",Nbin,min,max);
  TH1D * h_num = new TH1D(tag+"_hnum",TString("track candidate ")+xlabel+" distribution",Nbin,min,max);
  TH1D * h_bad = new TH1D(tag+"_hbad",TString("bad seed ")+xlabel+" distribution",Nbin,min,max);
  TH1D * h_ratio = new TH1D(tag+"_hratio","efficiency versus "+xlabel,Nbin,min,max);
  h_ratio->Sumw2();h_num->Sumw2();h_denom->Sumw2();

  h_denom->SetXTitle(xlabel);
  h_num->SetXTitle(xlabel);
  h_bad->SetXTitle(xlabel);
  h_ratio->SetXTitle(xlabel);

  h_denom->SetYTitle("Nb. of event");
  h_num->SetYTitle("Nb. of event");
  h_bad->SetYTitle("Nb. of event");
  h_ratio->SetYTitle("efficiency=("+num_label+")/("+denom_label+")");
  
  h_denom->SetLineWidth(2);
  h_num->SetLineWidth(2);
  h_bad->SetLineWidth(2);
  h_bad->SetLineColor(2);
  h_bad->SetLineStyle(3);
  h_num->SetLineStyle(2);
  
  TCanvas * c = new TCanvas(tag+"_AlgoEfficiency","efficiency versus "+xlabel,2*xblock,1*yblock);

  data->Draw(var+">>"+h_denom->GetName(),denom && globalCut);
  data->Draw(var+">>"+h_num->GetName(),num && globalCut);
  data->Draw(var+">>"+h_bad->GetName(),denom && !num && globalCut);
  data->Draw(var+">>"+h_ratio->GetName(),num && globalCut);
  h_ratio->Divide(h_denom);
  
  c->Clear();
  c->Divide(2,1);
  c->cd(1);
  TLegend * leg = new TLegend(0.7,0.8,0.95,0.95);
  leg->AddEntry(h_denom,"# of "+denom_label,"l");
  leg->AddEntry(h_num,"# of "+num_label,"l");
  leg->AddEntry(h_bad,"# of bad "+num_label,"l");

  h_denom->SetMinimum(0.1);
  h_denom->Draw("hist");
  h_num->Draw("same hist");
  h_bad->Draw("same hist");
  leg->Draw();

  c->cd(2)->SetGrid();
  h_ratio->SetMaximum(1.1);
  h_ratio->SetMinimum(0);
  
  h_ratio->Draw("e1");
  //  h_ratio->Fit(function,"","e1");
  return c;
}

TCanvas * doAlgoEff(TString tag,TString var, TString xlabel,int Nbin,double min,double max)
{ return doAlgoStepEff("hasBest && RoadCandidateArray[].best.nRoadHit>=3","best track","hasSim","seeds",tag+"_GlobalAlgoEff",var,xlabel,Nbin,min,max);}


TCanvas * doSeedToPCAEff(TString tag,TString var, TString xlabel,int Nbin,double min,double max)
{ return doAlgoStepEff("hasPca","valid PCA","hasSim","seeds",tag+"_SeedToPCA",var,xlabel,Nbin,min,max);}
TCanvas * doPCAToCPCAEff(TString tag,TString var, TString xlabel,int Nbin,double min,double max)
{ return doAlgoStepEff("hasCpca","valid combined PCA","hasPca","valid PCA",tag+"_PCAToCPCA",var,xlabel,Nbin,min,max);}
TCanvas * doCPCAToBestEff(TString tag,TString var, TString xlabel,int Nbin,double min,double max)
{ return doAlgoStepEff("hasBest","best track","hasCpca","valid combined PCA",tag+"_CPCAToBEST",var,xlabel,Nbin,min,max);}
TCanvas * doSimToCPCAEff(TString tag,TString var, TString xlabel,int Nbin,double min,double max)
{ return doAlgoStepEff("hasCpca","valid combined PCA","hasSim","seeds",tag+"_SimToCPCA",var,xlabel,Nbin,min,max);}



TList * doOneVar_Eff(TString tag,TString var, TString xlabel,int Nbin,double min,double max)
{ TList * canvas = new TList();
  canvas->AddLast(doAlgoEff(tag,var,xlabel,Nbin,min,max));
  //  canvas->AddLast(doSeedToPCAEff(tag,var,xlabel,Nbin,min,max));
  //  canvas->AddLast(doPCAToCPCAEff(tag,var,xlabel,Nbin,min,max));	
  canvas->AddLast(doCPCAToBestEff(tag,var,xlabel,Nbin,min,max));
  canvas->AddLast(doSimToCPCAEff(tag,var,xlabel,Nbin,min,max)); 
  return canvas;}



  const double rmax=130;
  const double zmax=300;
  const double xfraction = zmax/(rmax+zmax);
  const double yfraction = (rmax+zmax)/rmax;


/*
TCanvas * display(TString which,TString w_label,TString tag,TCut select="(1==1)")
{



  TCanvas * c = new TCanvas("display","display",yfraction*yblock,1*yblock);
  
  TH2D * h2_xy = new TH2D("h2_xy","x-y display",2*rmax,-rmax,rmax,2*rmax,-rmax,rmax);
  h2_xy->SetXTitle("x [cm]");  
  h2_xy->SetYTitle("y [cm]");

  TH2D * h2_rz = new TH2D("h2_rz","r-z display",2*zmax,-zmax,zmax,rmax,0,rmax);
  h2_rz->SetXTitle("z [cm]");
  h2_rz->SetYTitle("r [cm]");


  data->Draw(which+".y:"+which+".x>>"+h2_xy->GetName(),select && globalCut);
  data->Draw(which+".r():"+which+".z>>"+h2_rz->GetName(),select && globalCut);

  TString opt="";
  if (h2_rz->GetEntries()<20)
    {
      opt="box";
      h2_rz->SetFillColor(4);
      h2_xy->SetFillColor(4);
    }

  c->Clear();
  c->Divide(2,1);

  double xlow1,ylow1,xup1,yup1;
  double xlow2,ylow2,xup2,yup2;
  c->GetPad(1)->GetPadPar(xlow1,ylow1,xup1,yup1);
  c->GetPad(2)->GetPadPar(xlow2,ylow2,xup2,yup2);

  double xdiff=xlow2-xup1;

  xup1=xfraction;
  xlow2=xup1+xdiff;
  c->GetPad(1)->SetPad(xlow1,ylow1,xup1,yup1);
  c->GetPad(2)->SetPad(xlow2,ylow2,xup2,yup2);


  c->cd(1);
  h2_rz->Draw();

  c->cd(2);
  h2_xy->Draw();

return c;
}

*/



TCanvas * NRecHits(TString which, TString w_label, TString tag,TString var, TString xlabel,int Nbin,double min,double max,double Nmax,TCut is="1==1")
{
cout<<"NRecHits: "<<which<<" "<<w_label<<" "<<tag<<" "<<var<<" "<<xlabel<<endl;

  gStyle->SetOptStat(false);
  TCanvas * c = new TCanvas(tag+"nRecHit","number of "+w_label+" RecHits as a function of "+xlabel,xblock,yblock);

  TH2F * h2 = new TH2F(tag+"_h2d_nRecHit","# of "+w_label+" RecHits as a function of "+xlabel,Nbin,min,max,200,0,Nmax);
  TProfile * pf = new TProfile(tag+"_pf_nRecHit","# of "+w_label+" RecHits as a function of "+xlabel,Nbin,min,max,"S");

  pf->SetMarkerStyle(7);
  pf->SetMarkerColor(2);
  pf->SetLineColor(4);

  h2->SetXTitle(xlabel);
  pf->SetXTitle(xlabel);

  h2->SetYTitle("# of "+w_label+" RecHits");
  pf->SetYTitle("# of "+w_label+" RecHits");


  data->Draw(which+":"+var+">>"+h2->GetName(),is && globalCut);
  data->Draw(which+":"+var+">>"+pf->GetName(),is && globalCut);

  c->Clear();
  c->cd();
  h2->Draw();
  pf->Draw("same");
  return c;
}

TCanvas * nMuRecHits(TString tag,TString var, TString xlabel,int Nbin,double min,double max)
{  return NRecHits("RoadCandidateArray[].Nrechits","STA","STA_"+tag,var,xlabel,Nbin,min,max,50);}
TCanvas * nBestRecHits(TString tag ,TString var, TString xlabel,int Nbin,double min,double max)
{ return NRecHits("RoadCandidateArray[].best.nRoadHit","best track","best_"+tag,var,xlabel,Nbin,min,max,20,"hasBest");}
TCanvas * nOtherRecHits(TString tag ,TString var, TString xlabel,int Nbin,double min,double max)
{ return NRecHits("RoadCandidateArray[].CandidateArray[].nRoadHit","other tracks","other_"+tag,var,xlabel,Nbin,min,max,20,"hasBest");}
TCanvas * nBestTotalRecHits(TString tag ,TString var, TString xlabel,int Nbin,double min,double max)
{ return NRecHits("(RoadCandidateArray[].best.nRoadHit + RoadCandidateArray[].Nrechits)","total","total_"+tag,var,xlabel,Nbin,min,max,70,"hasBest");}

TList * doOneVar_nRecHits(TString tag,TString var, TString xlabel,int Nbin,double min,double max)
{ TList * canvas = new TList();
  canvas->Add(nMuRecHits(tag,var,xlabel,Nbin,min,max));
  canvas->Add(nBestRecHits(tag,var,xlabel,Nbin,min,max));
  canvas->Add(nOtherRecHits(tag,var,xlabel,Nbin,min,max));
  canvas->Add(nBestTotalRecHits(tag,var,xlabel,Nbin,min,max));
  return canvas;}

TList * doOneVar_EveryThing(TString tag,TString var, TString xlabel,int Nbin,double min,double max)
{ TList * canvas = new TList();
  canvas->AddAll(doOneVar_Eff(tag,var,xlabel,Nbin,min,max));
  canvas->AddAll(doOneVar_Res(tag,var,xlabel,Nbin,min,max));
  canvas->AddAll(doOneVar_ContentAndEfficiency(tag,var,xlabel,Nbin,min,max));
  canvas->AddAll(doOneVar_nRecHits(tag,var,xlabel,Nbin,min,max));
  return canvas;}


//how to write every histogram from a TDirectory into a directory
void plotAll(TDirectory * dir,TString prefix,int Nx,int Ny,TString what=".eps")
{
  TList * lc = dir->GetListOfKeys();
  TListIter iter(lc);
  TObject *o;
  int iP=1;
  int NperPage=Nx*Ny;
  TString cname="display_";
  int iC=1;
  TCanvas *c = new TCanvas(Form("%s%d",cname.Data(),iC),Form("display number %d",iC),Nx*xblock,Ny*yblock);
  while((bool)(o=iter()))
    {
      if (iP==NperPage+1)
	{ 
	  c->Print(prefix+"/"+c->GetName()+what);
	  iC++;
	  c = new TCanvas(Form("%s%d",cname.Data(),iC),Form("display number %d",iC),Nx*xblock,Ny*yblock);
	  iP=1;
	}

      if (o->InheritsFrom("TH1"))
	{
	  c->cd(iP);
	  o->Draw();
	}
    }
}

void saveAll(TString file, bool andCanvases=true){
  cout<<"save all objects to: "<<file<<endl;
  TDirectory * oldDir =gDirectory;
  TFile * f = TFile::Open(file,"RECREATE");
  f->cd();
  
  //  TList * lc = gROOT->GetListOfKeys();
  //it happens that they are not keys...
  TList * lc = gROOT->GetList();
  TListIter iter(lc);
  TObject * o;
  int n=0;
  while((bool)(o=iter())){o->Write();n++;}
  cout<<n<<" object(s) saved"<<endl;

  if (andCanvases){
    lc = (TList*)gROOT->GetListOfCanvases();
    TListIter iter(lc);
    TCanvas * c ;
    n=0;
    while((bool)(c=(TCanvas*)iter())) { c->Write();n++;}
    cout<<n<<" canvas(es) saved"<<endl;}
  f->Close();
  cout<<"going back to oldDir"<<endl;
  oldDir->cd();}


void plotAll(TDirectory * dir,TString prefix,TString what=".eps"){
  cout<<"print all canvases from: "<<dir->GetName()<<" to: "<<prefix<<"*"<<what<<endl;
  TList * lc = (TList*)dir->GetListOfKeys();
  TListIter iter(lc);
  TCanvas * c ;
  TKey * k;
  int n=0;
  if (gROOT->IsBatch()) {cout<<"root is in batch mode. cannot do that."; return;}
  while((bool)(k=(TKey*)iter())) { 
    if (TString(k->GetClassName())==TCanvas::Class_Name())
      { 
	c=(TCanvas*)dir->Get(k->GetName());
	cout <<"doing: "<<c->GetName()<<endl;
	double xw=c->GetVirtCanvas()->GetWw();
	double yw=c->GetVirtCanvas()->GetWh();

	c->Draw();
	c->SetWindowPosition(100,100);
	c->SetWindowSize(xw,yw);

	c->Print(prefix+c->GetName()+what);n++;}
  }
  cout<<n<<" canvas(es) printed"<<endl;}

//how to write every canvas in a directory
void plotAll(TString prefix,TString what=".eps"){
  cout<<"print all canvases to: "<<prefix<<"*"<<what<<endl;
  TList * lc = (TList*)gROOT->GetListOfCanvases();
  TListIter iter(lc);
  TCanvas * c ;
  int n=0;
  while((bool)(c=(TCanvas*)iter())) { c->Print(prefix+c->GetName()+what);n++;}
  cout<<n<<" canvas(es) printed"<<endl;}


enum cluster {SNU, UAF, LXPLUS};
//test purpose
void Monitor(cluster c=UAF,int Nf=10)
{
  if (c==UAF){
    //    gROOT->ProcessLine(".L /afs/fnal.gov/files/home/room2/vlimant/work/CMSSW_1_2_0_pre9/lib/slc3_ia32_gcc323/libRecoTrackerRoadSearchMuonSeededFinder.so");
    TString lib=gSystem->ExpandPathName("$CMSSW_BASE/lib/$SCRAM_ARCH/libRecoMuonL3MuonAnalyzer.so");
    std::cout<<" loading: "+lib<<endl;
    gROOT->ProcessLine(".L "+lib);

    //    open("monitor-test.root");
    //    openParametrisationFile("errorMatrix.root");
    //    openMany("l3muon-monitor.rootlist");

    //    open("Temp/run-ttb-1perModule.monitor.root");
    //    open("Temp/run-ttb-regular.monitor.root");
    //    open("Temp/run-ttb-dynamic.monitor.root");

    //    open("Temp/fromCERN/run-ttb-1perModule.monitor.root");
    //    open("Temp/fromCERN/run-ttb-regular.monitor.root");
    //    open("Temp/fromCERN/run-ttb-dynamic.monitor.root");

    //        openMany("/uscms/home/vlimant/work/CMSSW_1_2_0_pre9/src/RecoTracker/L3MuonTestEDProducer/test/TrajectorySeededTrackReco-monitor-many.rootlist",Nf);

  }
  if (c==SNU)    {
    gROOT->ProcessLine(".x Load.C");
    //    open("test/monitor-test.root");
    //    openParametrisationFile("test/errorMatrix.root");

    //open("/homes/vlimant/WorkDir/UAFmirror/root/1_2_0_pre9/TrajectorySeededTrackReco-monitor-0.root");
    //    openMany("../../L3MuonTestEDProducer/test/snu.rootlist",Nf);

    open("/homes/vlimant/WorkDir/UAFmirror/root/1_2_0_pre9/run-ttb-1perModule.monitor.root");
    //    open("/homes/vlimant/WorkDir/UAFmirror/root/1_2_0_pre9/run-ttb-dynamic.monitor.root");
    //    open("/homes/vlimant/WorkDir/UAFmirror/root/1_2_0_pre9/run-ttb-regular.monitor.root");
  }
  
  //  doMonitor(NULL,"/net/top/homes/vlimant/WorkDir/UAFmirror/root/1_2_0_pre9/plots/");

  //  doAllLayer_RsizeResiduPull(50,-10,10,50,-10,10,100,0,200,5,6);
  
  //  doOneVar_ContentAndEfficiency("nhitSeed","RoadCandidateArray[].Nrechits","# rechit on seed",71,0.5,70.5);
  
  //  doParametrisation(true,false);
  
}


//main macro
void doMonitor(char * rootfile="", TString dir="")
{
  if (rootfile)
    open(rootfile);

  gROOT->cd();
  //road, best track and other track composition
  //by layer
  //as a function of true_eta, true_pt
  //efficiency and purity of roads
  doOneVar_ContentAndEfficiency("true_pt","RoadCandidateArray[].sim.pt()","muon p_{T}^{true} [Gev]",100,0,200);
  doOneVar_ContentAndEfficiency("true_abs_eta","RoadCandidateArray[].sim.eta(0)","muon |#eta^{true}|",20,0,3);
  doOneVar_ContentAndEfficiency("true_eta","RoadCandidateArray[].sim.eta()","muon #eta^{true}",20,-3,3);
  doOneVar_ContentAndEfficiency("true_phi","RoadCandidateArray[].sim.phi()","muon #varphi^{true}",64,-TMath::Pi(),TMath::Pi());
  doOneVar_ContentAndEfficiency("nhitSeed","RoadCandidateArray[].Nrechits","# rechit on seed",71,-0.5,70.5);

  //pulls and residue 
  //chi2 of the rechit supposed to be on the seed
  doAllLayer_RsizeResiduPull(50,-10,10,50,-10,10,100,0,200);

  //sta resolution
  //as a function of eta, pt
  //track resolution
  //as a function of eta, pt
  doOneVar_Res("true_pt","RoadCandidateArray[].sim.pt()","muon p_{T}^{true} [Gev]",20,0,200);
  doOneVar_Res("true_abs_eta","RoadCandidateArray[].sim.eta(0)","muon |#eta^{true}|",20,0,2.5);
  doOneVar_Res("true_eta","RoadCandidateArray[].sim.eta()","muon #eta^{true}",20,-3,3);

  //comparison at different steps  
  doComp("pt","pt()","muon p_{T} [GeV]",100,0,500);
  doComp("p","p()","muon p [GeV]",100,0,750);
  doComp("px","px","muon p_{x} [GeV]",100,0,500);
  doComp("py","py","muon p_{y} [GeV]",100,0,500);
  doComp("pz","pz","muon p_{z} [GeV]",100,0,750);
  
  doComp("eta","eta()","muon #eta_{momentum}",100,-3,3);
  doComp("eta_p","eta_pos()","muon #eta_{position}",100,-3,3);
  doComp("abs_eta","eta(0)","muon |#eta|",50,0,3);
  doComp("phi","phi()","muon #varphi",50,-TMath::Pi(),TMath::Pi());
  doComp("rho","rho()","muon #rho [cm]",50,0,50);

  //algorithm efficiencies
  //as a function of whatever you want to show
  //step by step algo efficiencies (never done)
  //seed, hasPca, hasCpca, hasBest...
  doOneVar_Eff("true_pt","RoadCandidateArray[].sim.pt()","muon p_{T}^{true} [Gev]",20,0,200);
  doOneVar_Eff("true_abs_eta","RoadCandidateArray[].sim.eta(0)","muon |#eta^{true}|",20,0,2.5);
  doOneVar_Eff("true_eta","RoadCandidateArray[].sim.eta()","muon #eta^{true}",20,-3,3);
  doOneVar_Eff("true_phi","RoadCandidateArray[].sim.phi()","muon #varphi^{true}",40,-TMath::Pi(),TMath::Pi());


  doOneVar_Eff("Dpt_pca_true","RoadCandidateArray[].pca.pt()-RoadCandidateArray[].sim.pt()","muon (p_{T}^{PCA}-p_{T}^{true}) [GeV]",20,-200,200);
  doOneVar_Eff("DRelpt_pca_true","(RoadCandidateArray[].pca.pt()-RoadCandidateArray[].sim.pt())/RoadCandidateArray[].sim.pt()","muon (p_{T}^{PCA}-p_{T}^{true})/p_{T}^{true} [100%]",20,-2,2);

  doOneVar_Eff("Zpca","RoadCandidateArray[].pca.z","z_{PCA} [cm]",100,-50,50);
  doOneVar_Eff("Rpca","RoadCandidateArray[].pca.r()","r_{PCA} [cm]",20,0,20);
  doOneVar_Eff("Dpca","RoadCandidateArray[].pca.rho()","#rho_{PCA} [cm]",20,0,50);

  doOneVar_Eff("nhitSeed","RoadCandidateArray[].Nrechits","# rechit on seed",71,0.5,70.5);
  doOneVar_Eff("nseeds","Nseeds","# of seeds",8,-0.5,8.5);
  
  //error matrix parametrisation
  doParametrisation(true,true);

  saveAll(dir+"monitor-saved.root");

  plotAll(dir,".eps");

  if (!gROOT->IsBatch())
    plotAll(dir,".gif");
}

void LPC_14Dec06()
{

  doAllLayer_RsizeResiduPull(50,-10,10,50,-10,10,100,0,200);
  
  doOneVar_Eff("true_pt","RoadCandidateArray[].sim.pt()","muon p_{T}^{true} [Gev]",100,0,200);
  doOneVar_Eff("true_abs_eta","RoadCandidateArray[].sim.eta(0)","muon |#eta^{true}|",50,0,3);

  doOneVar_Eff("nhitSeed","RoadCandidateArray[].Nrechits","# rechit on seed",71,0.5,70.5);
  doOneVar_Eff("DRelpt_pca_true","(RoadCandidateArray[].pca.pt()-RoadCandidateArray[].sim.pt())/RoadCandidateArray[].sim.pt()","muon (p_{T}^{PCA}-p_{T}^{true})/p_{T}^{true} [100%]",20,-2,2);
  doOneVar_Eff("DReleta_pca_true","(RoadCandidateArray[].pca.eta()-RoadCandidateArray[].sim.eta())/RoadCandidateArray[].sim.eta()","muon (#eta^{PCA}-#eta^{true})/#eta^{true} [100%]",20,-0.1,0.1);
  
  //  doComp("pt","pt()","muon p_{T} [GeV]",100,0,500);
  
  //  doOneVar_Res("true_abs_eta","RoadCandidateArray[].sim.eta(0)","muon |#eta^{true}|",50,0,2.5);
  //  doOneVar_Res("true_pt","RoadCandidateArray[].sim.pt()","muon p_{T}^{true} [Gev]",100,0,200);

  plotAll("tracking_LPC_15Dec06/",".gif");
  plotAll("tracking_LPC_15Dec06/",".eps");

  //  plotAll("/net/top/homes/vlimant/Presentation/Tracking/tracking_LPC_15Dec06/",".gif");
  //  plotAll("/net/top/homes/vlimant/Presentation/Tracking/tracking_LPC_15Dec06/",".eps");

  doOneVar_nRecHits("true_abs_eta","RoadCandidateArray[].sim.eta(0)","muon |#eta^{true}|",50,0,3);
  doOneVar_nRecHits("true_pt","RoadCandidateArray[].sim.pt()","muon p_{T}^{true} [Gev]",100,0,200);

 }
