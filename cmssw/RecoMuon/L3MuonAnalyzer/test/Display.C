int _ie,_is;

TGraph * grab(TString gname,TString gtitle="")
{ TGraph * gr = new TGraph(*(TGraph*)gPad->GetPrimitive("Graph"));
  gr->SetName(gname);
  gr->SetTitle(gtitle);
  return gr;}

TGraph * generic(TString & who,TString &gname,int &mc,int &mst,float &ms,TString X,TString Y,TString label,bool select)
{ TString exp = "RoadCandidateArray["; exp+=_is; exp+="]."+who+Y+":RoadCandidateArray["; exp+=_is; exp+="]."+who+X;
  cout<<"drawing "<<exp<<endl;

  TCut sel="";
  if (!select)
    {TString ss="RoadCandidateArray[";ss+=_is;ss+="].RoadHitArray.matchcode<=1";sel=ss;}

  int N=data->Draw(exp,sel,"",1,_ie);
  TGraph * gr=NULL;
  if(N==0) {gr=new TGraph();gr->SetName(gname+label);}
  else {gr = grab(gname+label,gname+label);}

  gr->SetMarkerColor(mc);
  gr->SetMarkerStyle(mst);
  gr->SetMarkerSize(ms);
  return gr;}

TGraph * RZ(TString who,TString gname,int mc,int mst,float ms,bool select=true)
{return generic(who,gname,mc,mst,ms,".z",".r()","_rz",select);}
TGraph * XY(TString who,TString gname,int mc,int mst,float ms,bool select=true)
{return generic(who,gname,mc,mst,ms,".x",".y","_xy",select);}

EventDump * _ed = new EventDump();
TList * XYlist(TString gname,int mc,int mst,float ms){
  TList * result=new TList();
  data->SetBranchAddress("aDump",&_ed);
  data->GetEntry(_ie);
  RoadCandidate * rc = (RoadCandidate*) _ed->RoadCandidateArray->At(_is);
  std::vector<Candidate> cv = rc->CandidateArray;
  for (int ci=0;ci!=rc->nCandidate;++ci){
    TGraph * gr= new TGraph();
    gr->SetName(Form("%s_%d_xy",gname.Data(),ci));
    unsigned int ig=0;
    for (int rhi=0;rhi!=cv[ci].nRoadHit;++rhi){
      gr->SetPoint(ig++,cv[ci].RoadHitArray[rhi].rechit.x,cv[ci].RoadHitArray[rhi].rechit.y);}
    
    gr->SetMarkerColor(mc+ci);
    gr->SetMarkerStyle(mst);
    gr->SetMarkerSize(ms);
    
    result->Add(gr);}
  return result;
}

TList * RZlist(TString gname,int mc,int mst,float ms){
  TList * result=new TList();
  data->SetBranchAddress("aDump",&_ed);
  data->GetEntry(_ie);
  RoadCandidate * rc = (RoadCandidate*) _ed->RoadCandidateArray->At(_is);
  std::vector<Candidate> cv = rc->CandidateArray;
  for (unsigned int ci=0; ci != rc->nCandidate;++ci){
    TGraph * gr= new TGraph();
    gr->SetName(Form("%s_%d_rz",gname.Data(),ci));
    unsigned int ig=0;
    for (unsigned int rhi=0;rhi!=cv[ci].nRoadHit;++rhi){
      gr->SetPoint(ig++,cv[ci].RoadHitArray[rhi].rechit.z,cv[ci].RoadHitArray[rhi].rechit.r());}

    gr->SetMarkerColor(mc+ci);
    gr->SetMarkerStyle(mst);
    gr->SetMarkerSize(ms);
    
    result->Add(gr);}
  return result;
}


TGraph * XYmore(TString gname,int mc,int mst,float ms){
  unsigned int ig=0;
  TGraph * gr= new TGraph();
  gr->SetName(gname+"_xy");
  
  data->SetBranchAddress("aDump",&_ed);
  data->GetEntry(_ie);
  RoadCandidate * rc = (RoadCandidate*) _ed->RoadCandidateArray->At(_is);
  std::vector<Candidate> cv = rc->CandidateArray;
  for (unsigned int ci=0; ci != cv.size();++ci){
    for (unsigned int rhi=0;rhi!=cv[ci].RoadHitArray.size();++rhi){
      gr->SetPoint(ig++,cv[ci].RoadHitArray[rhi].rechit.x,cv[ci].RoadHitArray[rhi].rechit.y);
    }
  }

  gr->SetMarkerColor(mc); 
  gr->SetMarkerStyle(mst); 
  gr->SetMarkerSize(ms);    
  return gr;
}
TGraph * RZmore(TString gname,int mc,int mst,float ms){
  unsigned int ig=0;
  TGraph * gr= new TGraph();
  gr->SetName(gname+"_rz");

  data->SetBranchAddress("aDump",&_ed);
  data->GetEntry(_ie);
  RoadCandidate * rc = (RoadCandidate*) _ed->RoadCandidateArray->At(_is);
  std::vector<Candidate> cv = rc->CandidateArray;
  for (unsigned int ci=0; ci != cv.size();++ci){
    for (unsigned int rhi=0;rhi!=cv[ci].RoadHitArray.size();++rhi){
      gr->SetPoint(ig++,cv[ci].RoadHitArray[rhi].rechit.z,cv[ci].RoadHitArray[rhi].rechit.r());
    }}
  
  gr->SetMarkerColor(mc); 
  gr->SetMarkerStyle(mst); 
  gr->SetMarkerSize(ms);    
  return gr;
}


TArrow * generic(TString &who, int &ac, int &ast, TString X,TString Y, TString PX,TString PY)
{
  data->GetEntry(_ie);
  TTreeFormula px("px","RoadCandidateArray."+who+PX,data);
  TTreeFormula py("py","RoadCandidateArray."+who+PY,data);
  double pxv=px.EvalInstance(_is);
  double pyv=py.EvalInstance(_is);
  TTreeFormula x("x","RoadCandidateArray."+who+X,data);
  TTreeFormula y("y","RoadCandidateArray."+who+Y,data);
  double xv=x.EvalInstance(_is);
  double yv=y.EvalInstance(_is);
  TArrow * ar = new TArrow(xv,yv,xv+pxv,yv+pyv);
  ar->SetLineColor(ac);
  ar->SetLineStyle(ast);
  return ar;
  //  return new TArrow();
}

TArrow *arrowXY(TString who, int ac, int ast)
{return generic(who,ac,ast,".x",".y",".px",".py");}
TArrow *arrowRZ(TString who, int ac, int ast)
{return generic(who,ac,ast,".z",".r()",".pz",".pt()");}

TArrow * XYrotation(TArrow * ar, double phi)
{
  //  TArrow * r=new TArrow(*ar);
  TArrow * r=new TArrow();
  r->SetLineColor(ar->GetLineColor());
  r->SetLineStyle(ar->GetLineStyle());

  TVector2 v1(ar->GetX1(),ar->GetY1());
  TVector2 rv1=v1.Rotate(phi);
  r->SetX1(rv1.X());
  r->SetY1(rv1.Y());
  TVector2 v2(ar->GetX2(),ar->GetY2());
  TVector2 rv2=v2.Rotate(phi);
  r->SetX2(rv2.X());
  r->SetY2(rv2.Y());
  return r;
}

TGraph * XYrotation(const TGraph * gr ,double phi)
{
  TGraph * r = new TGraph(*gr);
  r->SetName(gr->GetName()+TString("_rotated"));
  r->SetTitle(gr->GetTitle()+TString(" (rotated)"));

  for (int i=0;i!=gr->GetN();++i)
    {double x,y;
      gr->GetPoint(i,x,y);
      TVector2 v(x,y);
      TVector2 rv=v.Rotate(phi);
      r->SetPoint(i,rv.X(),rv.Y());}
  return r;
}

double su2( TVector3 & v3, double phi)
{ double sx2=v3.X();
  double xy2=v3.Z();
  double V=v3.Y();

  double c=cos(phi);
  double c2=c*c;
  double s2=1-c2;
  double s=sqrt(s2);

  double result = c2/sx2 - 2*c*s*V/(sx2*sy2) + s2/sy2;
  return 1/result;}

TVector2 vector(TArrow * ar)
{return TVector2(ar->GetX2()-ar->GetX1(),ar->GetY2()-ar->GetY1());}

void HalfToThird(TCanvas * c)
{
  double xl,yl,xh,yh;
  c->GetPad(1)->GetPadPar(xl,yl,xh,yh);
  c->GetPad(1)->SetPad(xl,yl,.3,yh);
  c->GetPad(2)->GetPadPar(xl,yl,xh,yh);  
  c->GetPad(2)->SetPad(0.31,yl,xh,yh);
}

TH2F * hall_xy=NULL;
TH2F * hall_rz=NULL;
TGraph * all_xy=NULL;


TCut part(int start,int N){
  TString cs="";
  for (int i=start;i!=start+N;i++){
    cs+=Form("RoadCandidateArray.RoadSimHitArray.simhit.layer(%d)",i);
    if (i!=start+N-1) cs+=" ||";}
  return TCut(cs);}


void createHall(char * filename)
{
  std::cout<<"creating the background image into: "<<filename<<std::endl;
  TFile * f =TFile::Open(filename,"recreate");
  f->cd();
  hall_xy = new TH2F("hall_xy","All Hits (x,y) plane",Nr,-rmax_display,+rmax_display,Nr,-rmax_display,+rmax_display);
  hall_xy->SetMarkerColor(12);

  TCut TIB=part(0,4);
  TCut TOB=part(4,6);
  TCut TID=part(10,3);
  TCut TEC=part(13,9);
  TCut barrel=TIB || TOB;
  TCut endcap=TID || TEC;


  data->Draw("RoadCandidateArray.RoadSimHitArray.simhit.y:RoadCandidateArray.RoadSimHitArray.simhit.x>>hall_xy",barrel,"");
  hall_rz = new TH2F("hall_rz","All Hits (r,z) plane",Nz,-zmax_display,+zmax_display,Nr/2,0,+rmax_display);
  hall_rz->SetMarkerColor(12);
  data->Draw("RoadCandidateArray.RoadSimHitArray.simhit.r():RoadCandidateArray.RoadSimHitArray.simhit.z>>hall_rz","","");
  
  data->Draw("RoadCandidateArray.RoadSimHitArray.simhit.y:RoadCandidateArray.RoadSimHitArray.simhit.x",barrel,"");
  all_xy = grab("all_xy","All Hits (x,y) plane");
  all_xy->SetMarkerColor(12);

  f->cd();
  hall_xy->Write();
  hall_rz->Write();
  all_xy->Write();
  f->Close();
}

void retreiveHall(char * filename){
  std::cout<<"retreiving the background from: "<<filename<<std::endl;
  hall_xy= (TH2F *)TFile::Open(filename)->Get("hall_xy");
  hall_rz= (TH2F *)TFile::Open(filename)->Get("hall_rz");
  all_xy= (TGraph *)TFile::Open(filename)->Get("all_xy");
  all_xy->SetMarkerColor(12);
}

bool Has(TString what){
  data->GetEntry(_ie);
  return data->GetLeaf("RoadCandidateArray."+what)->GetValue(_is);}

bool HasBest(){ return Has("hasBest");}
//  data->GetEntry(_ie);
//  return data->GetLeaf("RoadCandidateArray.hasBest")->GetValue(_is);}

bool canDo()
{ data->GetEntry(_ie);
  int nc=data->GetLeaf("nRoadCandidate")->GetValue();
  return (_is<nc) ;}

  const double rmax_display=120;
  const double zmax_display=300;
  const int Nr=600;
  const int Nz=1200;

void Display(int ievent=0,int instance=0)
{
  static int count=0;
  gStyle->SetOptStat(0);
  _ie=ievent;_is=instance;
  bool hasBest=HasBest();
  if (!canDo()){ cout<<"can't do that, sorry"<<endl;return;}



  TCanvas * c =new TCanvas(Form("display_xy_E%d_I%d",ievent,instance),Form("basic event display (x,y) %d:%d",ievent,instance),800,800);  

  if (!hall_xy) { retreiveHall("/uscms/home/vlimant/work/CMSSW_1_3_0_pre5/src/RecoMuon/L3MuonAnalyzer/test/backgroundLayout.root");}
  if (!hall_xy) {cout<<"you must retreive a backgroud graphs somewhere."<<endl; return;}
  
  TGraph * sim_xy=NULL;
  TGraph * sim_rz=NULL;
  if(Has("hasSim")){
    sim_xy = XY("RoadSimHitArray.simhit","sim",8,25,1.5);
    sim_rz = RZ("RoadSimHitArray.simhit","sim",8,25,1.5);}
  
  TGraph * best_xy=NULL;
  TGraph * best_rz=NULL;
  TGraph * other_xy=NULL;
  TGraph * other_rz=NULL;
  TList * other_xy_l=NULL;
  TList * other_rz_l=NULL;

  //  if (hasBest){
  if (Has("hasBest")){
    best_xy = XY("best.RoadHitArray.rechit","best",4,28,0.8);
    best_rz = RZ("best.RoadHitArray.rechit","best",4,28,0.8);
    other_xy = XYmore("other",6,27,1.0);
    other_rz = RZmore("other",6,27,1.0);
    other_xy_l = XYlist("other",6,27,1.0);
    other_rz_l = RZlist("other",6,27,1.0);
  }


  TGraph * road_xy=NULL;
  road_xy = XY("RoadHitArray.rechit","road",2,4,1);
  TGraph * road_rz=NULL;
  road_rz = RZ("RoadHitArray.rechit","road",2,4,1);

  TGraph * road_w_xy=NULL;
  road_w_xy = XY("RoadHitArray.rechit","Wroad",2,5,1,false);
  TGraph * road_w_rz=NULL;
  road_w_rz = RZ("RoadHitArray.rechit","Wroad",2,5,1,false);

  TGraph * road_way_xy = XY("RoadSimHitArray.road","road_way",2,26,1.5);
  TGraph * road_way_rz = RZ("RoadSimHitArray.road","road_way",2,26,1.5);
  /*TGraph * road_way_xy=NULL;
    road_way_xy = XY("RoadHitArray.road","road_way",2,26,1.5);
    TGraph * road_way_rz=NULL;
    road_way_rz = RZ("RoadHitArray.road","road_way",2,26,1.5);*/



  TArrow *sim_xy_p=NULL;
  TArrow *sim_rz_p=NULL;
  if (Has("hasSim")){
    sim_xy_p= arrowXY("sim",8,1);
    sim_rz_p= arrowRZ("sim",8,1);}
  
  TArrow *pca_xy_p=NULL;
  TArrow *pca_rz_p=NULL;
  if (Has("hasPca")){
    pca_xy_p= arrowXY("pca",6,1);
    pca_rz_p= arrowRZ("pca",6,1);}

  TArrow *cpca_xy_p=NULL;
  TArrow *cpca_rz_p=NULL;
  if (Has("hasCpca")){
    cpca_xy_p= arrowXY("Cpca",2,1);
    cpca_rz_p= arrowRZ("Cpca",2,1);}

  TArrow *best_xy_p=NULL;
  TArrow *best_rz_p=NULL;
  if(hasBest){
    TArrow *best_xy_p= arrowXY("best.updated",4,1);
    TArrow *best_rz_p= arrowRZ("best.updated",4,1);}
  
  TLegend  * leg = new TLegend(0.7,0.7,0.99,0.99);
  leg->AddEntry(sim_xy,"Generated SimHit","p");
  if (sim_xy_p){leg->AddEntry(sim_xy_p,"Generated","l");}
  if(sim_xy_p){  leg->AddEntry(sim_xy_p,Form("P: %4.1f P_{T}: %4.1f",vector(sim_rz_p).Mod(),vector(sim_xy_p).Mod()),"");}
  if(pca_xy_p){  leg->AddEntry(pca_xy_p,"PCA","l");}
  if(pca_xy_p){  leg->AddEntry(pca_xy_p,Form("P: %4.1f P_{T}: %4.1f",vector(pca_rz_p).Mod(),vector(pca_xy_p).Mod()),"");}
  if(road_way_xy->GetN()!=0){  leg->AddEntry(road_way_xy,"Road way","p");}
  if(road_xy->GetN()!=0){  leg->AddEntry(road_xy,"Road RH","p");}
  if(road_w_xy->GetN()!=0){  leg->AddEntry(road_w_xy,"Road wrong RH","p");}
  if(cpca_xy_p){  leg->AddEntry(cpca_xy_p,"cPCA","l");}
  if(cpca_xy_p){  leg->AddEntry(cpca_xy_p,Form("P: %4.1f P_{T}: %4.1f",vector(cpca_rz_p).Mod(),vector(cpca_xy_p).Mod()),"");}
  if (best_xy) leg->AddEntry(best_xy,"RecoTrack RH","p");
  if (best_xy_p) {leg->AddEntry(best_xy_p,"RecoTrack","l");
    leg->AddEntry(best_xy_p,Form("P: %4.1f P_{T}: %4.1f",vector(best_rz_p).Mod(),vector(best_xy_p).Mod()),"");}
  if (other_xy_l) {TListIter it(other_xy_l);int ig=0;TGraph* gr;while(gr=(TGraph*)it()){leg->AddEntry(gr,Form("Other (%d)",ig++),"p");}
  }
  
  c->Clear();

  c->cd()->SetGrid();
  hall_xy->Draw();
  if(sim_xy){  sim_xy->Draw("same p");}
  if(best_xy) best_xy->Draw("same p");
  if (other_xy_l) { TListIter it(other_xy_l);TGraph* gr;while(gr=(TGraph*)it()){gr->Draw("same p");}}
  if(road_way_xy->GetN()!=0){  road_way_xy->Draw("same p");}
  if(road_xy->GetN()!=0){  road_xy->Draw("same p");}
  if(road_w_xy->GetN()!=0){  road_w_xy->Draw("same p");}
  if (best_xy_p) best_xy_p->Draw("same");
  if(sim_xy_p){  sim_xy_p->Draw("same");}
  if(pca_xy_p){  pca_xy_p->Draw("same"); }
  if(cpca_xy_p){  cpca_xy_p->Draw("same"); }

  leg->Draw("same");

  c =new TCanvas(Form("display_rz_E%d_I%d",ievent,instance),Form("basic event display (r,z) %d:%d",ievent,instance),1000,500);  
  c->cd()->SetGrid(); 

  hall_rz->Draw();
  if(sim_rz){  sim_rz->Draw("same p");}
  if(best_rz) best_rz->Draw("same p");
  if (other_rz_l) { TListIter it(other_rz_l);TGraph* gr;while(gr=(TGraph*)it()){gr->Draw("same p");}}
  if(road_way_rz->GetN()!=0){  road_way_rz->Draw("same p");}
  if(road_rz->GetN()!=0){  road_rz->Draw("same p");}
  if(road_w_rz->GetN()!=0){  road_w_rz->Draw("same p");}
  if(best_rz_p) best_rz_p->Draw("same");
  if(sim_rz_p){  sim_rz_p->Draw("same");}
  if(pca_rz_p){  pca_rz_p->Draw("same");}
  if(cpca_rz_p){  cpca_rz_p->Draw("same"); }

  TLegend  * rzleg = new TLegend(0.85,0.5,0.99,0.99);
  if(sim_rz){  rzleg->AddEntry(sim_rz,"Generated SimHit","p");}
  if(sim_xy_p){  rzleg->AddEntry(sim_xy_p,"Generated","l");}
  if(sim_xy_p){  rzleg->AddEntry(sim_xy_p,Form("P: %4.1f P_{T}: %4.1f",vector(sim_rz_p).Mod(),vector(sim_xy_p).Mod()),"");}
  if(pca_xy_p){  rzleg->AddEntry(pca_xy_p,"PCA","l");}
  if(pca_xy_p){  rzleg->AddEntry(pca_xy_p,Form("P: %4.1f P_{T}: %4.1f",vector(pca_rz_p).Mod(),vector(pca_xy_p).Mod()),"");}
  if(road_way_rz->GetN()!=0){  rzleg->AddEntry(road_way_rz,"Road way","p");}
  if(road_rz->GetN()!=0){  rzleg->AddEntry(road_rz,"Road RH","p");}
  if(road_w_rz->GetN()!=0){  rzleg->AddEntry(road_w_rz,"Road wrong RH","p");}
  if(cpca_xy_p){  rzleg->AddEntry(cpca_xy_p,"cPCA","l");}
  if(cpca_xy_p){  rzleg->AddEntry(cpca_xy_p,Form("P: %4.1f P_{T}: %4.1f",vector(cpca_rz_p).Mod(),vector(cpca_xy_p).Mod()),"");}
  if (best_rz) rzleg->AddEntry(best_rz,"RecoTrack RH","p");
  if (best_xy_p) {rzleg->AddEntry(best_xy_p,"RecoTrack","l");
    rzleg->AddEntry(best_xy_p,Form("P: %4.1f P_{T}: %4.1f",vector(best_rz_p).Mod(),vector(best_xy_p).Mod()),"");}
  if(other_rz_l){TListIter it(other_rz_l); int ig=0;TGraph* gr;while(gr=(TGraph*)it()){rzleg->AddEntry(gr,Form("Other (%d)",ig++),"p");}}

  rzleg->Draw();

  c =new TCanvas(Form("display_xyR_E%d_I%d",ievent,instance),Form("rotated event display %d:%d",ievent,instance),800,800);  
  double angle = TMath::Pi()/2. -TVector2(sim_xy_p->GetX2()-sim_xy_p->GetX1(),sim_xy_p->GetY2()-sim_xy_p->GetY1()).Phi();


  TGraph * all_xy_r= XYrotation(all_xy,angle);
  TGraph * sim_xy_r=NULL;
  if(sim_xy){  sim_xy_r = XYrotation(sim_xy,angle);}
  TGraph * best_xy_r=NULL;
  if (best_xy){ best_xy_r = XYrotation(best_xy,angle);}
  TList * other_xy_l_r;
  if (other_xy_l){TListIter it(other_xy_l); other_xy_l_r = new TList(); TGraph* gr;while(gr=(TGraph*)it()){other_xy_l_r->Add(XYrotation(gr,angle));}}

  TGraph * road_way_xy_r = XYrotation(road_way_xy,angle);
  TGraph * road_xy_r = XYrotation(road_xy,angle);
  TGraph * road_w_xy_r = XYrotation(road_w_xy,angle);

  TArrow *sim_xy_p_r=NULL;
  if (sim_xy_p){sim_xy_p_r = XYrotation(sim_xy_p,angle);}
  TArrow *best_xy_p_r=NULL;
  if (best_xy_p){best_xy_p_r = XYrotation(best_xy_p,angle);}
  TArrow *pca_xy_p_r=NULL;
  if (pca_xy_p){pca_xy_p_r= XYrotation(pca_xy_p,angle);}
  TArrow *cpca_xy_p_r=NULL;
  if (cpca_xy_p){cpca_xy_p_r= XYrotation(cpca_xy_p,angle);}

  TLegend * cleg = new TLegend(0.7,0.7,0.99,0.99);
  if(sim_xy_r){  cleg->AddEntry(sim_xy_r,Form("SimHit: %d",sim_xy_r->GetN()),"p");}
  if(road_way_xy_r->GetN()!=0){  cleg->AddEntry(road_way_xy_r,"Road way","p");}
  if(road_xy_r->GetN()!=0){  cleg->AddEntry(road_xy_r,Form("Road RH: %d",road_xy_r->GetN()),"p");}
  if(road_w_xy_r->GetN()!=0){  cleg->AddEntry(road_w_xy_r,Form("Road wrong RH: %d",road_w_xy_r->GetN()),"p");}
  if(best_xy) cleg->AddEntry(best_xy,Form("RecoTrack RH: %d",best_xy->GetN()),"p");
  if (other_xy_l){TListIter it(other_xy_l);int ig=0;TGraph* gr; while(gr=(TGraph*)it()){cleg->AddEntry(gr,Form("Other (%d) RH: %d",ig++,gr->GetN()),"p");}}
  

  c->Clear();
  c->cd()->SetGrid();
  all_xy_r->Draw("ap");
  if(sim_xy_r){  sim_xy_r->Draw("same p");}
  if (best_xy_r) best_xy_r->Draw("same p");
  if (other_xy_l_r){TListIter it(other_xy_l_r); TGraph* gr;while(gr=(TGraph*)it()){gr->Draw("same p");}}
  if(road_way_xy_r->GetN()!=0){  road_way_xy_r->Draw("same p");}
  if(road_xy_r->GetN()!=0){  road_xy_r->Draw("same p");}
  if(road_w_xy_r->GetN()!=0){  road_w_xy_r->Draw("same p");}
  if(sim_xy_p_r){  sim_xy_p_r->Draw("same");}
  if (best_xy_p_r) best_xy_p_r->Draw("same");
  if(pca_xy_p_r){  pca_xy_p_r->Draw("same"); }
  if(cpca_xy_p_r){  cpca_xy_p_r->Draw("same"); }
  cleg->Draw();

  count++;
}

void Display(TString displayLocation){
  int ne = data->GetEntries();
  for (int ie=0;ie!=ne;++ie){
    data->SetBranchAddress("aDump",&_ed);
    data->GetEntry(ie);
    int ns=_ed->nRoadCandidate;
    for (int is=0;is!=ns;++is){
      cout <<"doing: "<<ie<<" "<<is<<endl;
      Display(ie,is);
    }
  }
  
  plotAll(displayLocation,".eps");
  plotAll(displayLocation,".gif");
  plotAll(displayLocation,".root");
}
