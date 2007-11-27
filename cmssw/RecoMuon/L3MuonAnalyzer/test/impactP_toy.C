{
  TRandom2 r;

  const double transverse=15; //microns
  TH2D * h2_v = new TH2D("h2_v","",100,-2*transverse,2*transverse,100,-2*transverse,2*transverse);
  TH1D * h1_phi = new TH1D("h1_phi","",100,0,2*TMath::Pi());

  TH2D * h2_b = new TH2D("h2_b","",100,-2*transverse,2*transverse,100,-2*transverse,2*transverse);
  TH1D * h1_b = new TH1D("h1_b","",100,0,2*transverse);

  const int max=10;
  for (int i=0;i!=max;++i){
    double r=r.Gaus(0,transverse);
    double phi_v=2*TMath::Pi()*r.Rndm();
    x= r*cos(phi_v);
    y= r*sin(phi_v);
    h2_v->Fill(x,y);

    double phi =2*TMath::Pi()*r.Rndm();
    h1_phi->Fill(phi);
    
    TVector2 v(x,y);
    TVector2 u(cos(phi),sin(phi));

    TVector2 b= v - u.Dot(v)*u;
    h2_b->Fill(b.X(),b.Y());
    h1_b->Fill(b.Mod());
  }

  TCanvas * c = new TCanvas("c","***",800,800);
  c->Divide(2,2);
  c->cd(1)->SetGrid();
  h2_v->Draw();
  
  c->cd(2)->SetGrid();
  h1_phi->Draw();

  c->cd(3)->SetGrid();
  h2_b->Draw();

  c->cd(4)->SetGrid();
  h1_b->Draw();


}
