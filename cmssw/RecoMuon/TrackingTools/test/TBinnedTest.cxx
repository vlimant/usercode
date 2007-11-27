#include "RecoMuon/TrackingTools/interface/TBinned.h"
#include "TRandom2.h"
#include <fstream>

class cumulator{
public:
  cumulator():N(0),SUM(0),SUMSQUARE(0){}
  double n(){ return N;}
  double sum(){ return SUM;}
  double sum2(){ return SUMSQUARE;}

  double mean(){
    if (N==0) return 0;
    else return SUM/static_cast<double>(N);
  }
  double rms(){
    if (N<=1) return 0;
    else {
      double m=mean();
      double m2=SUMSQUARE/static_cast<double>(N);
      return sqrt(m2-m*m);
    }
  }
  void cumulate(double x){
    N++;
    SUM+=x;
    SUMSQUARE+=x*x;
  }
  void operator+=(double c){cumulate(c);}
private:
  uint N;
  double SUM;
  double SUMSQUARE;
};


ostream & operator<<(ostream & o, cumulator & c){
  o<<"mean: "<<c.mean()<<", rms: "<<c.rms()<<", entries: "<<c.n();
  return o;
}

int main()
{
  
  //  TBinned<double> binned(3, "aName", "aTitle");
  TBinned<cumulator> binned(3, "aName", "aTitle");

  //  binned.setWithUnderOverFlow();

  binned.setAxis(0,"p_{T}", 10,0,100);// 10 bins from 0 to 100 GeV
  double etaBins[5]={0, 0.3, 0.33, 0.4, 1.5};
  binned.setAxis(1, "#eta", 4, etaBins);//variable bins in second axis
  binned.setAxis(2, "#phi", 5, 0,1);//5 bins in phi
  //  binned.setAxis(3, "nHits", 10, 0.5,10.5);//10 bins in unpredefined axis

  double pT=14.0, eta=0.9, phi=0.5;
  double pTmax=120, etamax=2.0, phimax=1.2;
  
  double valuemax = 100;
  TRandom2 r;
  
  if (binned.check()){
    for (uint i=0; i!=100000; i++){
      pT = r.Rndm()*pTmax;
      eta = r.Rndm()*etamax;
      phi = r.Rndm()*phimax;
      
      double val = r.Rndm()*valuemax;
      binned[pT][eta][phi]() += val;

      if (i%1000 ==0) { 
	std::cout<<i<<":   "<<pT<<", "<<eta<<", "<<phi<<" = "<<val<<std::endl;
	binned[pT][eta][phi].getContent(std::cout);
      }

    }
    
    binned.print(std::cout);
    std::ofstream dFile("testDump.txt");
    binned.dump(dFile);
  }
  return 0;
}
