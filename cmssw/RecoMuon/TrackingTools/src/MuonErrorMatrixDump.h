#ifndef DUMPCLASS_H
#define DUMPCLASS_H

#include "TMath.h"

#ifndef CINT
#include "DataFormats/CLHEP/interface/AlgebraicObjects.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#endif

class MuonErrorMatrixDump {
 public:
  MuonErrorMatrixDump();
  virtual ~MuonErrorMatrixDump();
  void clean();

  unsigned int run;
  unsigned int event;

  bool hasSim;
  double sim_vector[5];
  double sim_x[3];
  double sim_p[3];

  bool hasReco;
  double reco_vector[5];
  double reco_error[15];
  double reco_x[3];
  double reco_p[3];
  int reco_Nhit;

  // accessors to momentum and position
  double P(bool s=false){ return (s? mom_(sim_p): mom_(reco_p));}
  double px(bool s=false){ return p_(0,s);}
  double py(bool s=false){ return p_(1,s);}
  double pz(bool s=false){ return p_(2,s);}
  double pT(bool s=false){ return (s? pT_(sim_p) : pT_(reco_p));}

  double eta(bool s=false){ return (s? eta_(sim_p) : eta_(reco_p));}
  double feta(bool s=false){ return (s? feta_(sim_p) : feta_(reco_p));}
  double phi(bool s=false){ return (s? phi_(sim_p) : phi_(reco_p));}

  double x(bool s=false){ return x_(0,s);}
  double y(bool s=false){ return x_(1,s);}
  double z(bool s=false){ return x_(2,s);}
  double d(bool s=false){ return (s? pT_(sim_x): pT_(reco_x));}

  // accessors to curvilinear parameters
  double QoverP(bool s=false){ return v_(0,s); }
  double lambda(bool s=false){ return v_(1,s); }
  double phi0(bool s=false){ return v_(2,s); }
  double xT(bool s=false){ return v_(3,s); }
  double yT(bool s=false){ return v_(4,s); }

  // accessors to curvilinear errors
  double QoverP_e(){ return e_(0,0);}
  double lambda_e(){ return e_(1,1);}
  double phi0_e(){ return e_(2,2);}
  double xT_e(){ return e_(3,3);}
  double yT_e(){ return e_(4,4);}
  double corr(int i, int j) { return e_(i,j);}
  double cov(int i, int j) { return mat_(i,j);}


  //accessor to curvilinear parameters residual (reco-gen)
  double QoverP_r(){ return r_(0); }
  double lambda_r(){ return r_(1); }
  double phi0_r(){ return r_(2); }
  double xT_r(){ return r_(3); }
  double yT_r(){ return r_(4); }

  // accessor to the pulls
  double QoverP_p(){ return pu_(0); }
  double lambda_p(){ return pu_(1); }
  double phi0_p(){ return pu_(2); }
  double xT_p(){ return pu_(3); }
  double yT_p(){ return pu_(4); }

#ifndef CINT
  // fill method. only visible when not compiled by root interactively
  void fillSim(const AlgebraicVector5 &v, const GlobalPoint &x, const GlobalVector &p);
  void fillReco(const AlgebraicVector5 &v, const AlgebraicSymMatrix55 & e, const GlobalPoint &x, const GlobalVector &p);
#endif

 private:
  int index_(int i,int j){  if (j>i) return index_(j,i); return i*5+(i-j);}
  double mom2_(double * p) { return pT2_(p)+p[2]*p[2];}
  double mom_(double * p) { return sqrt(mom2_(p));}
  double p_(int i,bool s=false){ return (s ?sim_p[i]: reco_p[i]);}
  double pT2_(double * p){ return p[0]*p[0]+p[1]*p[1];}
  double pT_(double * p){ return sqrt(pT2_(p));}
  double eta_(double * p) { double P=mom_(p); double z = p[2]; return 0.5*log((P+z)/(P-z));}
  double feta_(double * p) { return fabs(eta_(p));}
  double phi_(double *p) { return TMath::Sign(TMath::ACos(p[1]/pT_(p)),p[2]);}
  double x_(int i, bool s=false){ return (s?sim_x[i]: reco_x[i]);}
  double v_(int i, bool s=false){ return (s?sim_vector[i]:reco_vector[i]);}
  double r_(int i){ return v_(i,false) - v_(i,true);}
  double mat_(int i, int j){ return reco_error[index_(i,j)];}
  double e_(int i, int j){ return ((i==j)? sqrt(mat_(i,j)) : mat_(i,j)/sqrt(mat_(i,i)*mat_(j,j)));}
  double pu_(int i){ return (v_(i,false) - v_(i,true))/e_(i,i);}

#ifdef CINT    
  ClassDef(MuonErrorMatrixDump, 0)
#endif

};


#ifdef CINT
ClassImp(MuonErrorMatrixDump)
#endif




MuonErrorMatrixDump::MuonErrorMatrixDump(){}
MuonErrorMatrixDump::~MuonErrorMatrixDump(){clean();}
void MuonErrorMatrixDump::clean(){
  hasSim= false;
  hasReco= false;
  for (uint i=0; i!=5; i++){
    sim_vector[i]=0;
    reco_vector[i]=0;
    if (i<=2){
      sim_x[i]=0;
      sim_p[i]=0;
      reco_x[i]=0;
      reco_p[i]=0;
    }
    for (uint j=0; j!=5 ; j++){
      reco_error[index_(i,j)]=0;
    }
  }
  reco_Nhit=0;
}

#ifndef CINT
void MuonErrorMatrixDump::fillSim(const AlgebraicVector5 &v, const GlobalPoint &x, const GlobalVector &p){
  hasSim= true;
  for (uint i=0; i!=5; i++){sim_vector[i]=v(i);}
  sim_x[0]=x.x();
  sim_x[1]=x.y();
  sim_x[2]=x.z();

  sim_p[0]=p.x();
  sim_p[1]=p.y();
  sim_p[2]=p.z();
}

void MuonErrorMatrixDump::fillReco(const AlgebraicVector5 &v, const AlgebraicSymMatrix55 & e, const GlobalPoint &x, const GlobalVector &p){
  hasReco = true;
  for (uint i=0; i!=5; i++){
    reco_vector[i]=v(i);
    for (uint j=0; j!=5; j++){
      reco_error[index_(i,j)]=e(i,j);
    }}
  reco_x[0]=x.x();
  reco_x[1]=x.y();
  reco_x[2]=x.z();

  reco_p[0]=p.x();
  reco_p[1]=p.y();
  reco_p[2]=p.z();
}
#endif

#endif
