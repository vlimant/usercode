#ifndef DUMPCLASS_H
#define DUMPCLASS_H


#define SEED_VERSION


#include <TClonesArray.h>
#include <TClass.h>
#include <TObject.h>
#include <TMath.h>

#include <vector>
#include <iostream>
using namespace std;

enum Hit__Layer {TIB_1,TIB_2,TIB_3,TIB_4,
	    TOB_1,TOB_2,TOB_3,TOB_4,TOB_5,TOB_6,
	    TID_1,TID_2,TID_3,
	    TEC_1,TEC_2,TEC_3,TEC_4,TEC_5,TEC_6,TEC_7,TEC_8,TEC_9};

class Hit: public TObject {
 public:
  Hit();
  virtual ~Hit();
  
  bool _3D;
  double x,y,z;
  double v11,v12,v13,v22,v23,v33;
  double x_l,y_l;
  double v11_l,v12_l,v13_l;
  double chi2_road;
  double chi2_track;

  double r(){return sqrt(x*x+y*y);}
  double rho() {return sqrt(x*x+y*y+z*z);}
  double eta(bool s =true) {double L=rho();if (s) return 0.5*log((L+z)/(L-z)); else return 0.5*fabs(log((L+z)/(L-z))); }
  double phi() {return TMath::Sign(TMath::ACos(x/r()),y);}
  double r_l(){return sqrt(x_l*x_l+y_l*y_l);}
  bool layer(Hit__Layer l)
  {
    static double radiuses[10]={24,32,42,50,62,68,78,85,95,106};
    static double Zes[12]={77.5, 91 , 105,132 , 146 , 160 , 174 , 188 , 205, 225 , 245 , 267};
    static double Zes_w[12]={8, 8 , 8, 6 , 6 , 6 , 6 , 6 , 7, 7 , 7 , 7};

    double radius=radiuses[min((int)l,9)];
    double Z=Zes[max(0,(int)l-10)];
    double Zw=Zes_w[max(0,(int)l-10)];

    if (z<0) Z*=-1;

    if (l<10)
      {
	if(l<4)
	  { return ((fabs(r()-radius)<4) && (fabs(z)<Zes[0]-8)) ;}
	else
	  { return ((fabs(r()-radius)<4) && (fabs(z)<Zes[3]-8)) ;}
      }
    else if(l<22)
      {
	if (l<13)
	  { return ((fabs(z-Z)<Zw) && (r()<radiuses[4]-5));}
	else
	  { return (fabs(z-Z)<Zw) ;}
      }
    return false;
  }

  void clean();
  
  ClassDef(Hit, 1)
};


class State : public TObject {
 public:
  State();
  virtual ~State();
  void clean();

  double px,py,pz;
  double x,y,z;
  //curvilinear
  double v11,v12,v13,v14,v15,v22,v23,v24,v25,v33,v34,v35,v44,v45,v55;
  //cartesian
  double c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66;
  bool hasLocal;
  double x_l,y_l;
  double v11_l,v12_l,v22_l;

  double pt(){return sqrt(px*px+py*py);}
  double p(){return  sqrt(px*px+py*py+pz*pz);}
  double eta_pos(bool s =true)  {double P=rho(); if(s) return 0.5*log((P+z)/(P-z)); else return 0.5*fabs(log((P+z)/(P-z)));}
  double eta(bool s =true)  {double P=p();if (s) return 0.5*log((P+pz)/(P-pz)); else return 0.5*fabs(log((P+pz)/(P-pz))); }
  double phi() {return TMath::Sign(TMath::ACos(px/pt()),py);}
  double phi_pos() {return TMath::Sign(TMath::ACos(x/r()),y);}
  double r(){return sqrt(x*x+y*y);}
  double rho(){return sqrt(x*x+y*y+z*z);}

  ClassDef(State, 1)
};

class DressedHit : public TObject {
 public:
  DressedHit();
  virtual ~DressedHit();
  void clean();

  Hit rechit; //the rechit part
  Hit simhit; //the simhit part
  int matchcode; //the match code between the rechit and the simhit (0:Noise, 1:OtherTrack, 2:Primary, 3:Secondary)

  bool hasState; //a state is available
  State road; //state at the surface of the rechit

  bool noise(){return matchcode==0;}
  bool other(){return matchcode==1;}
  bool primary(){return matchcode==2;}
  bool secondary(){return matchcode==3;}
  bool own(){return (matchcode==2 || matchcode==3);}

  ClassDef(DressedHit, 1)
};


class Candidate : public TObject{
 public:
  Candidate();
  virtual ~Candidate();
  void clean();

  vector<DressedHit> RoadHitArray; //vector of hit on the track candidate
  Int_t     nRoadHit; //number of rechit on the track candidate
  State updated; //state of the track candidate at vertex

  DressedHit * add();
  inline DressedHit *  current(){return (DressedHit*) &(RoadHitArray)[nRoadHit-1];}
  
  inline const vector<DressedHit> & get() {return RoadHitArray;}

  int noise(){ int sum=0; for (vector<DressedHit>::iterator it= RoadHitArray.begin(); it!=RoadHitArray.end();it++){ if (it->noise()) sum++;} return sum;}
  int other(){ int sum=0; for (vector<DressedHit>::iterator it= RoadHitArray.begin(); it!=RoadHitArray.end();it++){ if (it->other()) sum++;} return sum;}
  int primary(){ int sum=0; for (vector<DressedHit>::iterator it= RoadHitArray.begin(); it!=RoadHitArray.end();it++){ if (it->primary()) sum++;} return sum;}
  int own(){ int sum=0; for (vector<DressedHit>::iterator it= RoadHitArray.begin(); it!=RoadHitArray.end();it++){ if (it->own()) sum++; } return sum;}
  int secondary(){ int sum=0; for (vector<DressedHit>::iterator it= RoadHitArray.begin(); it!=RoadHitArray.end();it++){ if (it->secondary()) sum++;} return sum;}
  
  ClassDef(Candidate,1)
};




class RoadCandidate : public TObject {
 public:
  RoadCandidate();
  virtual ~RoadCandidate();
  void clean();

  //contains the rechits collected in the road
  vector<DressedHit> RoadHitArray; //vector of hit within the road
  Int_t     nRoadHit;
  //contains the true simhits associated with this road candidate //dressed hits
  vector<DressedHit> RoadSimHitArray; //vector of hit associated with the simulated seed
  Int_t     nRoadSimHit;
  //contains other candidate than best
  vector<Candidate> CandidateArray; //vector of trackcandidate, containing the best one
  Int_t     nCandidate;

#ifdef SEED_VERSION
  Candidate * addSeed();
  //contains the trajectory seed for this roadcandidate
  vector<Candidate> SeedArray;
  Int_t nSeed;
  inline Candidate * currentSeed(){return (Candidate*) &(SeedArray)[nSeed-1];}
#else
  Candidate * addSeed(){return NULL;}//so that it will crash
#endif

  State  seed; //initial state in the muon detector
  int Nrechits; //number of rechits on this muon seed
  bool hasSim; //the seed has a simulated match
  State  sim; //simulated state matched to the seed. at vertex
  bool hasPca; //propagation to vertex succeeded
  State  pca; //state of the seed at vertex. no constraint
  bool hasCpca; //combination with IP succeeded
  State  Cpca; //state of the sedd at vertec after constraint
  bool hasBest; //a track candidate has been found for this seed
  Candidate best; //state+hitlist of track candidate found. at vertex
  bool hasOutOut; //OutOut state available
  State outout; //state of muon seed at outer barrel. IP constraint
  bool hasOutIn; //OutIn state available
  State outin; //state of muon seed at outer barrel. no IP constraint

  unsigned int nOHPM; //number of hit per module cut

  bool hasSeed; //the roadcandidate might not be attached to a seed

  inline DressedHit * currentRecHit() {return (DressedHit *) &(RoadHitArray)[nRoadHit-1];}
  inline DressedHit * currentSimHit() { return (DressedHit *) &(RoadSimHitArray)[nRoadSimHit-1];}
  inline Candidate * currentCandidate(){return(Candidate *) &(CandidateArray)[nCandidate-1];}

  DressedHit * addRecHit();
  DressedHit * addSimHit();
  Candidate * addCandidate();


  int noise(){ int sum=0; for (vector<DressedHit>::iterator it= RoadHitArray.begin(); it!=RoadHitArray.end();it++){ if (it->matchcode==0) sum++;} return sum;}
  int other(){ int sum=0; for (vector<DressedHit>::iterator it= RoadHitArray.begin(); it!=RoadHitArray.end();it++){ if (it->matchcode==1) sum++;} return sum;}
  int primary(){ int sum=0; for (vector<DressedHit>::iterator it= RoadHitArray.begin(); it!=RoadHitArray.end();it++){ if (it->matchcode==2) sum++;} return sum;}
  int secondary(){ int sum=0; for (vector<DressedHit>::iterator it= RoadHitArray.begin(); it!=RoadHitArray.end();it++){ if (it->matchcode==3) sum++;} return sum;}
  int own(){ int sum=0; for (vector<DressedHit>::iterator it= RoadHitArray.begin(); it!=RoadHitArray.end();it++){ if (it->matchcode==2) sum++; else if(it->matchcode==3) sum++; } return sum;}
  
#ifndef SEED_VERSION
  ClassDef(RoadCandidate,1)//added the seedarray
#else
    ClassDef(RoadCandidate,2)
#endif
};


class EventDump : public TObject {
 public:
  EventDump();
  virtual ~EventDump();
  void clean();

  
  TClonesArray  *RoadCandidateArray; //array of seed road candidate
  static TClonesArray *gRoadCandidateArray; 
  static int iRoadCandidate;
  Int_t nRoadCandidate;
  int run; //run number
  int evt; //event number
  int Nseeds; //number of seed for this event
  int Nsim; //number of simulated muon in the event
  int Nmatched; 
  int Nsta; //number of Stand-Alone muon in the event
  int Nl1; //number of L1 muon in the event
  int NtrackCandidate; //number of track candidate inserted in the event

  RoadCandidate * add();
  inline RoadCandidate * current()   {return (RoadCandidate *) (*RoadCandidateArray)[iRoadCandidate];}
  inline const TClonesArray * get(){return RoadCandidateArray;}
  inline void reset(){iRoadCandidate=-1;}
  inline RoadCandidate * next() {iRoadCandidate++; return current();}

  //only usefull within root interactive, because recurisivity cannot go up above 2 in branches... how dumb !
  std::vector<std::pair<float,float> > xy(unsigned int i);
  std::vector<std::pair<float,float> > rz(unsigned int i);

  ClassDef(EventDump, 1)
};

//template function to add an item in a array pointer or a vector
template <class A> A* ADD(std::vector<A> & vector, int & N)
{
  vector.push_back(A());
  N=vector.size();
  return &(vector[N-1]);
}


template <class A> A* ADD(TClonesArray  * array,int & N)
{
  if (N>array->GetEntriesFast()) { std::cout<<"cannot add more than: "<<array->GetEntriesFast()<<" objects to the TClonesArray, with class: "<<array->GetClass()->GetName()<<std::endl;}
  TClonesArray &objects = *array;
  A* object = new(objects[N++]) A();
  return object;
}

#ifndef CINT

class TrajectoryStateOnSurface;
class TransientTrackingRecHit;
class FreeTrajectoryState;

void operator <<(Hit & rhs,const TrajectoryStateOnSurface & tsos); //generic for PSimHit when it has been propagated to a surface
void operator <<(Hit & rhs ,const TransientTrackingRecHit & rechit);

void operator <<(State & rhs, const FreeTrajectoryState &fts);
void operator <<(State & rhs, const TrajectoryStateOnSurface & tsos);

void operator<<(Candidate & rhs, TrajectoryStateOnSurface & tsos);

#endif



#endif
