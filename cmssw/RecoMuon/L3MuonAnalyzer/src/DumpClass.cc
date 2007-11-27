#include <RecoMuon/L3MuonAnalyzer/interface/DumpClass.h>

ClassImp(Candidate);
ClassImp(Hit);
ClassImp(DressedHit);
ClassImp(State);
ClassImp(RoadCandidate);
ClassImp(EventDump);

TClonesArray  * EventDump::gRoadCandidateArray = 0;
int EventDump::iRoadCandidate=-1;

//-------State ---------
State::State(){clean();}
State::~State(){;}    

void State::clean(){
  hasLocal=false;
  px=0;py=0;pz=0;
  x=0;y=0;z=0;
  v11=0;v12=0;v13=0;v14=0;v15=0;v22=0;v23=0;v24=0;v25=0;v33=0;v34=0;v35=0;v44=0;v45=0;v55=0;
  c11=0;c12=0;c13=0;c14=0;c15=0;c16=0;c22=0;c23=0;c24=0;c25=0;c26=0;c33=0;c34=0;c35=0;c36=0;c44=0;c45=0;c46=0;c55=0;c56=0;c66=0;
  x_l=0;y_l=0;
  v11_l=0;v12_l=0;v22_l=0;
}
//--------------------

//------Candidate-------
Candidate::Candidate(){
  RoadHitArray.clear();
  nRoadHit=0;
  clean();
}
Candidate::~Candidate(){clean();}
DressedHit *Candidate::add(){return ADD<DressedHit>(RoadHitArray,nRoadHit);}
void Candidate::clean()
{
  RoadHitArray.clear();
  nRoadHit=0;
  updated.clean();
}
//--------------------

//------Hit------------
Hit::Hit(){clean();}
Hit::~Hit(){clean();}
void Hit::clean(){
  _3D=false;
  x=0;y=0;z=0;
  v11=0;v12=0;v13=0;v22=0;v23=0;v33=0;
  x_l=0;y_l=0;   
  v11_l=0;v12_l=0;v13_l=0;
  chi2_road=0; 
  chi2_track=0;
}
//--------------------

//---DressedHit--------
DressedHit::DressedHit(){clean();}
DressedHit::~DressedHit(){clean();}
void DressedHit::clean(){
  rechit.clean();
  simhit.clean();
  hasState=false;
  road.clean();
}
//--------------------

//---RoadCandidate-----
RoadCandidate::RoadCandidate()
{
  nRoadHit=0;
  RoadHitArray.clear();

  nRoadSimHit=0;
  RoadSimHitArray.clear();

  nCandidate=0;
  CandidateArray.clear();

  clean();
}
RoadCandidate::~RoadCandidate(){
  clean();
}
void RoadCandidate::clean()
{
  for (vector<Candidate>::iterator it = CandidateArray.begin();it!=CandidateArray.end();it++)
    {it->clean();}

  RoadHitArray.clear();
  RoadSimHitArray.clear();
  CandidateArray.clear();
#ifdef SEED_VERSION
  SeedArray.clear();
  nSeed=0;
#endif
  nRoadHit=0;
  nRoadSimHit=0;
  nCandidate=0;


  seed.clean();
  hasSeed=false;
  Nrechits=0;
  hasSim=false;
  sim.clean();
  hasPca=false;
  pca.clean();
  hasCpca=false;
  Cpca.clean();
  hasBest=false;
  best.clean();

  hasOutOut=false;
  outout.clean();
  hasOutIn=false;
  outin.clean();
}

DressedHit * RoadCandidate::addRecHit(){return ADD<DressedHit>(RoadHitArray,nRoadHit);}
DressedHit * RoadCandidate::addSimHit(){return ADD<DressedHit>(RoadSimHitArray,nRoadSimHit);}
Candidate * RoadCandidate::addCandidate() { return ADD<Candidate>(CandidateArray,nCandidate);}
#ifdef SEED_VERSION
Candidate * RoadCandidate::addSeed(){return ADD<Candidate>(SeedArray,nSeed);}
#endif

//--------------------

//-------EventDump----
EventDump::EventDump(){
  if(!gRoadCandidateArray) gRoadCandidateArray = new TClonesArray("RoadCandidate",10);
  RoadCandidateArray=gRoadCandidateArray;

  nRoadCandidate=0;
  clean();
  iRoadCandidate=-1;
}
EventDump::~EventDump(){clean();}
RoadCandidate *EventDump::add(){iRoadCandidate++;return ADD<RoadCandidate>(RoadCandidateArray,nRoadCandidate);}
void EventDump::clean()
{
  TIter it(RoadCandidateArray);
  RoadCandidate * r;
  while ((bool)(r=(RoadCandidate *) it.Next())) {r->clean();}
  RoadCandidateArray->Clear();
  nRoadCandidate=0;
  iRoadCandidate=-1;

  run=0;evt=0;
#ifdef SEED_VERSION
  Nseeds=0;
#endif

Nsim=0;Nmatched=0;Nsta=0;Nl1=0;NtrackCandidate=0;
}
std::vector<std::pair<float,float> > EventDump::xy(unsigned int i){
  std::vector<std::pair<float,float> > result;
  RoadCandidate * rc = (RoadCandidate*) RoadCandidateArray->At(i);
  std::vector<Candidate> cv = rc->CandidateArray;
  uint ci=0;
  for (; ci != cv.size();++ci){
    uint rhi=0;
    for (;rhi!=cv[ci].RoadHitArray.size();++rhi){
      result.push_back(std::make_pair(cv[ci].RoadHitArray[rhi].rechit.x, cv[ci].RoadHitArray[rhi].rechit.y));
    }
  }
  return result;
}

std::vector<std::pair<float,float> > EventDump::rz(unsigned int i){
  vector<std::pair<float,float> > result;
  RoadCandidate * rc = (RoadCandidate*) RoadCandidateArray->At(i);
  std::vector<Candidate> cv = rc->CandidateArray;
  uint ci=0;
  for (; ci != cv.size();++ci){
    uint rhi=0;
    for (;rhi!=cv[ci].RoadHitArray.size();++rhi){
      result.push_back(std::make_pair(cv[ci].RoadHitArray[rhi].rechit.r(), cv[ci].RoadHitArray[rhi].rechit.z));
    }
  }
  return result;
}
//--------------------


#ifndef CINT

#include <TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h>
#include <TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h>


void operator <<(Hit & rhs,const TrajectoryStateOnSurface & tsos) 
//generic for PSimHit when it has been propagated to a surface
//void operator <<(Hit & rhs,const PSimHit & simhit)
{
  GlobalPoint gX(0,0,0);
  LocalPoint X(0,0,0);
  if (tsos.isValid()) //this is not good for business !!!
    {
      gX=tsos.globalPosition();
      X=tsos.localPosition();
    }
  rhs.x=gX.x();rhs.y=gX.y();rhs.z=gX.z();
  rhs.v11=0;rhs.v12=0;rhs.v13=0;rhs.v22=0;rhs.v23=0;rhs.v33=0;//no error, this is supposed to be a psimhit
  rhs.x_l=X.x();
  rhs.y_l=X.y();
  rhs.v11_l=0;rhs.v12_l=0;rhs.v13_l=0;
}

void operator <<(Hit & rhs ,const TransientTrackingRecHit & rechit)
{
  GlobalPoint gX= rechit.globalPosition();
  rhs.x=gX.x();rhs.y=gX.y();rhs.z=gX.z();
  GlobalError gm=rechit.globalPositionError();
  rhs.v11=gm.cxx();rhs.v12=gm.cyx();rhs.v13=gm.czx();rhs.v22=gm.cyy();rhs.v23=gm.czy();rhs.v33=gm.czz();
  //  rhs.v11=0;rhs.v12=0;rhs.v13=0;rhs.v22=0;rhs.v23=0;rhs.v33=0;
  LocalPoint X=rechit.localPosition();
  rhs.x_l=X.x();
  rhs.y_l=X.y();
  LocalError m=rechit.localPositionError();
  rhs.v11_l=m.xx();rhs.v12_l=m.xy();rhs.v13_l=m.yy();
  //  rhs.v11_l=0;rhs.v12_l=0;rhs.v13_l=0;
}



void operator <<(State & rhs, const FreeTrajectoryState &fts)
{
  GlobalVector p = fts.momentum();
  rhs.px=p.x();rhs.py=p.y();rhs.pz=p.z();
  GlobalPoint X = fts.position();
  rhs.x=X.x();rhs.y=X.y();rhs.z=X.z();
  AlgebraicSymMatrix m = fts.curvilinearError().matrix();
  rhs.v11=m(1,1);
  rhs.v12=m(1,2);
  rhs.v13=m(1,3);
  rhs.v14=m(1,4);
  rhs.v15=m(1,5);


  rhs.v22=m(2,2);
  rhs.v23=m(2,3);
  rhs.v24=m(2,4);
  rhs.v25=m(2,5);

  rhs.v33=m(3,3);
  rhs.v34=m(3,4);
  rhs.v35=m(3,5);

  rhs.v44=m(4,4);
  rhs.v45=m(4,5);

  rhs.v55=m(5,5);

  AlgebraicSymMatrix mc = fts.cartesianError().matrix();

  rhs.c11=mc(1,1);
  rhs.c12=mc(1,2);
  rhs.c13=mc(1,3);
  rhs.c14=mc(1,4);
  rhs.c15=mc(1,5);
  rhs.c16=mc(1,6);

  rhs.c22=mc(2,2);
  rhs.c23=mc(2,3);
  rhs.c24=mc(2,4);
  rhs.c25=mc(2,5);
  rhs.c26=mc(2,6);

  rhs.c33=mc(3,3);
  rhs.c34=mc(3,4);
  rhs.c35=mc(3,5);
  rhs.c36=mc(3,6);

  rhs.c44=mc(4,4);
  rhs.c45=mc(4,5);
  rhs.c46=mc(4,6);

  rhs.c55=mc(5,5);
  rhs.c56=mc(5,6);

  rhs.c66=mc(6,6);
  
  rhs.hasLocal=false;
}

void operator <<(State & rhs, const TrajectoryStateOnSurface & tsos)
{
  rhs<<(*tsos.freeState());
  rhs.hasLocal=true;
  LocalPoint X= tsos.localPosition();
  rhs.x_l=X.x();rhs.y_l=X.y();
  LocalError E= tsos.localError().positionError();
  rhs.v11_l = E.xx();
  rhs.v12_l = E.xy();
  rhs.v22_l = E.yy();
}

void operator<<(Candidate & rhs, TrajectoryStateOnSurface & tsos)
{
  rhs.updated<<tsos;
}


#endif

