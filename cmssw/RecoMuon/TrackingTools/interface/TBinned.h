#include "TAxis.h"
#include "TString.h"
#include <sstream>
#include <iostream>
#include "FWCore/MessageLogger/interface/MessageLogger.h"

template <typename O> class TBinned : public TNamed {

  public:
  TBinned(const uint N, const TString & name,const TString & title ): TNamed(name,title), Naxis(N) /*,iCall(0)*/
    {
      axis = std::vector<TAxis*>(Naxis, 0);
      //      bins = std::vector<uint>(Naxis, 0);
      withOverUnderFlow=false;
      tellIndex=false;
    }


    class Node{
    public:
      Node():who(0),iCall(0){}
      Node(TBinned<O> & w): who(&w), iCall(0), bins(w.Naxis,0){}
      Node(TBinned<O> & w, uint i, std::vector<uint> & b): who(&w), iCall(i), bins(b){}

      TBinned<O> * who;
      uint iCall;
      std::vector<uint> bins;
      
      Node operator[] (double x){
	if (iCall==bins.size()){
	  std::cout<<"error"<<std::endl;
	  return Node();
	}
	else{
	  Node n(*who,iCall,bins);
	  n.bins[n.iCall]=n.who->findBin(x, n.who->axis[n.iCall]);
	  n.iCall++;
	  return n;
	}
      }

      O & operator()() {
	if (iCall%who->Naxis==0){
	  iCall=0; 
	  uint ec = who->encode(bins);
	  O & r=who->data[ec];
	  for (uint i=0;i!=who->Naxis;i++) {bins[i]=0;}
	  return r;
	}
	else{
	  edm::LogError("TBinned")<<"error, not enough [] operators. ("<<iCall<<")";
	  std::cout<<"error, not enough [] operators. ("<<iCall<<")"<<std::endl;
	  iCall=0;
	  return who->data[0];
	}
      }

      void getContent( ostream & o){
	if (iCall%who->Naxis==0){
	  iCall=0; 
	  uint ec = who->encode(bins);
	  O & r=who->data[ec];
	  o<<"value in bin: ";
	  for (uint i=0;i!=who->Naxis;i++){
	    o<<who->bin(who->axis[i], bins[i]);
	    if (i<who->Naxis-1) o<<", ";
	  }
	  o<<"= "<<r<<std::endl;
	  
	  for (uint i=0;i!=who->Naxis;i++) {bins[i]=0;}
	}
	else{
	  //	edm::LogError("TBinned")<<"error, not enough [] operators. ("<<iCall<<")";
	  std::cout<<"error, not enough [] operators. ("<<iCall<<")"<<std::endl;
	}      
      }

    };
    
    void setWithUnderOverFlow(bool b=true) { withOverUnderFlow=b;}
    void setTellIndex(bool b=true) { tellIndex=b;}

    uint getNbins(TAxis * a){
      if (withOverUnderFlow) return a->GetNbins()+2;
      else return a->GetNbins();
    }

    //find the right bin, taking care of over/under flow
    uint findBin(double x, TAxis * a){
      if (withOverUnderFlow) return a->FindBin(x); //0: underflow, 1..Nbins: bins, Nbins+1: overflow
      else{
	int fb = a->FindBin(x);
	if (fb==0) return 0;//0: under and first bin
	else if (fb==a->GetNbins()+1) return static_cast<uint>(a->GetNbins()-1);//Nbins-1: last and over flow
	else return static_cast<uint>(fb-1);//regular bin
      }
    }
      
    //access the different bins.
    //    needs  Naxis [] operator, then one () operator to access the actual value.
    Node operator[] (const double x){
      Node n(*this);
      return n[x];
    }

    /*
    TBinned & operator[] (const double x){
      if (iCall==Naxis) { 
	edm::LogError("TBinned")<<"error, not the right amount of []. "<<iCall<<" so far.";
	std::cout<<"error, not the right amount of []. "<<iCall<<" so far."<<std::endl;
	iCall=0;
      }
      else {
	bins[iCall] = findBin(x, axis[iCall]);
	iCall++;
      }
      return (*this);
    }
    */
    
    //get the value when sufficient number of [] operator where put in front
    /*    O & operator()() {
	  if (iCall%Naxis==0){
	iCall=0; 
	uint ec = encode(bins);
	O & r=data[ec];
	for (uint i=0;i!=Naxis;i++) {bins[i]=0;}
	return r;
      }
      else{
	edm::LogError("TBinned")<<"error, not enough [] operators. ("<<iCall<<")";
	std::cout<<"error, not enough [] operators. ("<<iCall<<")"<<std::endl;
	iCall=0;
	return data[0];
      }
    }
    */

    // define an axis from fixed binning.
    void setAxis(uint i,const TString title, uint nBin, double * xbins){
      if (i>=Naxis) { 
	edm::LogError("TBinned")<<"in setting axis: "<<i<<", dimension: "<<Naxis<<" allowed.";
	Naxis=i+1;
	resize();
      }
      if (axis[i]) { edm::LogError("TBinned")<<"redefining axis: "<<i; delete axis[i];}
      std::stringstream ss;
      for (uint j=0;j!=nBin+1;j++) { ss<<xbins[j];if (j<nBin) ss<<", "; }
      edm::LogInfo("TBinned")<<"setting axis: "<<i<<" "<<title<<" "<<nBin<<" variable bins: ["<<ss.str()<<"]";
      std::cout<<"setting axis: "<<i<<" "<<title<<" "<<nBin<<" variable bins: ["<<ss.str()<<"]"<<std::endl;
      axis[i] = new TAxis(nBin, xbins);
      axis[i]->SetTitle(title);
    }
    
    //define an axis from variable binning
    void setAxis(uint i, const TString title, uint nBin, double xmin, double xmax){
      if (i>=Naxis) { 
	edm::LogError("TBinned")<<"in setting axis: "<<i<<", dimension: "<<Naxis<<" allowed.";
	Naxis=i+1;
	resize();
      }
      if (axis[i]) {  edm::LogError("TBinned")<<"redefining axis: "<<i; delete axis[i];}
      edm::LogInfo("TBinned")<<"setting axis: "<<i<<" "<<title<<" "<<nBin<<" "<<xmin<<" "<<xmax;
      std::cout<<"setting axis: "<<i<<" "<<title<<" "<<nBin<<" "<<xmin<<" "<<xmax<<std::endl;
      axis[i] = new TAxis(nBin, xmin, xmax);
      axis[i]->SetTitle(title);
    }
    
    // check axis consistency and allocate memory
    bool check(){
      if (axis.size()!=Naxis) return false;
      //      if (bins.size()!=Naxis) return false;
      for (uint i=0;i!=Naxis;i++) {if (!axis[i]) return false;} return reserve();
    }

    bool resize(){
      axis.resize(Naxis, 0);
      //      bins.resize(Naxis, 0);
      return true;
    }

    // allocate memory to the object
    bool reserve() {
      data.clear();
      //calculate max index
      uint maxIndex=1;
      for (uint i=0;i!=Naxis;i++)  maxIndex*=getNbins(axis[i]);
      //allocate sufficient memory to data field
      dataSize = maxIndex;
      data.reserve(maxIndex);
      return true;
    }

    // clear memory for the OBJECT
    //    bool clearDelete(){ for (uint i=0;i!=dataSize;i++){  delete data[i];} ; return true;}
    
    // goes from arry of iBin to single index
    uint encode( std::vector<uint> array){
      if (array.size()!=Naxis) { std::cout <<"error. wrong size array."<<std::endl; return 0;}
      std::stringstream ss;
      ss<<" array[0]: "<<array[0];
      uint index=array[0];
      for (uint i=1; i<Naxis;i++){
	index*= getNbins(axis[i]);
	index+= array[i];
	ss<<" array["<<i<<"]: "<<array[i];
      }
      ss<<" encoded value: "<<index;
      if (tellIndex){
	edm::LogInfo("TBinned")<<ss.str();
	std::cout<<ss.str()<<std::endl;
      }
      return index; 
    }
    
    // goes from single index to array of iBin.
    void decode(uint index, std::vector<uint> & array){
      if (array.size()!=Naxis) { std::cout <<"error. wrong size array."<<std::endl; return;}
      std::stringstream ss;
      ss<<"global index: "<<index;
      for (uint iV=0;iV!=Naxis;iV++){
	uint i=Naxis-iV-1;
	array[i]=index%getNbins(axis[i]);
	index-=array[i];
	index/=getNbins(axis[i]);
	ss<<" array["<<i<<"]: "<<array[i];
      }
      if (tellIndex){
	edm::LogInfo("TBinned")<<ss.str();
	std::cout<<ss.str()<<std::endl;
      }
    }
    
    double center(TAxis *a ,uint i) {
      if(withOverUnderFlow){ 
	return a->GetBinCenter(i);
      }
      else{
	return a->GetBinCenter(i+1);
      }
    }

    TString bin(TAxis *a ,uint i)  {
      if (withOverUnderFlow){
	if (i==0) return TString(Form("%s: [underflow, %5.2f]", a->GetTitle(),a->GetBinUpEdge(i)));
	else if (i==static_cast<uint>(a->GetNbins()+1)) return TString(Form("%s: [%5.2f, overflow]", a->GetTitle(),a->GetBinLowEdge(i)));
	else return TString(Form("%s: [%5.2f, %5.2f]", a->GetTitle(), a->GetBinLowEdge(i), a->GetBinUpEdge(i)));
      }
      else{
	return TString(Form("%s: [%5.2f, %5.2f]", a->GetTitle(), a->GetBinLowEdge(i+1), a->GetBinUpEdge(i+1))); 
      }
    }
    
    /*    void getContent( ostream & o){
	  if (iCall%Naxis==0){
	  iCall=0; 
	  uint ec = encode(bins);
	  O & r=data[ec];
	  o<<"value in bin: ";
	  for (uint i=0;i!=Naxis;i++){
	  o<<bin(axis[i], bins[i]);
	  if (i<Naxis-1) o<<", ";
	  }
	  o<<"= "<<r<<std::endl;

	  for (uint i=0;i!=Naxis;i++) {bins[i]=0;}
	  }
      else{
      //	edm::LogError("TBinned")<<"error, not enough [] operators. ("<<iCall<<")";
      std::cout<<"error, not enough [] operators. ("<<iCall<<")"<<std::endl;
      }      
      }
    */

    //show the content of the TBinned object
    void print( ostream & o){
      std::vector<uint> b(Naxis, 0);
      //loop over the global index
      for (uint iI=0;iI!=dataSize;iI++){
	//decode the global index into bin numbers
	decode(iI, b);
	//get the label of each bin
	for (uint iA=0;iA!=Naxis;iA++){
  	  o<<bin(axis[iA], b[iA]);
	  if (iA<Naxis-1) o<<", ";
	}
	//get the value for that global bin
	o<<"= "<<data[iI]<<std::endl;
      }
    }

    //dump is a root friendly file with bin center as values
    void dump( ostream & o){
      std::vector<uint> b(Naxis, 0);
      //loop over the global index
      for (uint iI=0;iI!=dataSize;iI++){
	//decode the global index into bin numbers
	decode(iI, b);
	//get the label of each bin
	for (uint iA=0;iA!=Naxis;iA++){ o<< center(axis[iA], b[iA])<<" ";}
	o<<data[iI]<<std::endl;
      }
      
    }
  private:
    //number of axis
    uint Naxis;

    // list of axis
    std::vector<TAxis*> axis;

    // internal index for accessor
    //    mutable uint iCall;
    // internal array of iBin for accessor
    //    mutable std::vector<uint> bins;

    // dynamic arry of data
    uint dataSize;
    std::vector<O> data;

    bool withOverUnderFlow;
    bool tellIndex;
};



