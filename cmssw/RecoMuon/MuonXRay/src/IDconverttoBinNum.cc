#include "RecoMuon/MuonXRay/interface/IDconverttoBinNum.h"  
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "RecoMuon/MuonXRay/interface/DQMHelper.h"

void IDconverttoBinNum::splitH1ByCategory(const std::string identifier, DQMHelper & h, const char * format){
  const std::string id=identifier+"_assoc_ID";
  h.HBC[id].resize(this->size());
  std::string htitle,hname,binname;
  TH1F* h1 = h.h1<TH1F>(identifier);
  for(int iut = 0;iut<this->size();iut++)
    {
      binname = this->GetBinName(iut);
      htitle = Form(format, h1->GetTitle(), binname.c_str());
      hname = Form("%s_%i",id.c_str(),iut);
      //fixme. what about variable sized axis
      h.book1D(iut,id, hname,htitle,
	       h1->GetXaxis()->GetNbins(),
	       h1->GetXaxis()->GetXmin(),
	       h1->GetXaxis()->GetXmax(),
	       h1->GetXaxis()->GetTitle(),
	       h1->GetYaxis()->GetTitle());
    }
}


IDconverttoBinNum::IDconverttoBinNum(edm::ParameterSet pset){
  std::vector<edm::ParameterSet> ranges = pset.getParameter<std::vector<edm::ParameterSet> > ("ranges");

  Bins.resize(ranges.size());
  for (uint r=0; r!=ranges.size();++r){
    Bins[r].first = ranges[r].getParameter<std::string>("label"); 
    Bins[r].second = ranges[r].getParameter<std::vector<int> >("pdgIDs");
    for (std::vector<int>::iterator i = Bins[r].second.begin(); i!=Bins[r].second.end();++i){
      map[*i]=r;
    }
  }

}

IDconverttoBinNum::IDconverttoBinNum()
{
  //configured by hand as Finn did it
  Bins.resize(14);

  Bins[0].first = "ID=0";
  Bins[0].second.push_back(0);

  Bins[1].first = "#pi+/-";
  Bins[1].second.push_back(211);
  Bins[1].second.push_back(-211);

  Bins[2].first = "K^{+/-}";
  Bins[2].second.push_back(321);
  Bins[3].second.push_back(-321);

  Bins[3].first = "#tau";
  Bins[3].second.push_back(15);
  Bins[3].second.push_back(-15);

  Bins[4].first = "D^{+/-}";
  Bins[4].second.push_back(411);
  Bins[4].second.push_back(-411);

  Bins[5].first = "D^{0}";
  Bins[5].second.push_back(421);
  Bins[5].second.push_back(-421);

  Bins[6].first = "B^{+}";
  Bins[6].second.push_back(521);
  Bins[6].second.push_back(-521);

  Bins[7].first = "B^{0";
  Bins[7].second.push_back(511);
  Bins[7].second.push_back(-511);

  Bins[8].first = "B_{s}";
  Bins[8].second.push_back(531);
  Bins[8].second.push_back(-531);

  Bins[9].first = "#Lambda_{b}";
  Bins[9].second.push_back(5122);
  Bins[9].second.push_back(-5122);

  Bins[10].first = "J/#Psi";
  Bins[10].second.push_back(443);
  Bins[10].second.push_back(-443);

  Bins[11].first = "#Upsilon(nS)";
  Bins[11].second.push_back(553);
  Bins[11].second.push_back(-553);
  Bins[11].second.push_back(200553);
  Bins[11].second.push_back(-200553);
  Bins[11].second.push_back(100553);
  Bins[11].second.push_back(-100553);

  Bins[12].first = "W^{+/-}";
  Bins[12].second.push_back(24);
  Bins[12].second.push_back(-24);

  Bins[13].first = "Z^{0}";
  Bins[13].second.push_back(23);

  //and filling the map
  for (uint r=0;r!=Bins.size();++r){
    for (uint i=0;i!=Bins[r].second.size();++i)
      {map[Bins[r].second[i]] = r;}}
}


int IDconverttoBinNum::GetBinNum(int pdgID)
{
  std::map<int, uint>::iterator i=map.find(pdgID);
  //if registerd pdgid, return the bin
  if ( i!= map.end()){return i->second;}
  //else return overflow bin
  else {
    edm::LogError("IDconverttoBinNum")<<pdgID<<" is not a registered pdgID.";
    return maxIndex();}
}

std::string IDconverttoBinNum::GetBinName(uint BinNumvalue)
{
  if (BinNumvalue>Bins.size()) return "bin number is not valid";
  else if (BinNumvalue==Bins.size()) return other();
  else { return Bins[BinNumvalue].first; }
}
