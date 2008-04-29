#include "FWCore/ParameterSet/interface/ParameterSet.h"

class DQMHelper;

class IDconverttoBinNum {
public:
  //constructor by hand
  IDconverttoBinNum();
  //construtor from a PSet
  IDconverttoBinNum(edm::ParameterSet pset);
  
  int GetBinNum(int pdgID);
  std::string GetBinName(uint BinNumvalue);

  //total size: bins + overflow
  int size() { return Bins.size()+1;}

  //overflow bin
  int maxIndex() { return Bins.size();}

  //and its label
  std::string other(){ return "other";}

  void splitH1ByCategory(std::string identifier, DQMHelper & h, const char * format);
private:
  typedef std::pair< std::string, std::vector<int> > bin;
  std::vector<bin> Bins;
  //in which bin a pdgID goes
  std::map<int, uint> map;
};
