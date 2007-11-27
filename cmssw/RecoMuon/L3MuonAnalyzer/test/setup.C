void setup(char * rfile)
{
  gROOT->ProcessLine(".L /uscms/home/vlimant/work/CMSSW_1_3_1/src/RecoMuon/L3MuonAnalyzer/test/Monitor.C+");
  Monitor(UAF);
  gROOT->ProcessLine(".L /uscms/home/vlimant/work/CMSSW_1_3_1/src/RecoMuon/L3MuonAnalyzer/test/Display.C");
  gROOT->ProcessLine(".L /uscms/home/vlimant/work/CMSSW_1_3_1/src/RecoMuon/L3MuonAnalyzer/test/ShowMe.C+");

  open(rfile);

}
