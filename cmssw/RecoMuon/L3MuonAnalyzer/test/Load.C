{
  //  TString workArea="/afs/fnal.gov/files/home/room2/vlimant/work/CMSSW_1_2_0_pre9/src";
  //  TString workArea="${CMSSW_BASE}/src/";

  TString workArea = gSystem->ExpandPathName("$CMSSW_BASE/");
  TString classCode="RecoMuon/L3MuonAnalyzer/src/DumpClass.cc";

  /*
    TString addToPath="-I"+workArea+"src/";
    gSystem->AddIncludePath(addToPath);
    addToPath="-I"+TString(gSystem->ExpandPathName("$CMSSW_RELEASE_BASE/src/"));
    gSystem->AddIncludePath(addToPath);
    addToPath="-I"+TString(gSystem->ExpandPathName("$CMSSW_PATH/$SCRAM_ARCH/external/"));
    gSystem->AddIncludePath(addToPath);
  */

  //compile the class
  gROOT->ProcessLine(".L "+workArea+"src/"+classCode+"+");

  //the libreary directly
  //  gROOT->ProcessLine(".L "+workArea+"lib/slc3_ia32_gcc323/libRecoMuonL3MuonAnalyzer.so");
}

