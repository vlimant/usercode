{
  gROOT->ProcessLine(".L Monitor.C");
  Monitor(UAF);
  gROOT->ProcessLine(".L Display.C");
  //  open("../../Run/muRSanalyzer-test.root") ;
  //    open("../../Run/ckfPIXELanalyzer-test.root") ;
      open("../../Run/ckfRSanalyzer-test.root") ;
    
  //  Display(0,0);

  //  Display("/uscms/home/vlimant/nobackup/Temp/evDisplay/display_");
}
