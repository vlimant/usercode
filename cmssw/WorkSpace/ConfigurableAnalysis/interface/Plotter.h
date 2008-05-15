#ifndef ConfigurableAnalysis_Plotter_H
#define ConfigurableAnalysis_Plotter_H

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/UtilAlgos/interface/TFileDirectory.h"

#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"

#include "Workspace/ConfigurableAnalysis/interface/VariableHelper.h"
#include "Workspace/ConfigurableAnalysis/interface/CachingVariable.h"
#include "Workspace/ConfigurableAnalysis/interface/ConfigurableHisto.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

class Plotter {
 public:
  Plotter(edm::ParameterSet iConfig) : currentDir_("youDidNotSetDirectoryFirst") {
    //create the master copy, never filled, just to make copies

    //    make TH1
    edm::ParameterSet th1=iConfig.getParameter<edm::ParameterSet>("TH1s");
    std::vector<std::string> th1Names;
    th1.getParameterSetNames(th1Names);
    for (uint iH=0;iH!=th1Names.size();++iH){
      std::string hname = th1Names[iH];
      edm::ParameterSet hPset=th1.getParameter<edm::ParameterSet>(hname);
      bool split=hPset.exists("splitter") || hPset.exists("splitters");
      if (split)
	master_[hname]=new SplittingConfigurableHisto(ConfigurableHisto::h1, hname, hPset);
      else
	master_[hname]=new ConfigurableHisto(ConfigurableHisto::h1, hname, hPset);
    }

    //    make profiles
    edm::ParameterSet tprof=iConfig.getParameter<edm::ParameterSet>("TProfiles");
    std::vector<std::string> tprofNames;
    tprof.getParameterSetNames(tprofNames);
    for (uint iH=0;iH!=tprofNames.size();++iH){
      std::string hname = tprofNames[iH];
      edm::ParameterSet hPset=tprof.getParameter<edm::ParameterSet>(hname);
      bool split=hPset.exists("splitter") || hPset.exists("splitters");
      if (split)
	master_[hname]=new SplittingConfigurableHisto(ConfigurableHisto::prof, hname, hPset);
      else
	master_[hname]=new ConfigurableHisto(ConfigurableHisto::prof, hname, hPset);
    }
    
    //    make TH2
    edm::ParameterSet th2=iConfig.getParameter<edm::ParameterSet>("TH2s");
    std::vector<std::string> th2Names;
    th2.getParameterSetNames(th2Names);
    for (uint iH=0;iH!=th2Names.size();++iH){
      std::string hname = th2Names[iH];
      edm::ParameterSet hPset=th2.getParameter<edm::ParameterSet>(hname);
      bool split=hPset.exists("splitter") || hPset.exists("splitters");
      if (split)
	master_[hname]=new SplittingConfigurableHisto(ConfigurableHisto::h2, hname, hPset);
      else
	master_[hname]=new ConfigurableHisto(ConfigurableHisto::h2, hname, hPset);
    }
  }

  void setDir(std::string dir){
    //insert a new one
    Directory & insertedDirectory = directories_[dir];

    //create the actual directory in TFile: name is <dir>
    insertedDirectory.dir=new TFileDirectory(edm::Service<TFileService>()->mkdir(dir));
    insertedDirectory.dirName=dir;

    //make a dummy histo???
    //    insertedDirectory.dir->make<TH1F>((dir+"dummy").c_str(),"dummy",10,0.,1.);
    //    insertedDirectory.dir->mkdir("testSubDir");

    //remember which directory name this is
    currentDir_=dir;
  }
  
  void fill(std::string subDir,const edm::Event& iEvent){
    //what is the current directory
    Directory & currentDirectory= directories_[currentDir_];

    //what is the current set of sub directories for this
    SubDirectories & currentSetOfSubDirectories=currentDirectory.subDir;
    
    //find the subDirectory requested:
    SubDirectory * subDirectoryToUse=0;
    SubDirectories::iterator subDirectoryFindIterator=currentSetOfSubDirectories.find(subDir);

    //not found? insert a new directory with this name
    if (subDirectoryFindIterator==currentSetOfSubDirectories.end()){
      SubDirectory & insertedDir = currentSetOfSubDirectories[subDir];
      subDirectoryToUse = &insertedDir;
      insertedDir.dir=new TFileDirectory(currentDirectory.dir->mkdir(subDir));
      insertedDir.dirName=subDir;

      //create a copy from the master copy
      DirectoryHistos::iterator masterHistogramIterator=master_.begin();
      DirectoryHistos::iterator masterHistogramIterator_end=master_.end();
      for (; masterHistogramIterator!=masterHistogramIterator_end;++masterHistogramIterator)
	{
	  //clone does not book histogram
	  insertedDir.histos[masterHistogramIterator->first]=masterHistogramIterator->second->clone();
	}
      
      //book all copies of the histos
      DirectoryHistos::iterator clonedHistogramIterator=insertedDir.histos.begin();
      DirectoryHistos::iterator clonedHistogramIterator_end=insertedDir.histos.end();
      for (; clonedHistogramIterator!=clonedHistogramIterator_end;++clonedHistogramIterator)
	{
	  clonedHistogramIterator->second->book(insertedDir.dir);
	}
    }
    else{
      subDirectoryToUse=&subDirectoryFindIterator->second;
    }
    
    //now that you have the subdirectory: fill histograms for this sub directory
    DirectoryHistos::iterator histogramIterator=subDirectoryToUse->histos.begin();
    DirectoryHistos::iterator histogramIterator_end=subDirectoryToUse->histos.end();
    for(; histogramIterator!=histogramIterator_end;++histogramIterator)
      { histogramIterator->second->fill(iEvent); }
  }
  
 private:
  typedef std::map<std::string, ConfigurableHisto *> DirectoryHistos;
  DirectoryHistos master_;

  class SubDirectory {
  public:
    SubDirectory() : dir(0){}
    //    SubDirectory(TFileDirectory d) : dir(new TFileDirectory(d)){}
    std::string dirName;
    DirectoryHistos histos;
    TFileDirectory * dir;
  };
  typedef std::map<std::string, SubDirectory> SubDirectories;

  class Directory {
  public:
    Directory() : dir(0){}
    std::string dirName;
    SubDirectories subDir;
    TFileDirectory * dir;
  };
  typedef std::map<std::string, Directory> Directories;

  std::string currentDir_;
  Directories directories_;
};

#endif
