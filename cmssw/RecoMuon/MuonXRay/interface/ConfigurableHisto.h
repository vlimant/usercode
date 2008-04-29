class ConfigurableAxis {
 public:
  ConfigurableAxis(){}
  ConfigurableAxis(edm::ParameterSet par){
    if (par.exists("Label"))
      Label_=par.getParameter<std::string>("Label");
    else Label_="";

    if (par.exists("nBins")){
      nBin_=par.getParameter<uint>("nBins");
      Min_=par.getParameter<double>("Min");
      Max_=par.getParameter<double>("Max");
    }else{
      vBins_=par.getParameter<std::vector<double> >("vBins");
    }
  }
  bool variableSize(){return vBins_.size()!=0;}
  uint nBin(){if (variableSize()) return vBins_.size()-1;else return nBin_;}
  double Min(){if (variableSize()) return vBins_.front(); else return Min_;}
  double Max(){if (variableSize()) return vBins_.back(); else return Max_;}
  const std::string & Label(){ return Label_;}

 private:
  std::vector<double> vBins_;
  uint nBin_;
  double Min_;
  double Max_;
  std::string Label_;
};

class ConfigurableH1 {
 public:
  ConfigurableH1(edm::ParameterSet & par, DaqMonitorBEInterface*dbe){
    name_=par.getParameter<std::string>("name");
    std::string title=par.getParameter<std::string>("title");
    ConfigurableAxis a;
    //one does not have to specify PSet xAxis =... for TH1
    if (par.exists("xAxis")){      a=ConfigurableAxis(par.getParameter<edm::ParameterSet>("xAxis"));
    }else{       a=ConfigurableAxis(par);}

    me_=dbe->book1D(name_,title,a.nBin(),a.Min(),a.Max());
    //set axis label
    TH1 * h=MEA1<TH1>(me_);
    h->GetXaxis()->SetTitle(a.Label().c_str());
    std::string yLabel=par.getUntrackedParameter<std::string>("yLabel",std::string());
    h->GetYaxis()->SetTitle(yLabel.c_str());
  }
  const std::string & name() { return name_;}
  MonitorElement * me() { return me_;}

 private:
  std::string name_;
  MonitorElement * me_;

};

class ConfigurableH2 {
 public:
  ConfigurableH2(edm::ParameterSet & par, DaqMonitorBEInterface*dbe){
    name_=par.getParameter<std::string>("name");
    std::string title=par.getParameter<std::string>("title");
    ConfigurableAxis ax(par.getParameter<edm::ParameterSet>("xAxis"));
    ConfigurableAxis ay(par.getParameter<edm::ParameterSet>("yAxis"));
    me_=dbe->book2D(name_,title,
		    ax.nBin(),ax.Min(),ax.Max(),
		    ay.nBin(),ay.Min(),ay.Max());
    //set axis label
    TH2 * h=MEA2<TH2>(me_);
    h->GetXaxis()->SetTitle(ax.Label().c_str());
    h->GetYaxis()->SetTitle(ay.Label().c_str());
  }

  const std::string & name() { return name_;}
  MonitorElement * me() { return me_;}
 private :
  std::string name_;
  MonitorElement * me_;
};

/*
class ConfigurableStackH1 {
 public:
  ConfigurableStackH1(edm::ParameterSet & par){
    identifier=par.getParameter<std::string>("identifier");

    std::string title=par.getParameter<std::string>("title");
    uint nBinsX=par.getParameter<uint>("nBinsX");
    double xMin=par.getParameter<double>("xMin");
    double xMax=par.getParameter<double>("xMax");
    uint nBinxY=par.getParameter<uint>("nBinsY");
    double yMin=par.getParameter<double>("yMin");
    double yMax=par.getParameter<double>("yMax");
    
  }
  void push_back(std::string specificLabel){
    me_=dbe->Book1D
  }
  std::string identifier;
  std::vector<MonitorElement *> me
};
*/


class ConfigurableHistogram {
 public:
  enum TYPE { H1, H2 };
  ConfigurableHistogram(edm::ParameterSet & par, DaqMonitorBEInterface*dbe){
    bool returnToDir=false;
    std::string oldDir=dbe->pwd();
    if (par.exists("directory")){
      dbe->setCurrentFolder(par.getParameter<std::string>("directory"));
      returnToDir=true;
    }else if (par.exists("subdirectory")){
      dbe->setCurrentFolder(oldDir+std::string("/")+par.getParameter<std::string>("subdirectory"));
      returnToDir=true;
    }
      
	
    name_=par.getParameter<std::string>("name");
    std::string title=par.getParameter<std::string>("title");
    
    //the presence of PSet yAxis={... tells it is a TH2 or TProfile

    if (!par.exists("yAxis")){
      //TH1
      type_=H1;
      ConfigurableAxis a;
      //one does not have to specify PSet xAxis =... for TH1
      if (par.exists("xAxis")){      a=ConfigurableAxis(par.getParameter<edm::ParameterSet>("xAxis"));
      }else{       a=ConfigurableAxis(par);}
      
      me_=dbe->book1D(name_,title,a.nBin(),a.Min(),a.Max());
      //set axis label
      TH1 * h=MEA1<TH1>(me_);
      h->GetXaxis()->SetTitle(a.Label().c_str());
      std::string yLabel=par.getUntrackedParameter<std::string>("yLabel",std::string());
      h->GetYaxis()->SetTitle(yLabel.c_str());
    }else{
      //TH2
      type_=H2;
      ConfigurableAxis ax(par.getParameter<edm::ParameterSet>("xAxis"));
      ConfigurableAxis ay(par.getParameter<edm::ParameterSet>("yAxis"));
      me_=dbe->book2D(name_,title,
		      ax.nBin(),ax.Min(),ax.Max(),
		      ay.nBin(),ay.Min(),ay.Max());
      //set axis label
      TH2 * h=MEA2<TH2>(me_);
      h->GetXaxis()->SetTitle(ax.Label().c_str());
      h->GetYaxis()->SetTitle(ay.Label().c_str());
    }

    if (returnToDir){
      dbe->setCurrentFolder(oldDir);
    }
  }

  const std::string & name() { return name_;}
  MonitorElement * me() { return me_;}
  bool isH1(){ return type_==H1;}
  bool isH2(){return type_==H2;}
  TYPE type(){ return type_;}
 private :
  TYPE type_;
  std::string name_;
  MonitorElement * me_;
};
