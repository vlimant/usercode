#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "UserCode/MStoye/interface/SUSYHelper.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include <TMath.h>
#include "DataFormats/TrackReco/interface/Track.h"

using std::vector;

float  SUSYHelper::SUSYEMETIndirectLepton(const edm::Event& iEvent) const
{
  using namespace edm;
  using reco::TrackCollection;
  Handle<TrackCollection> tracks;
  iEvent.getByLabel("ctfWithMaterialTracks",tracks);
  if(tracks->size()==0) return -2;

  int leadingTrack;

  float apt=0;
  int i=0;
  for(TrackCollection::const_iterator itTrack = tracks->begin();
      itTrack != tracks->end();++itTrack) {
 
    if(apt<itTrack->pt()) {
      apt=itTrack->pt();
      leadingTrack=i;
      // edm::LogInfo("Demo")<<leadingTrack<< " "<< itTrack->pt()<<" "<< apt;
   }
    i++;
  } 
   
  //  edm::LogInfo("Demo")<< i;
 // return (*tracks)[i].pt();
 if((*tracks)[leadingTrack].pt()<15) return -1;

  float ptinR=0;
  for(TrackCollection::const_iterator itTrack = tracks->begin();
      itTrack != tracks->end();++itTrack) {
    // edm::LogInfo("Demo")<< " track i pt:: "<<(*tracks)[leadingTrack].pt()<< " track pt:: "<<itTrack->pt();
 
    float dR =  sqrt( deltaPhi((*tracks)[0].phi(),itTrack->phi())*deltaPhi((*tracks)[0].phi(),itTrack->phi())+((*tracks)[0].eta()-itTrack->eta())*((*tracks)[0].eta()-itTrack->eta()));
    if(dR<0.35)
      {  
	ptinR= ptinR+ itTrack->pt(); 
      }
  }

  return (*tracks)[leadingTrack].pt();
  //  return ptinR/(*tracks)[].pt();
}



float  SUSYHelper::SUSYEMETCorMuons(const pat::MET&  myMET,const vector< pat::Muon>& myMuons  ) const
{
reco::Particle::Vector METmometum(0,0,0);
 
 for( vector< pat::Muon >::const_iterator myMuonIt = myMuons.begin();myMuonIt!=myMuons.end();myMuonIt++)
   {
     METmometum = METmometum + (myMuonIt->momentum());
   }
 
 METmometum = 2*METmometum;
 reco::Particle::Vector final =  myMET.momentum()-METmometum;
 return (sqrt(final.perp2()));
}

float  SUSYHelper::SUSYEChFrack( const vector< pat::Jet >& myJets) const
{
  float aChFrack=0;

  for ( vector< pat::Jet >::const_iterator myJetIt = myJets.begin();myJetIt!=myJets.end();myJetIt++)
    {
      //   if (fabs((*myJetIt).eta())<2.4) continue; // in tracker region!!
      //  if (((*myJetIt).associatedTracks()).size()<4 )continue; 
      aChFrack = aChFrack+ SUSYEChFrack(*myJetIt);
    }
    return aChFrack/float(myJets.size());
}

float  SUSYHelper::SUSYEChFrack( const  pat::Jet &  aJet) const
{
 float trackEt=0;
 const reco::TrackRefVector jetTracks = aJet.associatedTracks();
 for (TrackRefVector::const_iterator jetTrackIt = jetTracks.begin();jetTrackIt!= jetTracks.end();jetTrackIt++)
   {
     trackEt = trackEt + (*jetTrackIt)->pt();
   } 
 return trackEt/aJet.et();
}

float  SUSYHelper::SUSYETEMFrack( const vector< pat::Jet >& myJets) const
{
  float ET_EM_Frack=0;// odd & hardcoded!!
  float myEtSum=0;
  float my_ET_EM_Frack_Sum = 0;
  for ( vector< pat::Jet >::const_iterator myJetIt = myJets.begin();myJetIt!=myJets.end();myJetIt++)
    {
      float et = myJetIt->et();
      my_ET_EM_Frack_Sum =  my_ET_EM_Frack_Sum + ( et* (myJetIt->emEnergyFraction()));
      myEtSum= myEtSum+et;
      ///    edm::LogInfo("Demo")<< " em: "<<myJetIt->emEnergyFraction ();
    }
  ET_EM_Frack = my_ET_EM_Frack_Sum/myEtSum;
  return ET_EM_Frack;
}

// this method is odd but tries to stick to the TDR SUSY Analysis 
float  SUSYHelper::SUSYHT( const vector< pat::Jet >& myJets,const pat::MET&  myMET) const
{
  unsigned int noJets = 4;// odd & hardcoded!!
  float myHT = myMET.et();
  for (unsigned int i = 1;i<noJets;i++)
    {
      // funnily the first jet gets skipped in the TDR  
      if(i<myJets.size())  
	{
	  //  if(myJets[i].recJet().et()>30)
	  //  myHT = myHT + myJets[i].recJet().et();
	  myHT = myHT + myJets[i].et();

	}
    }
  return myHT;
}

float SUSYHelper::SUSYMeff( const vector< pat::Jet >& myJets,const pat::MET&  myMET) const
{
  float myMeff = myMET.et();
  for( vector< pat::Jet >::const_iterator myJetIt = myJets.begin();myJetIt!=myJets.end();myJetIt++)
    {
        myMeff = myMeff + myJetIt->et();
    }
  return myMeff;
}

float SUSYHelper::SUSYMETIso(const vector< pat::Jet >& myJets,const pat::MET&  myMET) const
{
  float myMETIso = 100.0;
  int i =0;
  double deltaPhiAbs =0;
  for( vector< pat::Jet >::const_iterator myJetIt = myJets.begin();myJetIt!=myJets.end();myJetIt++)
    { i++; if(i>4) return myMETIso;

    //  edm::LogInfo("Demo")<< "frac: " << myJetIt->recJet().et() << " " <<myJetIt->et() << " " <<myJetIt->noCorrJet().et()  ;
    //if(fabs(myJetIt->eta())>3||myJetIt->recJet().et()<30) continue;     
     deltaPhiAbs = fabs(deltaPhi(myJetIt->phi(),myMET.phi()));
      if(myMETIso>deltaPhiAbs) myMETIso = deltaPhiAbs;
      
    }
  return myMETIso;
}

float SUSYHelper::SUSYMETRXY(const vector< pat::Jet >& myJets,const pat::MET&  myMET,unsigned int x,unsigned int y) const
{
  if (myJets.size()<3) return -999;
  if (myJets.size()<x||myJets.size()<y){  edm::LogInfo("Demo")<< "SUSYHelper::SUSYMETRXY:: Error ";return -999;}
  double deltax = fabs(deltaPhi(myJets[x].phi(),myMET.phi()));
  double deltay = fabs(deltaPhi(myJets[y].phi(),myMET.phi()))-TMath::Pi();
  float Rxy = sqrt( deltax*deltax+deltay*deltay );
  // LogInfo("Demo")<< " Rx "<<Rxy;
  return Rxy;
}

reco::Particle::Vector SUSYHelper::SUSYAllRecoilMET(const vector< pat::Jet >& myJets,const vector< pat::Electron >& myElectrons,const vector< pat::Muon>& myMuons ,const vector< pat::Tau>&  myTaus) const
{
reco::Particle::Vector METmometum(0,0,0);
for( vector< pat::Jet >::const_iterator myJetIt = myJets.begin();myJetIt!=myJets.end();myJetIt++)
    {
      METmometum = METmometum + myJetIt->momentum();
    }

for( vector< pat::Muon >::const_iterator myMuonIt = myMuons.begin();myMuonIt!=myMuons.end();myMuonIt++)
    {
      METmometum = METmometum + myMuonIt->momentum();
    }

/*for( vector< pat::Electron >::const_iterator myElectronIt = myElectrons.begin();myElectronIt!=myElectrons.end();myElectronIt++)
    {
      METmometum = METmometum + myElectronIt->momentum();
    }

for( vector< pat::Tau >::const_iterator myTauIt = myTaus.begin();myTauIt!=myTaus.end();myTauIt++)
    {
      METmometum = METmometum + myTauIt->momentum();
    }
*/
return -METmometum;
}

float SUSYHelper::SUSYAllGenRecoilMET(const vector< pat::Jet >& myJets,const vector< pat::Muon>& myMuons ,const pat::MET&  myMET) const
{
  reco::Particle::Vector METmometum(0,0,0);
  for( vector< pat::Jet >::const_iterator myJetIt = myJets.begin();myJetIt!=myJets.end();myJetIt++)
    {
      METmometum = METmometum + myJetIt->momentum() - myJetIt->genJet()->momentum();
    }
  
  for( vector< pat::Muon >::const_iterator myMuonIt = myMuons.begin();myMuonIt!=myMuons.end();myMuonIt++)
    {
      METmometum = METmometum +  myMuonIt->momentum() - myMuonIt->genLepton()->momentum();
    }
  const reco::Particle* myGenMET = myMET.genMET();

  METmometum  = METmometum + myGenMET->momentum();
  return sqrt(METmometum.perp2());
  
  /*for( vector< pat::Electron >::const_iterator myElectronIt = myElectrons.begin();myElectronIt!=myElectrons.end();myElectronIt++)
    {
      METmometum = METmometum + myElectronIt->momentum();
    }

for( vector< pat::Tau >::const_iterator myTauIt = myTaus.begin();myTauIt!=myTaus.end();myTauIt++)
    {
      METmometum = METmometum + myTauIt->momentum();
    }
*/
}

reco::Particle::Vector SUSYHelper::MuonCorrection( reco::Particle::Vector METmometum, const vector< pat::Muon>& myMuons) const
{

 for( vector< pat::Muon >::const_iterator myMuonIt = myMuons.begin();myMuonIt!=myMuons.end();myMuonIt++)
    {
      METmometum = METmometum +  myMuonIt->momentum() - myMuonIt->genLepton()->momentum();
    }
 return METmometum;
}

reco::Particle::Vector SUSYHelper::SUSYRecoilMET(const vector< pat::Jet >& myJets) const
{ reco::Particle::Vector METmometum(0,0,0);
    
  
  for( vector< pat::Jet >::const_iterator myJetIt = myJets.begin();myJetIt!=myJets.end();myJetIt++)
    {
      METmometum = METmometum + myJetIt->momentum();
      //    edm::LogInfo("Demo") << "METmometum " << sqrt(METmometum.perp2())<< " pt jet" <<myJetIt->et() <<METmometum ;
    }
  //  edm::LogInfo("Demo") << "METmometum met " << METmometum.r();
  return -METmometum;
}

reco::Particle::Vector SUSYHelper::SUSYRecoilMETCutted(const vector< pat::Jet >& myJets, float ptMin = 30, float etaMax = 3) const
{ 
  reco::Particle::Vector METmometum(0,0,0);
  for( vector< pat::Jet >::const_iterator myJetIt = myJets.begin();myJetIt!=myJets.end();myJetIt++)
    {
      // some cuts to avoid fake jets
      if ( myJetIt->et()>ptMin && fabs(myJetIt->eta())<etaMax )   METmometum = METmometum + myJetIt->momentum();
    }
  return -METmometum;
}
