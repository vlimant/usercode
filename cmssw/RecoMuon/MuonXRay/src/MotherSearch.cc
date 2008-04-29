#include "RecoMuon/MuonXRay/interface/MotherSearch.h"

MotherSearch::MotherSearch(const SimTrack * isimtk,
			   edm::Handle<edm::SimTrackContainer> & SimTk,
			   edm::Handle<edm::SimVertexContainer> & SimVtx,
			   edm::Handle<edm::HepMCProduct> & hepmc):
  Sim_vertex(0), Sim_mother(0), Gen_vertex(0), Gen_mother(0) 
{
  const std::string theCategory = "MotherSearch::MotherSearch";
  
  simtrack = isimtk;
  const HepLorentzVector momentum = isimtk->momentum();
  int selfID = isimtk->type();

  useGen=false;
  
  if(isimtk->vertIndex()>=0 && isimtk->vertIndex()<SimVtx->size()){
    //access the vertex
    Sim_vertex = &(*SimVtx)[isimtk->vertIndex()];
    if (!Sim_vertex->noParent()){
      LogDebug(theCategory)<<"I am here 3";
      
      const HepLorentzVector position = Sim_vertex->position();
      edm::LogVerbatim(theCategory)<<"This sim vertex position is :"<<position.v()
				   <<"\nThis sim vertex magnitude is :"<<position.v().mag()
				   <<"\nThis sim vertex magnitude is :"<<position.mag();
      
      LogDebug(theCategory)<<"I am here 4";      
      
      int parentTrkNum = Sim_vertex->parentIndex();
      edm::LogVerbatim(theCategory)<<"This track ID is: "<<isimtk->type()
				   <<"\nThe mother track number is: "<<parentTrkNum; 
      int num_matches = 0; 
      LogDebug(theCategory)<<"I am here 5";
      //now loop the simtracks and find the parent
      for (std::vector<SimTrack>::const_iterator isimtk_parent = SimTk->begin();isimtk_parent!=SimTk->end();++isimtk_parent){
	if(isimtk_parent->trackId()==parentTrkNum){
	  num_matches++;
	  int parentID = isimtk_parent->type();
	  if(abs(parentID)!=selfID){
	    edm::LogError(theCategory)<<"The mother ID is: "<<isimtk_parent->type();
	    Sim_mother = &(*isimtk_parent);
	  }
	  else{edm::LogError(theCategory)<<"The mother of the SimTrack is a SimTrack of same type. skipping";break;}
	  LogDebug(theCategory)<<"I am here 6";      
	}//matching trackId
      }//second SimTrack loop
    }
    else {//the SimVertex has no parent
      useGen=true;
    }

  }//the simtrack has no vertex
  else{useGen=true;}

  if (useGen){ 
    LogDebug(theCategory)<<"I am here 7";      
    //get corresponding gen particle
    int IndexGenPart = isimtk->genpartIndex();
    edm::LogVerbatim(theCategory)<<"The Gen Part index is: "<<IndexGenPart; 
    //do all the stuff to get the parent of mu from HEPMCproduct
    
    const HepMC::GenEvent *evt = hepmc->GetEvent();
    LogDebug(theCategory)<<"I am here 8";      
    
    edm::LogVerbatim(theCategory)<<"The gen particle index is: "<<IndexGenPart;      
    
    //skip it if not a valid GenIndex
    if(IndexGenPart<0){
      edm::LogError(theCategory)<<"The IndexGenPart is: "<<IndexGenPart
				<<"\n the sim momentum of this "<<selfID<<" is :"<<momentum.perp();return;}
    
    const HepMC::GenParticle *part = evt->barcode_to_particle(IndexGenPart);
    gentrack = part;
    if(!part){
      edm::LogError(theCategory)<<"My gen particle pointer is null";return;}
    
    int ipdg = part->pdg_id();
    edm::LogVerbatim(theCategory)<<"The Gen Part ID is: "<<ipdg; 
    
    if(ipdg==selfID)
      {
	HepMC::FourVector momentum_MC = part->momentum();
	edm::LogVerbatim(theCategory)<<"The Gen Part momentum is: "<<momentum_MC.perp(); 
	LogDebug(theCategory)<<"I am here 8.1";      
	
	if(!part->production_vertex()){
	  edm::LogError(theCategory)<<"there is no vertex to this Gen "<<selfID;
	  return;}
	
	const HepMC::GenVertex * gvertex = part->production_vertex();
	if(part->production_vertex()->particles_in_size()==0){
	  edm::LogError(theCategory)<<"there is no incoming partticle to this "<<selfID<<" Gen vertex.";
	  return;}
	
	const HepMC::GenParticle *mother = *(part->production_vertex()->particles_in_const_begin());
	
	LogDebug(theCategory)<<"I am here 8.2";		   
	while (abs(mother->pdg_id())==selfID)
	  {
	    if(!mother->production_vertex()){
	      edm::LogError(theCategory)<<"there is no vertex to this Gen "<<selfID<<". while looking recursively.";mother=0;break;}
	    gvertex = mother->production_vertex();
	    if(mother->production_vertex()->particles_in_size()==0){
	      edm::LogError(theCategory)<<"there is no incoming partticle to this muon Gen vertex. while looking recursively.";mother=0;break;}
	    mother = *(mother->production_vertex()->particles_in_const_begin());
	  }
	if (!mother){
	  edm::LogError(theCategory)<<"could not get a proper mother to this Gen "<<selfID; 
	  return;}
	
	Gen_vertex =gvertex;
	Gen_mother=mother;
	edm::LogVerbatim(theCategory)<<"the mother Id is: "<<mother->pdg_id();
	
      }// the gen associated to a tyep is of right type
    else{
      edm::LogError(theCategory)<<"the gen associated to a:"<< selfID<<" is: "<<ipdg;}
    LogDebug(theCategory)<<"I am here 9";      

  }
}

reco::Particle MotherSearch::particle(){
  /*
  if (SimIsValid()){
    reco::Particle::Point x(Sim_vertex->position().x(),
			    Sim_vertex->position().y(),
			    Sim_vertex->position().z());
    reco::Particle::LorentzVector p(simtrack->momentum().x(),
				    simtrack->momentum().y(),
				    simtrack->momentum().z());
    reco::Particle::Charge c(simtrack->charge());
    return reco::Particle(c, p, x, simtrack->type());
  }
    else if(GenIsValid()){
      reco::Particle::Point x(Gen_vertex->position().x(),
			      Gen_vertex->position().y(),
			      Gen_vertex->position().z());
      
      reco::Particle::LorentzVector p();
      //FIXME
      reco::Particle::Charge c(0);
      return reco::Particle(c, p, x, gentrack->pdg_id());
    }
    else*/{
      //default
      return reco::Particle();
    }
}

reco::Particle MotherSearch::mother(){
  /*
  if (SimIsValid()){
    reco::Particle::Point x(Sim_vertex->position().x(),
			    Sim_vertex->position().y(),
			    Sim_vertex->position().z());
    
    reco::Particle::LorentzVector p();
    reco::Particle::Charge c(Sim_mother->charge());
    return reco::Particle(c, p, x, Sim_mother->type());
  }
  else if(GenIsValid()){
    
    reco::Particle::Point x(Gen_vertex->position().x(),
			    Gen_vertex->position().y(),
			    Gen_vertex->position().z());

    reco::Particle::LorentzVector p();
    //FIXME
    reco::Particle::Charge c(0);
    return reco::Particle(c, p, x, Gen_mother->pdg_id());
  }
  else*/{
    //default
    return reco::Particle();
  }
}

