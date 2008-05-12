#ifndef NTupler_H
#define NTupler_H


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include <FWCore/Framework/interface/ProducerBase.h>

#include "FWCore/Framework/interface/Event.h"
#include "TTree.h"

/*
 * Description:
 * placeholder for common ntuplizer tools
 *
 */

//base generic class

class NTupler {
 public:
  NTupler() : useTFileService_(false){}
  virtual ~NTupler(){}

  virtual uint registerleaves(edm::ProducerBase * producer) =0;
  virtual void fill(edm::Event& iEvent)=0;
 protected:
  bool useTFileService_;
  TTree * tree_;
};


#endif
