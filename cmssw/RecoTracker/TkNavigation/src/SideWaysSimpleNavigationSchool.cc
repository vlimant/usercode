#include "RecoTracker/TkNavigation/interface/SideWaysSimpleNavigationSchool.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "RecoTracker/TkNavigation/interface/SimpleBarrelNavigableLayer.h"
#include "RecoTracker/TkNavigation/interface/SimpleForwardNavigableLayer.h"
#include "RecoTracker/TkNavigation/interface/SimpleNavigableLayer.h"
#include "RecoTracker/TkNavigation/interface/DiskLessInnerRadius.h"
#include "RecoTracker/TkNavigation/interface/SymmetricLayerFinder.h"

#include "TrackingTools/DetLayers/interface/BarrelDetLayer.h"
#include "TrackingTools/DetLayers/interface/ForwardDetLayer.h"
#include "TrackingTools/DetLayers/src/DetBelowZ.h"
#include "TrackingTools/DetLayers/src/DetLessZ.h"
#include "TrackingTools/DetLayers/interface/NavigationSetter.h"

#include "DataFormats/GeometrySurface/interface/BoundCylinder.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"

#include "Utilities/General/interface/CMSexception.h"

#include <functional>
#include <algorithm>
#include <map>
#include <cmath>

#ifndef COM
#define COM(var) std::cout<<"### "<<var<<std::endl;
//#define COM(var) ;
#endif

using namespace std;

SideWaysSimpleNavigationSchool::SideWaysSimpleNavigationSchool(const GeometricSearchTracker* theInputTracker,
					       const MagneticField* field) : 
  theBarrelLength(0),theField(field), theTracker(theInputTracker)
{
  // Get barrel layers
  /*sideways does not need barrels*/
  //but let's try it
  vector<BarrelDetLayer*> blc = theTracker->barrelLayers(); 
  for ( vector<BarrelDetLayer*>::iterator i = blc.begin(); i != blc.end(); i++) {
    theBarrelLayers.push_back( (*i) );
  }

  // get forward layers
  vector<ForwardDetLayer*> flc = theTracker->forwardLayers(); 
  for ( vector<ForwardDetLayer*>::iterator i = flc.begin(); i != flc.end(); i++) {
    theForwardLayers.push_back( (*i) );
  }
  
  FDLI middle = find_if( theForwardLayers.begin(), theForwardLayers.end(),
			 not1(DetBelowZ(0)));
  theLeftLayers  = FDLC( theForwardLayers.begin(), middle);
  theRightLayers = FDLC( middle, theForwardLayers.end());
  
  SymmetricLayerFinder symFinder( theForwardLayers);

  // only work on positive Z side; negative by mirror symmetry later
  /*sideways does not need barrels*/
  //but let's try it
  linkBarrelLayers( symFinder);

  linkForwardLayers( symFinder);
  COM("inverse relation");
  establishInverseRelations();

  //add the necessary inward links to end caps
  COM("linkOtherEndLayer");
  linkOtherEndLayers( symFinder);

}

SideWaysSimpleNavigationSchool::StateType 
SideWaysSimpleNavigationSchool::navigableLayers() const
{
  StateType result;
  for ( vector< SimpleBarrelNavigableLayer*>::const_iterator 
	  ib = theBarrelNLC.begin(); ib != theBarrelNLC.end(); ib++) {
    result.push_back( *ib);
  }
  for ( vector< SimpleForwardNavigableLayer*>::const_iterator 
	  ifl = theForwardNLC.begin(); ifl != theForwardNLC.end(); ifl++) {
    result.push_back( *ifl);
  }
  return result;
}

void SideWaysSimpleNavigationSchool::
linkBarrelLayers( SymmetricLayerFinder& symFinder)
{
  // Link barrel layers outwards
  for ( BDLI i = theBarrelLayers.begin(); i != theBarrelLayers.end(); i++) {
    BDLC reachableBL;
    FDLC leftFL;
    FDLC rightFL;

    // always add next barrel layer first
    if ( i+1 != theBarrelLayers.end()) reachableBL.push_back(*(i+1));
 
    // Add closest reachable forward layer (except for last BarrelLayer)
    if (i != theBarrelLayers.end() - 1) {
      linkNextForwardLayer( *i, rightFL);
    }

    // Add next BarrelLayer with length larger than the current BL
    if ( i+2 < theBarrelLayers.end()) {
      linkNextLargerLayer( i, theBarrelLayers.end(), reachableBL);
    }

    theBarrelNLC.push_back( new 
       SimpleBarrelNavigableLayer( *i, reachableBL,
				   symFinder.mirror(rightFL),
				   rightFL,theField, 5.));
  }
}

void SideWaysSimpleNavigationSchool::linkNextForwardLayer( BarrelDetLayer* bl, 
						   FDLC& rightFL)
{
  // find first forward layer with larger Z and larger outer radius
  float length = bl->surface().bounds().length() / 2.;
  float radius = bl->specificSurface().radius();
  for ( FDLI fli = theRightLayers.begin();
	fli != theRightLayers.end(); fli++) {
    if ( length < (**fli).position().z() &&
	 radius < (**fli).specificSurface().outerRadius()) {
      rightFL.push_back( *fli);
      return;
    }
  }
}

void SideWaysSimpleNavigationSchool::linkNextLargerLayer( BDLI bli, BDLI end,
						  BDLC& reachableBL)
{
  // compare length of next layer with length of following ones
  float length = (**(bli+1)).surface().bounds().length();
  float epsilon = 0.1;

  for ( BDLI i = bli+2; i < end; i++) {
    if ( length + epsilon < (**i).surface().bounds().length()) {
      reachableBL.push_back( *i);
      return;
    }
  }
}

void SideWaysSimpleNavigationSchool::
linkForwardLayers( SymmetricLayerFinder& symFinder)
{

  // handle right side first, groups are only on the right 
  vector<FDLC> groups = splitForwardLayers();

  LogDebug("TkNavigation") << "SideWaysSimpleNavigationSchool, Forward groups size = " << groups.size() ;
  for (vector<FDLC>::iterator g = groups.begin(); g != groups.end(); g++) {
    LogDebug("TkNavigation") << "group " << g - groups.begin() << " has " 
			     << g->size() << " layers " ;
  }

  for ( vector<FDLC>::iterator group = groups.begin();
	group != groups.end(); group++) {

    for ( FDLI i = group->begin(); i != group->end(); i++) {

      BDLC reachableBL;
      FDLC reachableFL;
 
      // Always connect to next barrel layer first, if exists
      linkNextBarrelLayer( *i, reachableBL);

      // Then always connect to next forward layer of "same" size, 
      // and layers of larger inner Radius
      linkNextLayerInGroup( i, *group, reachableFL);

      // Then connect to next N fw layers of next size
      if ( group+1 != groups.end()) {
	linkOuterGroup( *i, *(group+1), reachableFL);
      }

      // or connect within the group if outer radius increases
      linkWithinGroup( i, *group, reachableFL);

      theForwardNLC.push_back( new SimpleForwardNavigableLayer( *i,reachableBL,
								reachableFL,
								theField,
								5.));
      theForwardNLC.push_back( new SimpleForwardNavigableLayer( symFinder.mirror(*i),
								reachableBL,
								symFinder.mirror(reachableFL),
								theField,
								5.));

    }
  }

//    // now the left side by symmetry
//    for ( FDLI ileft = theLeftLayers.begin(); 
//  	ileft != theLeftLayers.end(); ileft++) {
//      ForwardDetLayer* right = symFinder.mirror( *ileft);
    
//      theForwardNLC.push_back( new 
//         SimpleForwardNavigableLayer( *ileft , right->nextBarrelLayers(),
//  	                      symFinder.mirror(right->nextForwardLayers())));
//    }
}

void SideWaysSimpleNavigationSchool::linkNextBarrelLayer( ForwardDetLayer* fl,
						  BDLC& reachableBL)
{
  if ( fl->position().z() > barrelLength()) return;

  float outerRadius = fl->specificSurface().outerRadius();
  float zpos        = fl->position().z();
  for ( BDLI bli = theBarrelLayers.begin(); bli != theBarrelLayers.end(); bli++) {
    if ( outerRadius < (**bli).specificSurface().radius() &&
	 zpos        < (**bli).surface().bounds().length() / 2.) {
      reachableBL.push_back( *bli);
      return;
    }
  }
}


void SideWaysSimpleNavigationSchool::linkNextLayerInGroup( FDLI fli,
						   const FDLC& group,
						   FDLC& reachableFL)
{
  // Always connect to next forward layer of "same" size, if exists
  if ( fli+1 != group.end()) {
    reachableFL.push_back( *(fli+1));
    // If that layer has an inner radius larger then the current one
    // also connect ALL next disks of same radius.
    float innerRThis = (**fli).specificSurface().innerRadius();
    float innerRNext =  (**(fli+1)).specificSurface().innerRadius();
    const float epsilon = 2.f;

    if (innerRNext > innerRThis + epsilon) {
      // next disk is smaller, so it doesn't cover fully subsequent ones
      // of same radius

      int i = 2;
      while ( (fli+i) != group.end()) {
	if ( (**(fli+i)).specificSurface().innerRadius() < 
	     innerRNext + epsilon) {
	  // following disk has not increased in ineer radius 
	  reachableFL.push_back( *(fli+i));
	  i++;
	} else {
	  break;
	}
      }
    }
  }
}


void SideWaysSimpleNavigationSchool::linkOuterGroup( ForwardDetLayer* fl,
					     const FDLC& group,
					     FDLC& reachableFL)
{

  // insert N layers with Z grater than fl

  ConstFDLI first = find_if( group.begin(), group.end(), 
			     not1( DetBelowZ( fl->position().z())));
  if ( first != group.end()) {

    // Hard-wired constant!!!!!!
    ConstFDLI last = min( first + 7, group.end());

    reachableFL.insert( reachableFL.end(), first, last);
  }
}

void SideWaysSimpleNavigationSchool::linkWithinGroup( FDLI fl,
					      const FDLC& group,
					      FDLC& reachableFL)
{
  ConstFDLI biggerLayer = outerRadiusIncrease( fl, group);
  if ( biggerLayer != group.end() && biggerLayer != fl+1) {
    reachableFL.push_back( *biggerLayer);
  }
}

SideWaysSimpleNavigationSchool::ConstFDLI
SideWaysSimpleNavigationSchool::outerRadiusIncrease( FDLI fl, const FDLC& group)
{
  const float epsilon = 5.f;
  float outerRadius = (**fl).specificSurface().outerRadius();
  while ( ++fl != group.end()) {
    if ( (**fl).specificSurface().outerRadius() > outerRadius + epsilon) {
      return fl;
    }
  }
  return fl;
}

vector<SideWaysSimpleNavigationSchool::FDLC> 
SideWaysSimpleNavigationSchool::splitForwardLayers() 
{
  // only work on positive Z side; negative by mirror symmetry later

  FDLC myRightLayers( theRightLayers);
  FDLI begin = myRightLayers.begin();
  FDLI end   = myRightLayers.end();

  // sort according to inner radius
  sort ( begin, end, DiskLessInnerRadius()); 

  // partition in cylinders
  vector<FDLC> result;
  FDLC current;
  current.push_back( *begin);
  for ( FDLI i = begin+1; i != end; i++) {

    LogDebug("TkNavigation") << "(**i).specificSurface().innerRadius()      = "
			     << (**i).specificSurface().innerRadius() << endl
			     << "(**(i-1)).specificSurface().outerRadius()) = "
			     << (**(i-1)).specificSurface().outerRadius() ;

    // if inner radius of i is larger than outer radius of i-1 then split!
    if ( (**i).specificSurface().innerRadius() > 
	 (**(i-1)).specificSurface().outerRadius()) {

      LogDebug("TkNavigation") << "found break between groups" ;

      // sort layers in group along Z
      sort ( current.begin(), current.end(), DetLessZ());

      result.push_back(current);
      current.clear();
    }
    current.push_back(*i);
  }
  result.push_back(current); // save last one too 

  // now sort subsets in Z
  for ( vector<FDLC>::iterator ivec = result.begin();
	ivec != result.end(); ivec++) {
    sort( ivec->begin(), ivec->end(), DetLessZ());
  }

  return result;
}

float SideWaysSimpleNavigationSchool::barrelLength() 
{
  if ( theBarrelLength < 1.) {
    for (BDLI i=theBarrelLayers.begin(); i!=theBarrelLayers.end(); i++) {
      theBarrelLength = max( theBarrelLength,
			     (**i).surface().bounds().length() / 2.f);
    }

    LogDebug("TkNavigation") << "The barrel length is " << theBarrelLength ;
  }
  return theBarrelLength;
}



void SideWaysSimpleNavigationSchool::establishInverseRelations() {
  COM("setter");
  NavigationSetter setter(*this);

  // find for each layer which are the barrel and forward
  // layers that point to it
  typedef map<const DetLayer*, vector<BarrelDetLayer*>, less<const DetLayer*> > BarrelMapType;
  typedef map<const DetLayer*, vector<ForwardDetLayer*>, less<const DetLayer*> > ForwardMapType;
  
  
  BarrelMapType reachedBarrelLayersMap;
  ForwardMapType reachedForwardLayersMap;
  
  COM("do barrels");
  for ( BDLI bli = theBarrelLayers.begin();
        bli!=theBarrelLayers.end(); bli++) {
    DLC reachedLC = (**bli).nextLayers( insideOut);
    for ( DLI i = reachedLC.begin(); i != reachedLC.end(); i++) {
      reachedBarrelLayersMap[*i].push_back( *bli);
    }
  }

  COM("do disks");
  for ( FDLI fli = theForwardLayers.begin();
        fli!=theForwardLayers.end(); fli++) {
    DLC reachedLC = (**fli).nextLayers( insideOut);
    for ( DLI i = reachedLC.begin(); i != reachedLC.end(); i++) {
      reachedForwardLayersMap[*i].push_back( *fli);
    }
    }
  
  
  COM("set inward links");
  vector<DetLayer*> lc = theTracker->allLayers();
  for ( vector<DetLayer*>::iterator i = lc.begin(); i != lc.end(); i++) {
    SimpleNavigableLayer* navigableLayer =
      dynamic_cast<SimpleNavigableLayer*>((**i).navigableLayer());
    if (!navigableLayer) COM("you'r in trouble man ! but hey no worry. that's normal");
    if (navigableLayer)
      navigableLayer->setInwardLinks( reachedBarrelLayersMap[*i],reachedForwardLayersMap[*i] );
  }
  
}


void SideWaysSimpleNavigationSchool::
linkOtherEndLayers(  SymmetricLayerFinder& symFinder){
  COM("reachable from horizontal");
  //generally, on the right side, what are the forward layers reachable from the horizontal
  FDLC reachableFL= reachableFromHorizontal();

  //even simpler navigation from end to end.
  //for each of them
  for (FDLI fl=reachableFL.begin();fl!=reachableFL.end();fl++)
    {
      COM("adding inward from right");
      //link it inward to the mirror reachable from horizontal
      addInward((DetLayer*)*fl,symFinder.mirror(*fl));
      
      COM("adding inward from mirror of right (left?)");
      addInward((DetLayer*)symFinder.mirror(*fl),*fl);
    }

  /* this is not enough to set reachable from each of them: too few links
     //this is enough in the end
     //for each of them
     for (FDLI fl=reachableFL.begin();fl!=reachableFL.end();fl++)
     {
     COM("adding inward from right");
     //link it inward to the mirror reachable from horizontal
     addInward((DetLayer*)*fl,symFinder.mirror(reachableFL));
     
     COM("adding inward from mirror of right (left?)");
     //do the same from the the mirrored layer to the reachable from horizontal
     addInward((DetLayer*)symFinder.mirror(*fl),reachableFL);
     }*/




  /* what about linking every not masked layer in each group.
     except for within the same group
     
     vector<FDLC> groups splitForwardLayers();
     FDLC reachable;
     
     for ( vector<FDLC>::iterator group = groups.begin();
     group != groups.end(); group++) {
     //for each group
     
     for ( FDLI i = group->begin(); i != group->end(); i++) {
     
     for ( vector<FDLC>::iterator other_group = groups.begin();
     other_group != groups.end(); other_group++) {
     //for each other group
     
     if (other_group==group && i==group->begin()){
     //other_group is the same as group and dealing with the first layer of the group
     //link the first of each group
     reachable.push_back(other_group.front());
     continue;}
     
     //now dealing as if other_group is different than group
     for ( FDLI other_i = other_group->begin(); other_i != other_group->end(); other_i++) {
     //for each of other group
     //is the layer in the other group "masking" this one
     //inner radius smaller OR outer radius bigger
     if ((**other_i).specificSurface().innerRadius() < (**i).specificSurface().innerRadius() ||
     (**other_i).specificSurface().outerRadius() > (**i).specificSurface().outerRadius())
     {         //not masked
     reachableFL.push_back(*other_i);
     }


     }

	

     //do something special with the first of each group
     //do somethign special with layers in its own group

     }
      
     }


     }

  */



   /* this is too much links between layers
  FDLC myRightLayers( theRightLayers);
  FDLI begin = myRightLayers.begin();
  FDLI end   = myRightLayers.end();
  
  //  for each of the right layers
  for (FDLI fl = begin;fl!=end;++fl)
    {
      //get the navigable layer for this DetLayer
      SimpleNavigableLayer* navigableLayer =
	dynamic_cast<SimpleNavigableLayer*>((*fl)->navigableLayer());
      
      COM("retreive the next layers outsidein");
      //get the OUTward reachable layers.
      DLC inwardsLayers(navigableLayer->nextLayers(insideOut));

      //what is reachable horizontaly
      FDLC thisReachableFL (reachableFL);

      COM("combine");
      //combine the two vector with a conversion to forward layer
      for (DLI i=inwardsLayers.begin();i!=inwardsLayers.end();++i)
	{
	  ForwardDetLayer* fd=dynamic_cast<ForwardDetLayer*>(const_cast<DetLayer*>(*i));
	  //	  ForwardDetLayer* fd=const_cast<ForwardDetLayer*>(*i);
	  if (fd){
	    //	    if (thisReachableFL.find(fd)==thisReachableFL.end())
	    //	      {//no duplicate. insert it
		thisReachableFL.push_back(fd);
		//}
	  }
	  else COM("investigate");
	}
      
      //please no duplicate !!!
      COM("no duplicate");
      FDLI new_end =unique(thisReachableFL.begin(),thisReachableFL.end());
      thisReachableFL.erase(new_end,thisReachableFL.end());

      //then set the inwards links
      COM("adding inward from right");
      //link it inward to the mirror reachable from horizontal
      addInward((DetLayer*)*fl,symFinder.mirror(thisReachableFL));
      
      COM("adding inward from mirror of right (left?)");
      //do the same from the the mirrored layer to the reachable from horizontal
      addInward((DetLayer*)symFinder.mirror(*fl),thisReachableFL);
    }
   */


}

void SideWaysSimpleNavigationSchool::
addInward(DetLayer * det, ForwardDetLayer * newF){
  //get the navigable layer for this DetLayer
  SimpleNavigableLayer* navigableLayer =
    dynamic_cast<SimpleNavigableLayer*>((*det).navigableLayer());

  COM("retreive the nextlayer outsidein");
  //get the inward reachable layers.
  DLC inwardsLayers(navigableLayer->nextLayers(outsideIn));
  
  COM("split them barrel/forward");
  // split barrel and forward layers
  BDLC inwardsBarrel;
  FDLC inwardsForward;
  for ( DLC::iterator dli=inwardsLayers.begin();dli!=inwardsLayers.end();dli++)
    {
      if ((**dli).location()==GeomDetEnumerators::barrel)
        inwardsBarrel.push_back((BarrelDetLayer*)*dli);
      else
        inwardsForward.push_back((ForwardDetLayer*)*dli);
    }
  COM("add the new ones");
  //add the other forward layers provided
  inwardsForward.push_back(newF);

  COM("no duplicate please");
  FDLI new_end =unique(inwardsForward.begin(),inwardsForward.end());
  inwardsForward.erase(new_end,inwardsForward.end());

  COM("set back the inward links (no duplicate)");
  //  set them back to the navigable layer
  navigableLayer->setInwardLinks( inwardsBarrel, inwardsForward);
}

void SideWaysSimpleNavigationSchool::
addInward(DetLayer * det, FDLC news){
  //get the navigable layer for this DetLayer
  SimpleNavigableLayer* navigableLayer =
    dynamic_cast<SimpleNavigableLayer*>((*det).navigableLayer());

  COM("retreive the nextlayer outsidein");
  //get the inward reachable layers.
  DLC inwardsLayers(navigableLayer->nextLayers(outsideIn));

  COM("split them barrel/forward");
  // split barrel and forward layers
  BDLC inwardsBarrel;
  FDLC inwardsForward;
  for ( DLC::iterator dli=inwardsLayers.begin();dli!=inwardsLayers.end();dli++)
    {
      if ((**dli).location()==GeomDetEnumerators::barrel)
	inwardsBarrel.push_back((BarrelDetLayer*)*dli);	
      else
	inwardsForward.push_back((ForwardDetLayer*)*dli);
    }
  
  COM("add the new ones");
  //add the other forward layers provided
  inwardsForward.insert( inwardsForward.end(), news.begin(), news.end());

  COM("no duplicate please");
  FDLI new_end =unique(inwardsForward.begin(),inwardsForward.end());
  inwardsForward.erase(new_end,inwardsForward.end());

  COM("set back the inward links (no duplicate)");
  //  set them back to the navigable layer
  navigableLayer->setInwardLinks( inwardsBarrel, inwardsForward);
}

SideWaysSimpleNavigationSchool::FDLC
SideWaysSimpleNavigationSchool::reachableFromHorizontal()
{    
  //determine which is the list of forward layers that can be reached from inside-out
  //at horizontal direction

  FDLC myRightLayers( theRightLayers);
  FDLI begin = myRightLayers.begin();
  FDLI end   = myRightLayers.end();

  //sort along Z to be sure
  sort(begin, end, DetLessZ());

  FDLC reachableFL;

  begin = myRightLayers.begin();
  end   = myRightLayers.end();

  //the first one is always reachable
  reachableFL.push_back(*begin);
  FDLI current = begin;
  for (FDLI i = begin+1; i!= end; i++)
    {
      //is the previous layer NOT masking this one
      //inner radius smaller OR outer radius bigger
      if ((**i).specificSurface().innerRadius() < (**current).specificSurface().innerRadius() ||
	  (**i).specificSurface().outerRadius() > (**current).specificSurface().outerRadius())
	{	  //not masked
	  reachableFL.push_back(*i);
	  current=i;
	}
    }
  return reachableFL;
}
