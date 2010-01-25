// -*- C++ -*-
//
// Package:   EcalTimePi0Tree
// Class:     EcalTimePi0Tree
//
/**\class EcalTimePi0Tree EcalTimePi0Tree.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
 */
//
// Skeleton Derived from an example by:  F. DE GUIO C. DOGLIONI P. MERIDIANI
// Authors:                              Seth Cooper, Giovanni Franzoni (UMN)
//         Created:  Mo Jul 14 5:46:22 CEST 2008
// $Id: EcalTimePi0Tree.cc,v 1.1 2010/01/25 09:45:55 franzoni Exp $
//
//

 
#include "ECALTime/EcalTimePi0/plugins/EcalTimePi0Tree.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

//#include "CaloOnlineTools/EcalTools/interface/calcEcalDeposits.h"

#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h"

#include "CondFormats/EcalObjects/interface/EcalADCToGeVConstant.h"
#include "CondFormats/DataRecord/interface/EcalADCToGeVConstantRcd.h"


#include <vector>
#include <iostream>
#include <iterator>

using namespace cms ;
using namespace edm ;
using namespace std ;



EcalTimePi0Tree::EcalTimePi0Tree (const edm::ParameterSet& iConfig) :
  barrelEcalRecHitCollection_              (iConfig.getParameter<edm::InputTag> ("barrelEcalRecHitCollection")),
  endcapEcalRecHitCollection_              (iConfig.getParameter<edm::InputTag> ("endcapEcalRecHitCollection")),
  barrelBasicClusterCollection_            (iConfig.getParameter<edm::InputTag> ("barrelBasicClusterCollection")),
  endcapBasicClusterCollection_            (iConfig.getParameter<edm::InputTag> ("endcapBasicClusterCollection")),
  barrelSuperClusterCollection_            (iConfig.getParameter<edm::InputTag> ("barrelSuperClusterCollection")),
  endcapSuperClusterCollection_            (iConfig.getParameter<edm::InputTag> ("endcapSuperClusterCollection")),
  muonCollection_                          (iConfig.getParameter<edm::InputTag> ("muonCollection")),
  l1GMTReadoutRecTag_   (iConfig.getUntrackedParameter<std::string> ("L1GlobalReadoutRecord","gtDigis")),
  runNum_               (iConfig.getUntrackedParameter<int> ("runNum")),
  fileName_             (iConfig.getUntrackedParameter<std::string> ("fileName", std::string ("EcalTimePi0Tree"))),
  naiveId_ (0)              

{
  // TrackAssociator parameters // gfwork: can we remove this? 
  edm::ParameterSet trkParameters = iConfig.getParameter<edm::ParameterSet> ("TrackAssociatorParameters") ;
  trackParameters_.loadParameters ( trkParameters ) ;
  //trackAssociator_.useDefaultPropagator () ;

  // Create File
  fileName_ += "_"+intToString (runNum_)+".root" ;
  file_ = new TFile (fileName_.c_str () , "RECREATE") ;
  file_->cd () ;

  // Initialize Tree
  tree_ = new TTree ( "EcalTimePi0Analysis","EcalTimePi0Analysis" ) ;

  // GFdoc EcalCosmicsTreeContent builds branches for group of variables (trg, muon, ecal, hcal, tkass, tpg, l1)
  // GFdoc after setting branches, setBranchAddresses will be called 

  setBranches (tree_, myTreeVariables_) ;

}


// -----------------------------------------------------------------------------------------


EcalTimePi0Tree::~EcalTimePi0Tree ()
{
}


// -----------------------------------------------------------------------------------------


void EcalTimePi0Tree::analyze (const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  ++naiveId_ ;
  //std::cout << ("EcalTimePi0Tree") << "event " << iEvent.id ().event () << "   "
  //	    << "naiveEvent " <<naiveId_ << "\n" << std::endl;  
  

  // Geometry
  edm::ESHandle<CaloGeometry> pGeometry ;
  iSetup.get<CaloGeometryRecord> ().get (pGeometry) ;
  const CaloGeometry * theGeometry = pGeometry.product () ;
  // Topology
  edm::ESHandle<CaloTopology> pCaloTopology ;
  iSetup.get<CaloTopologyRecord> ().get (pCaloTopology) ;
  const CaloTopology * theCaloTopology = pCaloTopology.product () ;
  
  

  // Ecal barrel RecHits 
  edm::Handle<EcalRecHitCollection> pBarrelEcalRecHits ;
  iEvent.getByLabel (barrelEcalRecHitCollection_, pBarrelEcalRecHits) ;
  const EcalRecHitCollection* theBarrelEcalRecHits = pBarrelEcalRecHits.product () ;   

  if (! (pBarrelEcalRecHits.isValid ()) )
    {
      LogWarning ("EcalTimePi0Tree") << barrelEcalRecHitCollection_ 
                                     << " not available" ;
      return ;
    }

  // Ecal endcap RecHits
  edm::Handle<EcalRecHitCollection> pEndcapEcalRecHits ;
  iEvent.getByLabel (endcapEcalRecHitCollection_, pEndcapEcalRecHits) ;
  const EcalRecHitCollection* theEndcapEcalRecHits = pEndcapEcalRecHits.product () ;   
  
  if (! (pEndcapEcalRecHits.isValid ()))
    {
      LogWarning ("EcalTimePi0Tree") << endcapEcalRecHitCollection_ 
                                     << " not available" ;
      return ;
    }
  
  
  
  
  // GFdo switch to appropriate clusters here (basic instead of super) 
  // Barrel SuperClusters
  edm::Handle<reco::SuperClusterCollection> pBarrelSuperClusters ;
  iEvent.getByLabel (barrelSuperClusterCollection_, pBarrelSuperClusters) ;
  const reco::SuperClusterCollection* theBarrelSuperClusters = pBarrelSuperClusters.product () ;

  if (! (pBarrelSuperClusters.isValid ()) )
    {
      LogWarning ("EcalTimePi0Tree") << barrelSuperClusterCollection_ 
                                     << " not available" ;
      return ;
    }

  // Endcap SuperClusters
  edm::Handle<reco::SuperClusterCollection> pEndcapSuperClusters ;
  iEvent.getByLabel (endcapSuperClusterCollection_, pEndcapSuperClusters) ;
  const reco::SuperClusterCollection* theEndcapSuperClusters = pEndcapSuperClusters.product () ;
  
  if (! (pEndcapSuperClusters.isValid ()) )
    {
      LogWarning ("EcalTimePi0Tree") << endcapSuperClusterCollection_ 
                                     << " not available" ;
      return ;
    }
  
  
  // Barrel BasicClusters
  edm::Handle<reco::BasicClusterCollection> pBarrelBasicClusters ;
  iEvent.getByLabel (barrelBasicClusterCollection_, pBarrelBasicClusters) ;
  const reco::BasicClusterCollection* theBarrelBasicClusters = pBarrelBasicClusters.product () ;

  if (! (pBarrelBasicClusters.isValid ()) )
    {
      LogWarning ("EcalTimePi0Tree") << barrelBasicClusterCollection_ 
                                     << " not available" ;
      return ;
    }

  // Endcap BasicClusters
  edm::Handle<reco::BasicClusterCollection> pEndcapBasicClusters ;
  iEvent.getByLabel (endcapBasicClusterCollection_, pEndcapBasicClusters) ;
  const reco::BasicClusterCollection* theEndcapBasicClusters = pEndcapBasicClusters.product () ;
  
  if (! (pEndcapBasicClusters.isValid ()) )
    {
      LogWarning ("EcalTimePi0Tree") << endcapBasicClusterCollection_ 
                                     << " not available" ;
      return ;
    }
  
  
  
  
  
  // ClusterShapes
  EcalClusterLazyTools* lazyTools = new EcalClusterLazyTools(iEvent, iSetup, barrelEcalRecHitCollection_, endcapEcalRecHitCollection_);


  // Xtal - TkLength map
  std::map<int,float> XtalMap ;
  std::map<int,float> XtalMapCurved ;
  std::map<int,float> XtalMapCurved_high ;
  std::map<int,float> XtalMapCurved_low ;

  // GFdoc initialize variables to 0/false
  initializeBranches(tree_, myTreeVariables_);

  myTreeVariables_.runId = iEvent.id ().run () ;
  myTreeVariables_.eventId = iEvent.id ().event () ;
  myTreeVariables_.eventNaiveId = naiveId_ ;
  myTreeVariables_.timeStampLow = ( 0xFFFFFFFF & iEvent.time ().value () ) ;
  myTreeVariables_.timeStampHigh = ( iEvent.time ().value () >> 32 ) ;

  // GFdoc GT information 
  dump3Ginfo(iEvent, iSetup, myTreeVariables_) ;

  dumpBarrelClusterInfo(theGeometry, theCaloTopology,
                   theBarrelEcalRecHits, theBarrelBasicClusters, theBarrelSuperClusters, lazyTools, XtalMap, XtalMapCurved, myTreeVariables_) ;
  dumpEndcapClusterInfo(theGeometry, theCaloTopology,
                   theEndcapEcalRecHits, theEndcapBasicClusters, theEndcapSuperClusters, lazyTools, XtalMap, XtalMapCurved, myTreeVariables_) ;

  /*
  dumpHcalInfo(theGeometry, iEvent, iSetup, myTreeVariables_);
  */

  tree_ -> Fill();
}






// -----------------------------------------------------------------------------------------

void EcalTimePi0Tree::endJob ()
{
  using namespace std ;
  file_->cd () ;
  tree_->Write () ;
  file_->Close () ;

}






// -----------------------------------------------------------------------------------------

std::string EcalTimePi0Tree::intToString (int num)
{
    using namespace std ;
    ostringstream myStream ;
    myStream << num << flush ;
    return (myStream.str ()) ; //returns the string form of the stringstream object
}






// -----------------------------------------------------------------------------------------

void EcalTimePi0Tree::dumpBarrelClusterInfo (const CaloGeometry * theGeometry,
                                        const CaloTopology * theCaloTopology,
                                        const EcalRecHitCollection* theBarrelEcalRecHits,
                                        const reco::BasicClusterCollection* theBarrelBasicClusters,
                                        const reco::SuperClusterCollection* theBarrelSuperClusters,
                                        EcalClusterLazyTools* lazyTools,
                                        const std::map<int,float> & XtalMap, //GFdoc unclear
                                        const std::map<int,float> & XtalMapCurved, //GFdoc unclear
                                        EcalTimePi0TreeContent & myTreeVariables_)
{            
  // get number of of objects already present in the tree
  int numberOfSuperClusters = myTreeVariables_.nSuperClusters;
  int numberOfClusters = myTreeVariables_.nClusters ;
  int numberOfXtals = myTreeVariables_.nXtals ;
  
  
  //number of superClusters in event (collection = vector!)
  myTreeVariables_.nSuperClusters += theBarrelSuperClusters->size () ;
  myTreeVariables_.nBarrelSuperClusters = theBarrelSuperClusters->size () ;
  
  /////////////////////////////////////////////////////////////////////////////////////////
  //loop on all superclusters in event
  for (reco::SuperClusterCollection::const_iterator sclus = theBarrelSuperClusters->begin () ; 
       sclus != theBarrelSuperClusters->end () ; 
       ++sclus) 
    {
  
      int numberOfXtalsInSuperCluster = 0 ;//counter for all xtals in supercluster 
      
      myTreeVariables_.nClustersInSuperCluster[numberOfSuperClusters] = sclus -> clustersSize () ;
      
      myTreeVariables_.superClusterType[numberOfSuperClusters] = 0 ;
      
      myTreeVariables_.superClusterEta[numberOfSuperClusters] = sclus -> position ().eta () ;
      myTreeVariables_.superClusterPhi[numberOfSuperClusters] = sclus -> position ().phi () ;
      myTreeVariables_.superClusterX[numberOfSuperClusters] = sclus -> position ().x () ;
      myTreeVariables_.superClusterY[numberOfSuperClusters] = sclus -> position ().y () ;
      myTreeVariables_.superClusterZ[numberOfSuperClusters] = sclus -> position ().z () ;
      myTreeVariables_.superClusterRawEnergy[numberOfSuperClusters] = sclus -> rawEnergy () ;
      myTreeVariables_.superClusterPhiWidth[numberOfSuperClusters] = sclus -> phiWidth () ;
      myTreeVariables_.superClusterEtaWidth[numberOfSuperClusters] = sclus -> etaWidth () ;
      myTreeVariables_.clusterIndexInSuperCluster[numberOfSuperClusters] = numberOfClusters ;
      myTreeVariables_.xtalIndexInSuperCluster[numberOfSuperClusters] = numberOfXtals ;
      //float & energySum = myTreeVariables_.superClusterEnergySum[numberOfSuperClusters];
      //energySum = 0.;


      ///////////////////////////////////////////////////////////////////////////////////////
      // independent loop on barrel basic clusters
      float & energySum = myTreeVariables_.superClusterEnergySum[numberOfSuperClusters];
      energySum = 0.;
      for (reco::BasicClusterCollection::const_iterator clus = theBarrelBasicClusters->begin () ; 
	   clus != theBarrelBasicClusters->end () ;  
	   ++clus) // loop on barrel Bclusters
	{        
	  double energy = (clus)->energy () ;
	  double phi    = (clus)->phi () ;
	  double eta    = (clus)->eta () ;
	  double sinTheta         = fabs( sin( 2 *atan( exp(-1*(clus)->eta()) ) ) );
	  double transverseEnergy = (clus)->energy () * sinTheta;
	  double time = -1000.0 ; // gfdoc: work on this to provide a combination of crystals?
	  double ampli = 0. ;
	  double secondMin = 0. ;
	  double secondTime = -1000. ;
	  int numberOfXtalsInCluster=0 ;//counter for all xtals in cluster
	  int numXtalsinCluster = 0 ;//xtals in cluster above 3sigma noise  
	  EBDetId maxDet ;
	  EBDetId secDet ;
	  
	  // GFdoc clusterDetIds holds crystals that participate to this basic cluster 
	  std::vector<std::pair<DetId, float> > clusterDetIds = (clus)->hitsAndFractions() ; //get these from the cluster
	  //std::cout << "gf size of eb Bcluster: " << clusterDetIds.size() << std::endl;	  

	  myTreeVariables_.xtalIndexInCluster[numberOfClusters] = numberOfXtals ;
	  
	  //////////////////////////////////////////////////////////////////////////////////////
	  //loop on xtals in cluster
	  for (std::vector<std::pair<DetId, float> >::const_iterator detitr = clusterDetIds.begin () ; 
	       detitr != clusterDetIds.end () ; 
	       ++detitr)// loop on rechics of barrel basic clusters
	    {
	      //Here I use the "find" on a digi collection... I have been warned...   (GFdoc: ??)
	      // GFdoc: check if DetId belongs to ECAL; if so, find it among those if this basic cluster
	      if ( (detitr -> first).det () != DetId::Ecal) 
		{ 
		  std::cout << " det is " << (detitr -> first).det () << " (and not DetId::Ecal)" << std::endl ;
		  continue ;
		}
	      if ( (detitr -> first).subdetId () != EcalBarrel) 
		{
		  std::cout << " subdet is " << (detitr -> first).subdetId () << " (and not EcalBarrel)" << std::endl ; 
		  continue ; 
		}
	      
	      // GFdoc now find it!
	      EcalRecHitCollection::const_iterator thishit = theBarrelEcalRecHits->find ( (detitr -> first) ) ;
	      if (thishit == theBarrelEcalRecHits->end ()) 
		{
		  continue ;
		}
	      //The checking above should no longer be needed...... 
	      //as only those in the cluster would already have rechits..
	      
	      // GFdoc this is one crystal in the basic cluster
	      EcalRecHit myhit = (*thishit) ;
	      
	      // GFdoc do we want this cut here???/////////////////////////////////////////////////
	      double thisamp = myhit.energy () ;
	      if (thisamp > 0.027) //cut on energy->number of crystals in cluster above 3sigma noise
		{ 
		  numXtalsinCluster++ ; 
		}
	      if (thisamp > secondMin) 
		{
		  secondMin = thisamp ; 
		  secondTime = myhit.time () ; 
		  secDet = (EBDetId) (detitr -> first) ;
		}
	      if (secondMin > ampli) 
		{
		  std::swap (ampli, secondMin) ; 
		  std::swap (time, secondTime) ; 
		  std::swap (maxDet, secDet) ;
		}
    
	      myTreeVariables_.xtalEnergy[numberOfXtals] = (float) thisamp ; //gf: add uncalRecHits here?
	      myTreeVariables_.xtalTime[numberOfXtals] = (float) (myhit.time ()) ;
	      myTreeVariables_.xtalHashedIndex[numberOfXtals] = EBDetId (detitr -> first).hashedIndex () ;

	      // xtal variables inside a barrel basic cluster 
	      myTreeVariables_.xtalInBCEnergy[numberOfClusters][numberOfXtalsInCluster]=      (float) thisamp;
	      myTreeVariables_.xtalInBCTime[numberOfClusters][numberOfXtalsInCluster]=        (float) (myhit.time ()); 
	      myTreeVariables_.xtalInBCHashedIndex[numberOfClusters][numberOfXtalsInCluster]= EBDetId (detitr -> first).hashedIndex () ; 
	      myTreeVariables_.xtalInBCIeta[numberOfClusters][numberOfXtalsInCluster]=        EBDetId((detitr -> first)).ieta(); 
	      myTreeVariables_.xtalInBCIphi[numberOfClusters][numberOfXtalsInCluster]=        EBDetId((detitr -> first)).iphi();
	      myTreeVariables_.xtalInBCIx[numberOfClusters][numberOfXtalsInCluster]=          0; 
	      myTreeVariables_.xtalInBCIy[numberOfClusters][numberOfXtalsInCluster]=          0; 
	      myTreeVariables_.xtalInBCFlag[numberOfClusters][numberOfXtalsInCluster]=        myhit.recoFlag(); 

	      energySum += (float) thisamp ; // GFdoc incrementing energy of SC
	      
	      //MF Lenght evaluation in XTals
	      int raw = (detitr -> first).rawId () ;
	      
	      // GFdoc: unclear
	      if (XtalMap.find (raw) != XtalMap.end ())
		myTreeVariables_.xtalTkLength[numberOfXtals] = XtalMap.find (raw)->second ;
	      else
		myTreeVariables_.xtalTkLength[numberOfXtals] = -1. ;
              
	      // GFdoc: unclear
	      if (XtalMapCurved.find (raw) != XtalMapCurved.end ())
		myTreeVariables_.xtalTkLengthCurved[numberOfXtals] = XtalMapCurved.find (raw)->second ;
	      else
		myTreeVariables_.xtalTkLengthCurved[numberOfXtals] = -1. ;
              
	      numberOfXtals++ ; // GFdoc increment number of crystals in basic cluster 
	      numberOfXtalsInCluster++ ;
	      numberOfXtalsInSuperCluster++ ;
	      
	    } //end loop on rechics of barrel basic clusters
	  
	  float E2 = (float)ampli + (float)secondMin ;
	  
         //Set some more values for this cluster
         myTreeVariables_.clusterEnergy[numberOfClusters] = (float) energy ;
         myTreeVariables_.clusterTransverseEnergy[numberOfClusters] = (float) transverseEnergy ;
         myTreeVariables_.clusterTime[numberOfClusters] = (float) time ;
         myTreeVariables_.clusterEta[numberOfClusters] = (float) eta ;
         myTreeVariables_.clusterPhi[numberOfClusters] = (float) phi ;
         myTreeVariables_.clusterE1[numberOfClusters] = (float) ampli ;
         myTreeVariables_.clusterE2[numberOfClusters] = (float) E2 ;
         myTreeVariables_.clusterXtalsAbove3Sigma[numberOfClusters] = numXtalsinCluster ;
         myTreeVariables_.clusterXtals[numberOfClusters] = clusterDetIds.size () ;
         myTreeVariables_.clusterMaxId[numberOfClusters] =  maxDet.rawId () ;
         myTreeVariables_.clusterMaxId[numberOfClusters] =  secDet.rawId () ;
         myTreeVariables_.nXtalsInCluster[numberOfClusters]= numberOfXtalsInCluster ;    


	 // (basic) cluster shapes for barrel
         myTreeVariables_.clusterE2x2[numberOfClusters] = lazyTools -> e2x2(*(clus));
         myTreeVariables_.clusterE3x2[numberOfClusters] = lazyTools -> e3x2(*(clus));
         myTreeVariables_.clusterE3x3[numberOfClusters] = lazyTools -> e3x3(*(clus));
         myTreeVariables_.clusterE4x4[numberOfClusters] = lazyTools -> e4x4(*(clus));
         myTreeVariables_.clusterE5x5[numberOfClusters] = lazyTools -> e5x5(*(clus));

         myTreeVariables_.clusterE2x5Right[numberOfClusters] = lazyTools -> e2x5Right(*(clus));
         myTreeVariables_.clusterE2x5Left[numberOfClusters] = lazyTools -> e2x5Left(*(clus));
         myTreeVariables_.clusterE2x5Top[numberOfClusters] = lazyTools -> e2x5Top(*(clus));
         myTreeVariables_.clusterE2x5Bottom[numberOfClusters] = lazyTools -> e2x5Bottom(*(clus));
         //myTreeVariables_.clusterE3x2Ratio[numberOfClusters] = lazyTools -> e3x2Ratio(*(clus));

         myTreeVariables_.clusterCovEtaEta[numberOfClusters] = (lazyTools -> covariances(*(clus)))[0];
         myTreeVariables_.clusterCovPhiPhi[numberOfClusters] = (lazyTools -> covariances(*(clus)))[2];
         myTreeVariables_.clusterCovEtaPhi[numberOfClusters] = (lazyTools -> covariances(*(clus)))[1];

         myTreeVariables_.clusterLat[numberOfClusters] = (lazyTools -> lat(*(clus)))[2];
         myTreeVariables_.clusterEtaLat[numberOfClusters] = (lazyTools -> lat(*(clus)))[0];
         myTreeVariables_.clusterPhiLat[numberOfClusters] = (lazyTools -> lat(*(clus)))[1];

         myTreeVariables_.clusterZernike20[numberOfClusters] = lazyTools -> zernike20(*(clus));
         myTreeVariables_.clusterZernike42[numberOfClusters] = lazyTools -> zernike42(*(clus));


    
         numberOfClusters++ ;

       } //end (BASIC)cluster loop

      myTreeVariables_.nXtalsInSuperCluster[numberOfSuperClusters]=numberOfXtalsInSuperCluster ;
      numberOfSuperClusters++ ;
  
    } //end supercluster loop

  //number of clusters / xtals for the whole event
  myTreeVariables_.nClusters += numberOfClusters ;
  myTreeVariables_.nXtals += numberOfXtals ;

  return ;
} // dumpBarrelClusterInfo  



void EcalTimePi0Tree::dumpEndcapClusterInfo (const CaloGeometry * theGeometry,
                                        const CaloTopology * theCaloTopology,
                                        const EcalRecHitCollection* theEndcapEcalRecHits,
                                        const reco::BasicClusterCollection* theEndcapBasicClusters,
                                        const reco::SuperClusterCollection* theEndcapSuperClusters,
                                        EcalClusterLazyTools* lazyTools,
                                        const std::map<int,float> & XtalMap,
                                        const std::map<int,float> & XtalMapCurved,
                                        EcalTimePi0TreeContent & myTreeVariables_)
{
  int numberOfSuperClusters = myTreeVariables_.nSuperClusters;
  int numberOfClusters = myTreeVariables_.nClusters ;
  int numberOfXtals = myTreeVariables_.nXtals ;
  
  
  
  //number of superClusters in event (collection = vector!)
  myTreeVariables_.nSuperClusters += theEndcapSuperClusters->size () ;
  myTreeVariables_.nEndcapSuperClusters = theEndcapSuperClusters->size () ;
  
  //loop on all endcap superclusters in event
  for (reco::SuperClusterCollection::const_iterator sclus = theEndcapSuperClusters->begin () ; 
       sclus != theEndcapSuperClusters->end () ; 
       ++sclus) 
    {
  
      int numberOfXtalsInSuperCluster = 0 ;//counter for all xtals in supercluster 
  
      myTreeVariables_.nClustersInSuperCluster[numberOfSuperClusters] = sclus -> clustersSize () ;
      
      if(sclus -> position ().eta () > 0.) 
        myTreeVariables_.superClusterType[numberOfSuperClusters] = 1 ;
      else if(sclus -> position ().eta () < 0.) 
        myTreeVariables_.superClusterType[numberOfSuperClusters] = -1 ;
      

      myTreeVariables_.superClusterEta[numberOfSuperClusters] = sclus -> position ().eta () ;
      myTreeVariables_.superClusterPhi[numberOfSuperClusters] = sclus -> position ().phi () ;
      myTreeVariables_.superClusterX[numberOfSuperClusters] = sclus -> position ().x () ;
      myTreeVariables_.superClusterY[numberOfSuperClusters] = sclus -> position ().y () ;
      myTreeVariables_.superClusterZ[numberOfSuperClusters] = sclus -> position ().z () ;
      myTreeVariables_.superClusterRawEnergy[numberOfSuperClusters] = sclus -> rawEnergy () ;
      myTreeVariables_.superClusterPhiWidth[numberOfSuperClusters] = sclus -> phiWidth () ;
      myTreeVariables_.superClusterEtaWidth[numberOfSuperClusters] = sclus -> etaWidth () ;
      myTreeVariables_.clusterIndexInSuperCluster[numberOfSuperClusters] = numberOfClusters ;
      myTreeVariables_.xtalIndexInSuperCluster[numberOfSuperClusters] = numberOfXtals ;
      //float & energySum = myTreeVariables_.superClusterEnergySum[numberOfSuperClusters];
      //energySum = 0.;
      
      myTreeVariables_.nXtalsInSuperCluster[numberOfSuperClusters]=numberOfXtalsInSuperCluster ;
      numberOfSuperClusters++ ;
  
    } //end endcap supercluster loop

  

  ///////////////////////////////////////////////////////////////////////////////////////
  // independent loop on endcap basic clusters
  float & energySum = myTreeVariables_.superClusterEnergySum[numberOfSuperClusters];
  energySum = 0.;
  for (reco::BasicClusterCollection::const_iterator clus = theEndcapBasicClusters->begin () ; 
       clus != theEndcapBasicClusters->end () ;  
       ++clus) // loop on endcap Bclusters
    {        

         double energy = (clus)->energy () ;
         double phi    = (clus)->phi () ;
         double eta    = (clus)->eta () ;
	 double sinTheta         = fabs( sin( 2 *atan( exp(-1*(clus)->eta()) ) ) );
	 double transverseEnergy = (clus)->energy () * sinTheta;
         double time = -1000.0 ;  // gfdoc: work on this to provide a combination of crystals?
         double ampli = 0. ;
         double secondMin = 0. ;
         double secondTime = -1000. ;
         int numberOfXtalsInCluster=0 ;//counter for all xtals in cluster
         int numXtalsinCluster = 0 ;//xtals in cluster above 3sigma noise  
         EEDetId maxDet ;
         EEDetId secDet ;
         
         std::vector<std::pair<DetId, float> > clusterDetIds = (clus)->hitsAndFractions() ; //get these from the cluster
	 //std::cout << "gf size of ee Bcluster: " << clusterDetIds.size() << std::endl;
         myTreeVariables_.xtalIndexInCluster[numberOfClusters] = numberOfXtals ;


        //loop on xtals in cluster
         for (std::vector<std::pair<DetId, float> >::const_iterator detitr = clusterDetIds.begin () ; 
              detitr != clusterDetIds.end () ; 
              ++detitr)// loop on rechics of endcap basic clusters
           {
             //Here I use the "find" on a digi collection... I have been warned...
             if ( (detitr -> first).det () != DetId::Ecal) 
               { 
                 std::cout << " det is " << (detitr -> first).det () << std::endl ;
                 continue ;
               }
             if ( (detitr -> first).subdetId () != EcalEndcap) 
               {
                 std::cout << " subdet is " << (detitr -> first).subdetId () << std::endl ; 
                 continue ; 
               }
             EcalRecHitCollection::const_iterator thishit = theEndcapEcalRecHits->find ( (detitr -> first) ) ;
             if (thishit == theEndcapEcalRecHits->end ()) 
               {
                 continue ;
               }
             //The checking above should no longer be needed...... 
             //as only those in the cluster would already have rechits..
             
             EcalRecHit myhit = (*thishit) ;
             
             double thisamp = myhit.energy () ;
             if (thisamp > 0.027) //cut on energy->number of crystals in cluster above 3sigma noise
               { 
                 numXtalsinCluster++ ; //xtals in cluster above 3sigma noise  
               }
	     
	     // GFdoc: 'second' keeps track of current cry if it has larger energy than previous
             if (thisamp > secondMin) 
               {
                 secondMin = thisamp ; 
                 secondTime = myhit.time () ; 
                 secDet = (EEDetId) (detitr -> first) ;
               }
             if (secondMin > ampli) 
               {
                 std::swap (ampli, secondMin) ; 
                 std::swap (time, secondTime) ; 
                 std::swap (maxDet, secDet) ;
               }
    
             myTreeVariables_.xtalEnergy[numberOfXtals] = (float) thisamp ; //gf: add uncalRecHits here?
             myTreeVariables_.xtalTime[numberOfXtals] = (float) (myhit.time ()) ;
             myTreeVariables_.xtalHashedIndex[numberOfXtals] = EEDetId (detitr -> first).hashedIndex () ;
	     energySum += (float) thisamp ;
             //MF Lenght evaluation in XTals
             int raw = (detitr -> first).rawId () ;

	     // xtal variables inside an endcap basic cluster 
	      myTreeVariables_.xtalInBCEnergy[numberOfClusters][numberOfXtalsInCluster]=      (float) thisamp;
	      myTreeVariables_.xtalInBCTime[numberOfClusters][numberOfXtalsInCluster]=        (float) (myhit.time ()); 
	      myTreeVariables_.xtalInBCHashedIndex[numberOfClusters][numberOfXtalsInCluster]= EEDetId (detitr -> first).hashedIndex () ; 
	      myTreeVariables_.xtalInBCIeta[numberOfClusters][numberOfXtalsInCluster]=0;
	      myTreeVariables_.xtalInBCIphi[numberOfClusters][numberOfXtalsInCluster]=0;
	      myTreeVariables_.xtalInBCIx[numberOfClusters][numberOfXtalsInCluster]=          EEDetId((detitr -> first)).ix();
	      myTreeVariables_.xtalInBCIy[numberOfClusters][numberOfXtalsInCluster]=          EEDetId((detitr -> first)).iy();
	      myTreeVariables_.xtalInBCFlag[numberOfClusters][numberOfXtalsInCluster]=         myhit.recoFlag(); 

             if (XtalMap.find (raw) != XtalMap.end ())
               myTreeVariables_.xtalTkLength[numberOfXtals] = XtalMap.find (raw)->second ;
             else
               myTreeVariables_.xtalTkLength[numberOfXtals] = -1. ;
              

             if (XtalMapCurved.find (raw) != XtalMapCurved.end ())
               myTreeVariables_.xtalTkLengthCurved[numberOfXtals] = XtalMapCurved.find (raw)->second ;
             else
               myTreeVariables_.xtalTkLengthCurved[numberOfXtals] = -1. ;

              
             numberOfXtals++ ;
             numberOfXtalsInCluster++ ;
             //numberOfXtalsInSuperCluster++ ;
    
           } //end loop on rechics of endcap basic clusters
         
         float E2 = (float)ampli + (float)secondMin ;
         
         //Set some more values for this cluster
         myTreeVariables_.clusterEnergy[numberOfClusters] = (float) energy ;
         myTreeVariables_.clusterTransverseEnergy[numberOfClusters] = (float) transverseEnergy ;
         myTreeVariables_.clusterTime[numberOfClusters] = (float) time ;
         myTreeVariables_.clusterEta[numberOfClusters] = (float) eta ;
         myTreeVariables_.clusterPhi[numberOfClusters] = (float) phi ;
         myTreeVariables_.clusterE1[numberOfClusters] = (float) ampli ;
         myTreeVariables_.clusterE2[numberOfClusters] = (float) E2 ;
         myTreeVariables_.clusterXtalsAbove3Sigma[numberOfClusters] = numXtalsinCluster ;
         myTreeVariables_.clusterXtals[numberOfClusters] = clusterDetIds.size () ;
         myTreeVariables_.clusterMaxId[numberOfClusters] =  maxDet.rawId () ;
         myTreeVariables_.clusterMaxId[numberOfClusters] =  secDet.rawId () ;
         myTreeVariables_.nXtalsInCluster[numberOfClusters]=numberOfXtalsInCluster ;    
	 // std::cout << "gf numberOfXtalsInCluster in ee Bcluster: " << numberOfXtalsInCluster << std::endl;


	 // (basic) cluster shapes for endcap
         myTreeVariables_.clusterE2x2[numberOfClusters] = lazyTools -> e2x2(*(clus));
         //myTreeVariables_.clusterE2x2[numberOfClusters] = lazyTools -> e2x2(*(clus));
         myTreeVariables_.clusterE3x2[numberOfClusters] = lazyTools -> e3x2(*(clus));
         myTreeVariables_.clusterE3x3[numberOfClusters] = lazyTools -> e3x3(*(clus));
         myTreeVariables_.clusterE4x4[numberOfClusters] = lazyTools -> e4x4(*(clus));
         myTreeVariables_.clusterE5x5[numberOfClusters] = lazyTools -> e5x5(*(clus));

         myTreeVariables_.clusterE2x5Right[numberOfClusters] = lazyTools -> e2x5Right(*(clus));
         myTreeVariables_.clusterE2x5Left[numberOfClusters] = lazyTools -> e2x5Left(*(clus));
         myTreeVariables_.clusterE2x5Top[numberOfClusters] = lazyTools -> e2x5Top(*(clus));
         myTreeVariables_.clusterE2x5Bottom[numberOfClusters] = lazyTools -> e2x5Bottom(*(clus));
         //myTreeVariables_.clusterE3x2Ratio[numberOfClusters] = lazyTools -> e3x2Ratio(*(clus));

         myTreeVariables_.clusterCovEtaEta[numberOfClusters] = (lazyTools -> covariances(*(clus)))[0];
         myTreeVariables_.clusterCovPhiPhi[numberOfClusters] = (lazyTools -> covariances(*(clus)))[2];
         myTreeVariables_.clusterCovEtaPhi[numberOfClusters] = (lazyTools -> covariances(*(clus)))[1];

         myTreeVariables_.clusterLat[numberOfClusters] = (lazyTools -> lat(*(clus)))[2];
         myTreeVariables_.clusterEtaLat[numberOfClusters] = (lazyTools -> lat(*(clus)))[0];
         myTreeVariables_.clusterPhiLat[numberOfClusters] = (lazyTools -> lat(*(clus)))[1];

         myTreeVariables_.clusterZernike20[numberOfClusters] = lazyTools -> zernike20(*(clus));
         myTreeVariables_.clusterZernike42[numberOfClusters] = lazyTools -> zernike42(*(clus));
    
         numberOfClusters++ ;

       } //end endcap basic cluster loop
     

  //number of clusters / xtals for the whole event
  myTreeVariables_.nClusters += numberOfClusters ;
  myTreeVariables_.nXtals += numberOfXtals ;

  return ;
} // dumpEndcapClusterInfo  






//// --------------------------------------------------------------------------------------------------------
//// GFdoc muon reco and track length in ECAL
//void EcalTimePi0Tree::dumpMUinfo (const edm::Event& iEvent, const edm::EventSetup& iSetup,
//                                  const CaloGeometry * theGeometry,
//                                  const CaloTopology * theCaloTopology,
//                                  const EcalRecHitCollection* theBarrelEcalRecHits,
//                                  const reco::MuonCollection* theMuons,
//                                  std::map<int,float>& XtalMap,
//                                  std::map<int,float>& XtalMapCurved,
//                                  std::map<int,float>& XtalMapCurved_high,
//                                  std::map<int,float>& XtalMapCurved_low,
//                                  std::map<int,float>& muonCrossedXtalMap,
//                                  std::map<int,float>& muonCrossedXtalCurvedMap,
//                                  EcalTimePi0TreeContent & myTreeVariables_) 
//{
//  //myTreeVariables_.nRecoMuons = theMuons -> size () ;
//  int numberOfRecoMuons = 0 ;
//
//
//  
//  //--------FIXME-----------
//  TrajectoryStateTransform tsTransform;
//  edm::ESHandle<MagneticField> bField;
//  iSetup.get<IdealMagneticFieldRecord>().get(bField);
//     
//  edm::ESHandle<CaloGeometry> pGeometry ;
//  iSetup.get<CaloGeometryRecord> ().get (pGeometry) ;
//  //  const CaloGeometry * theGeometry = pGeometry.product () ;
//  //--------FIXME-----------
//      
//  
//  
//  // begin loop on muons
//  for (reco::MuonCollection::const_iterator MUit = theMuons->begin () ; 
//       MUit != theMuons->end () ; ++MUit)
//    {
//      if(!MUit->isGlobalMuon()) continue;
//
//
//
//      const reco::TrackRef muonGlobalTrack = MUit -> globalTrack () ;
//      const reco::TrackRef muonInnerTrack = MUit -> innerTrack () ;
//      const reco::TrackRef muonOuterTrack = MUit -> outerTrack () ;
//           
//      myTreeVariables_.muonX[numberOfRecoMuons] = (float)muonGlobalTrack->vx () ;
//      myTreeVariables_.muonY[numberOfRecoMuons] = (float)muonGlobalTrack->vy () ;
//      myTreeVariables_.muonZ[numberOfRecoMuons] = (float)muonGlobalTrack->vz () ;
//      myTreeVariables_.muonPhi[numberOfRecoMuons] =  (float)muonGlobalTrack->phi () ;
//      myTreeVariables_.muonEta[numberOfRecoMuons] =  (float)muonGlobalTrack->eta () ;
//      myTreeVariables_.muond0[numberOfRecoMuons] =  (float)muonGlobalTrack->d0 () ;
//      myTreeVariables_.muondz[numberOfRecoMuons] =  (float)muonGlobalTrack->dz () ;
//      myTreeVariables_.muonPx[numberOfRecoMuons] = (float)muonGlobalTrack->px () ;
//      myTreeVariables_.muonPy[numberOfRecoMuons] = (float)muonGlobalTrack->py () ;
//      myTreeVariables_.muonPz[numberOfRecoMuons] = (float)muonGlobalTrack->pz () ;
//      myTreeVariables_.muonP[numberOfRecoMuons] =  (float)muonGlobalTrack->p () ;
//      myTreeVariables_.muonPt[numberOfRecoMuons] =  (float)muonGlobalTrack->pt () ;
//      myTreeVariables_.muonPtError[numberOfRecoMuons] =  (float)muonGlobalTrack->ptError () ;
//      myTreeVariables_.muonCharge[numberOfRecoMuons] =  muonGlobalTrack->charge () ;
//      myTreeVariables_.muonQOverP[numberOfRecoMuons] =  (float)muonGlobalTrack->qoverp () ;
//      myTreeVariables_.muonQOverPError[numberOfRecoMuons] =  (float)muonGlobalTrack->qoverpError () ;
//      myTreeVariables_.muonNChi2[numberOfRecoMuons] =  (float)muonGlobalTrack->normalizedChi2 () ;
//      myTreeVariables_.muonNDof[numberOfRecoMuons] =  muonGlobalTrack->ndof () ;
//      myTreeVariables_.muonNHits[numberOfRecoMuons] = (float)muonGlobalTrack->numberOfValidHits () ;
//      
//      myTreeVariables_.muonInnerHitX[numberOfRecoMuons] = (float)muonGlobalTrack->innerPosition ().x () ;
//      myTreeVariables_.muonInnerHitY[numberOfRecoMuons] = (float)muonGlobalTrack->innerPosition ().y () ;
//      myTreeVariables_.muonInnerHitZ[numberOfRecoMuons] = (float)muonGlobalTrack->innerPosition ().z () ;
//      myTreeVariables_.muonInnerHitEta[numberOfRecoMuons] = (float)muonGlobalTrack->innerPosition ().eta () ;
//      myTreeVariables_.muonInnerHitPhi[numberOfRecoMuons] = (float)muonGlobalTrack->innerPosition ().phi () ;
//      myTreeVariables_.muonInnerHitPx[numberOfRecoMuons] = (float)muonGlobalTrack->innerMomentum ().x () ;
//      myTreeVariables_.muonInnerHitPy[numberOfRecoMuons] = (float)muonGlobalTrack->innerMomentum ().y () ;
//      myTreeVariables_.muonInnerHitPz[numberOfRecoMuons] = (float)muonGlobalTrack->innerMomentum ().z () ;
//      myTreeVariables_.muonInnerHitP[numberOfRecoMuons] = sqrt((float)muonGlobalTrack->innerMomentum ().mag2 ()) ;
//      myTreeVariables_.muonInnerHitPt[numberOfRecoMuons] = sqrt((float)muonGlobalTrack->innerMomentum ().perp2 ()) ;
//
//      myTreeVariables_.muonOuterHitX[numberOfRecoMuons] = (float)muonGlobalTrack->outerPosition ().x () ;
//      myTreeVariables_.muonOuterHitY[numberOfRecoMuons] = (float)muonGlobalTrack->outerPosition ().y () ;
//      myTreeVariables_.muonOuterHitZ[numberOfRecoMuons] = (float)muonGlobalTrack->outerPosition ().z () ;
//      myTreeVariables_.muonOuterHitEta[numberOfRecoMuons] = (float)muonGlobalTrack->outerPosition ().eta () ;
//      myTreeVariables_.muonOuterHitPhi[numberOfRecoMuons] = (float)muonGlobalTrack->outerPosition ().phi () ;
//      myTreeVariables_.muonOuterHitPx[numberOfRecoMuons] = (float)muonGlobalTrack->outerMomentum ().x () ;
//      myTreeVariables_.muonOuterHitPy[numberOfRecoMuons] = (float)muonGlobalTrack->outerMomentum ().y () ;
//      myTreeVariables_.muonOuterHitPz[numberOfRecoMuons] = (float)muonGlobalTrack->outerMomentum ().z () ;
//      myTreeVariables_.muonOuterHitP[numberOfRecoMuons] = (float)muonGlobalTrack->outerP () ;
//      myTreeVariables_.muonOuterHitPt[numberOfRecoMuons] = (float)muonGlobalTrack->outerPt () ;
//
//      myTreeVariables_.muonInnTkInnerHitX[numberOfRecoMuons] = (float)muonInnerTrack->innerPosition ().x () ;
//      myTreeVariables_.muonInnTkInnerHitY[numberOfRecoMuons] = (float)muonInnerTrack->innerPosition ().y () ;
//      myTreeVariables_.muonInnTkInnerHitZ[numberOfRecoMuons] = (float)muonInnerTrack->innerPosition ().z () ;
//      myTreeVariables_.muonInnTkInnerHitEta[numberOfRecoMuons] = (float)muonInnerTrack->innerPosition ().eta () ;
//      myTreeVariables_.muonInnTkInnerHitPhi[numberOfRecoMuons] = (float)muonInnerTrack->innerPosition ().phi () ;
//      myTreeVariables_.muonInnTkInnerHitPx[numberOfRecoMuons] = (float)muonInnerTrack->innerMomentum ().x () ;
//      myTreeVariables_.muonInnTkInnerHitPy[numberOfRecoMuons] = (float)muonInnerTrack->innerMomentum ().y () ;
//      myTreeVariables_.muonInnTkInnerHitPz[numberOfRecoMuons] = (float)muonInnerTrack->innerMomentum ().z () ;
//      myTreeVariables_.muonInnTkInnerHitP[numberOfRecoMuons] = sqrt((float)muonInnerTrack->innerMomentum ().mag2 ()) ;
//      myTreeVariables_.muonInnTkInnerHitPt[numberOfRecoMuons] = sqrt((float)muonInnerTrack->innerMomentum ().perp2 ()) ;
//
//      myTreeVariables_.muonInnTkOuterHitX[numberOfRecoMuons] = (float)muonInnerTrack->outerPosition ().x () ;
//      myTreeVariables_.muonInnTkOuterHitY[numberOfRecoMuons] = (float)muonInnerTrack->outerPosition ().y () ;
//      myTreeVariables_.muonInnTkOuterHitZ[numberOfRecoMuons] = (float)muonInnerTrack->outerPosition ().z () ;
//      myTreeVariables_.muonInnTkOuterHitEta[numberOfRecoMuons] = (float)muonInnerTrack->outerPosition ().eta () ;
//      myTreeVariables_.muonInnTkOuterHitPhi[numberOfRecoMuons] = (float)muonInnerTrack->outerPosition ().phi () ;
//      myTreeVariables_.muonInnTkOuterHitPx[numberOfRecoMuons] = (float)muonInnerTrack->outerMomentum ().x () ;
//      myTreeVariables_.muonInnTkOuterHitPy[numberOfRecoMuons] = (float)muonInnerTrack->outerMomentum ().y () ;
//      myTreeVariables_.muonInnTkOuterHitPz[numberOfRecoMuons] = (float)muonInnerTrack->outerMomentum ().z () ;
//      myTreeVariables_.muonInnTkOuterHitP[numberOfRecoMuons] = (float)muonInnerTrack->outerP () ;
//      myTreeVariables_.muonInnTkOuterHitPt[numberOfRecoMuons] = (float)muonInnerTrack->outerPt () ;
//
//      myTreeVariables_.muonOutTkInnerHitX[numberOfRecoMuons] = (float)muonOuterTrack->innerPosition ().x () ;
//      myTreeVariables_.muonOutTkInnerHitY[numberOfRecoMuons] = (float)muonOuterTrack->innerPosition ().y () ;
//      myTreeVariables_.muonOutTkInnerHitZ[numberOfRecoMuons] = (float)muonOuterTrack->innerPosition ().z () ;
//      myTreeVariables_.muonOutTkInnerHitEta[numberOfRecoMuons] = (float)muonOuterTrack->innerPosition ().eta () ;
//      myTreeVariables_.muonOutTkInnerHitPhi[numberOfRecoMuons] = (float)muonOuterTrack->innerPosition ().phi () ;
//      myTreeVariables_.muonOutTkInnerHitPx[numberOfRecoMuons] = (float)muonOuterTrack->innerMomentum ().x () ;
//      myTreeVariables_.muonOutTkInnerHitPy[numberOfRecoMuons] = (float)muonOuterTrack->innerMomentum ().y () ;
//      myTreeVariables_.muonOutTkInnerHitPz[numberOfRecoMuons] = (float)muonOuterTrack->innerMomentum ().z () ;
//      myTreeVariables_.muonOutTkInnerHitP[numberOfRecoMuons] = sqrt((float)muonOuterTrack->innerMomentum ().mag2 ()) ;
//      myTreeVariables_.muonOutTkInnerHitPt[numberOfRecoMuons] = sqrt((float)muonOuterTrack->innerMomentum ().perp2()) ;
//
//      myTreeVariables_.muonOutTkOuterHitX[numberOfRecoMuons] = (float)muonOuterTrack->outerPosition ().x () ;
//      myTreeVariables_.muonOutTkOuterHitY[numberOfRecoMuons] = (float)muonOuterTrack->outerPosition ().y () ;
//      myTreeVariables_.muonOutTkOuterHitZ[numberOfRecoMuons] = (float)muonOuterTrack->outerPosition ().z () ;
//      myTreeVariables_.muonOutTkOuterHitEta[numberOfRecoMuons] = (float)muonOuterTrack->outerPosition ().eta () ;
//      myTreeVariables_.muonOutTkOuterHitPhi[numberOfRecoMuons] = (float)muonOuterTrack->outerPosition ().phi () ;
//      myTreeVariables_.muonOutTkOuterHitPx[numberOfRecoMuons] = (float)muonOuterTrack->outerMomentum ().x () ;
//      myTreeVariables_.muonOutTkOuterHitPy[numberOfRecoMuons] = (float)muonOuterTrack->outerMomentum ().y () ;
//      myTreeVariables_.muonOutTkOuterHitPz[numberOfRecoMuons] = (float)muonOuterTrack->outerMomentum ().z () ;
//      myTreeVariables_.muonOutTkOuterHitP[numberOfRecoMuons] = (float)muonOuterTrack->outerP () ;
//      myTreeVariables_.muonOutTkOuterHitPt[numberOfRecoMuons] = (float)muonOuterTrack->outerPt () ;
//      
//      
//      
//      // decide whether muon is up or down leg by looking at the outerTrack  
//      if ( (myTreeVariables_.muonOutTkInnerHitY[numberOfRecoMuons] > 0.) &&
//           (myTreeVariables_.muonOutTkOuterHitY[numberOfRecoMuons] > 0.) )
//        myTreeVariables_.muonLeg[numberOfRecoMuons] = 1;
//
//      else if ( (myTreeVariables_.muonOutTkInnerHitY[numberOfRecoMuons] < 0.) &&
//                (myTreeVariables_.muonOutTkOuterHitY[numberOfRecoMuons] < 0.) )
//        myTreeVariables_.muonLeg[numberOfRecoMuons] = -1;
//  
//      else
//        myTreeVariables_.muonLeg[numberOfRecoMuons] = 0;      
//     
//      
//      
//      // define point and direction to look for intersection of muon track with ecal
//      float pointX, pointY, pointZ;
//      float directionX, directionY, directionZ;
//
//      if (myTreeVariables_.muonLeg[numberOfRecoMuons] > 0)
//        {
//          pointX = myTreeVariables_.muonInnTkInnerHitX[numberOfRecoMuons] ;
//          pointY = myTreeVariables_.muonInnTkInnerHitY[numberOfRecoMuons] ;
//          pointZ = myTreeVariables_.muonInnTkInnerHitZ[numberOfRecoMuons] ;
//          
//          directionX = - myTreeVariables_.muonInnTkInnerHitPx[numberOfRecoMuons] ;
//          directionY = - myTreeVariables_.muonInnTkInnerHitPy[numberOfRecoMuons] ;
//          directionZ = - myTreeVariables_.muonInnTkInnerHitPz[numberOfRecoMuons] ;
//        }
//
//      else 
//        {
//          pointX = myTreeVariables_.muonInnTkOuterHitX[numberOfRecoMuons] ;
//          pointY = myTreeVariables_.muonInnTkOuterHitY[numberOfRecoMuons] ;
//          pointZ = myTreeVariables_.muonInnTkOuterHitZ[numberOfRecoMuons] ;
//          
//          directionX = myTreeVariables_.muonInnTkOuterHitPx[numberOfRecoMuons] ;
//          directionY = myTreeVariables_.muonInnTkOuterHitPy[numberOfRecoMuons] ;
//          directionZ = myTreeVariables_.muonInnTkOuterHitPz[numberOfRecoMuons] ;
//        }
//
//      GlobalPoint point (pointX, pointY, pointZ) ;
//      GlobalVector direction (directionX, directionY, directionZ) ;      
//      
//      
//      
//      
//      
//      
//      // compute approximate muon track length in ecal
//      // myTreeVariables_.muonTkLengthInEcalApprox[numberOfRecoMuons] = 
//      //  EcalTrackLengthMeasurement::getApproximateTrackLength (point, direction) ;
//      // this requires tools developed on top of track associator
//      
//      
//      
//      double totalLengthStraight = 0. ;
//      GlobalPoint internalPoint ;
//      GlobalPoint externalPoint ;
//      EcalTrackLengthMeasurement::getDetailedTrackLengthInXtals (XtalMap,
//                                                                 muonCrossedXtalMap,
//                                                                 totalLengthStraight,
//                                                                 internalPoint,
//								 externalPoint,
//								 & (*theGeometry), 
//                                                                 & (*theCaloTopology), 
//                                                                 point, direction, 0.01) ;
//      myTreeVariables_.muonTkLengthInEcalDetail[numberOfRecoMuons] = totalLengthStraight ;
//
//
//      myTreeVariables_.muonTkInternalPointInEcalX[numberOfRecoMuons] = internalPoint.x () ;
//      myTreeVariables_.muonTkInternalPointInEcalY[numberOfRecoMuons] = internalPoint.y () ;
//      myTreeVariables_.muonTkInternalPointInEcalZ[numberOfRecoMuons] = internalPoint.z () ;
//
//      myTreeVariables_.muonTkExternalPointInEcalX[numberOfRecoMuons] = externalPoint.x () ;
//      myTreeVariables_.muonTkExternalPointInEcalY[numberOfRecoMuons] = externalPoint.y () ;
//      myTreeVariables_.muonTkExternalPointInEcalZ[numberOfRecoMuons] = externalPoint.z () ;
//
//
//
//      // compute curved muon track length in ecal and xtal
//
//      FreeTrajectoryState glbInnerState = tsTransform.innerFreeState(*muonGlobalTrack, &*bField);
//      FreeTrajectoryState glbOuterState = tsTransform.outerFreeState(*muonGlobalTrack, &*bField);
//      FreeTrajectoryState innInnerState = tsTransform.innerFreeState(*muonInnerTrack, &*bField);
//      FreeTrajectoryState innOuterState = tsTransform.outerFreeState(*muonInnerTrack, &*bField);
//     
//
//      std::vector<SteppingHelixStateInfo> neckLace;
//      neckLace = calcEcalDeposit(iSetup,
//                                 &glbInnerState, &glbOuterState,
//                                 &innInnerState, &innOuterState,
//                                 trackParameters_, myTreeVariables_.muonLeg[numberOfRecoMuons]);
//      //std::vector<SteppingHelixStateInfo> neckLace_high;
//      //neckLace_high = calcEcalDeposit(iSetup,
//      //                                &glbInnerState, &glbOuterState,
//      //                                &innInnerState, &innOuterState,
//      //                                trackParameters_, myTreeVariables_.muonLeg[numberOfRecoMuons], 1);
//      //std::vector<SteppingHelixStateInfo> neckLace_low;
//      //neckLace_low = calcEcalDeposit(iSetup,
//      //                               &glbInnerState, &glbOuterState,
//      //                               &innInnerState, &innOuterState,
//      //                               trackParameters_, myTreeVariables_.muonLeg[numberOfRecoMuons], -1);
//
//
//
//      
//      //myTreeVariables_.muonNecklaceSize[numberOfRecoMuons] = neckLace.size();      
//      //for(unsigned int pointIt = 0; pointIt < neckLace.size(); ++pointIt)
//      //{
//      //  myTreeVariables_.muonNecklaceX[numberOfRecoMuons][pointIt] = neckLace.at(pointIt).position().x();      
//      //  myTreeVariables_.muonNecklaceY[numberOfRecoMuons][pointIt] = neckLace.at(pointIt).position().y();      
//      //  myTreeVariables_.muonNecklaceZ[numberOfRecoMuons][pointIt] = neckLace.at(pointIt).position().z();      
//      //}
//      
//
//
//      double totalLengthCurved = 0.;
//      GlobalPoint internalPointCurved (0., 0., 0.) ;
//      GlobalPoint externalPointCurved (0., 0., 0.) ;
//      GlobalVector internalMomentumCurved (0., 0., 0.) ;
//      GlobalVector externalMomentumCurved (0., 0., 0.) ;
//      //double totalLengthCurved_high = 0.;
//      //GlobalPoint internalPointCurved_high (0., 0., 0.) ;
//      //GlobalPoint externalPointCurved_high (0., 0., 0.) ;
//      //GlobalVector internalMomentumCurved_high (0., 0., 0.) ;
//      //GlobalVector externalMomentumCurved_high (0., 0., 0.) ;
//      //double totalLengthCurved_low = 0.;
//      //GlobalPoint internalPointCurved_low (0., 0., 0.) ;
//      //GlobalPoint externalPointCurved_low (0., 0., 0.) ;
//      //GlobalVector internalMomentumCurved_low (0., 0., 0.) ;
//      //GlobalVector externalMomentumCurved_low (0., 0., 0.) ;
//
//      if(neckLace.size() > 1)
//      {
//        EcalTrackLengthMeasurement::getDetailedTrackLengthInXtals (XtalMapCurved,
//                                                                   muonCrossedXtalCurvedMap,
//                                                                   totalLengthCurved,
//                                                                   internalPointCurved,
//                                                                   externalPointCurved,
//                                                                   internalMomentumCurved,
//                                                                   externalMomentumCurved,
//                                                                   & (*theGeometry), 
//                                                                   & (*theCaloTopology), 
//                                                                   neckLace) ;
//        myTreeVariables_.muonTkLengthInEcalDetailCurved[numberOfRecoMuons] = totalLengthCurved ;
//
//        myTreeVariables_.muonTkInternalPointInEcalCurvedX[numberOfRecoMuons] = internalPointCurved.x () ;
//        myTreeVariables_.muonTkInternalPointInEcalCurvedY[numberOfRecoMuons] = internalPointCurved.y () ;
//        myTreeVariables_.muonTkInternalPointInEcalCurvedZ[numberOfRecoMuons] = internalPointCurved.z () ;
//
//        myTreeVariables_.muonTkExternalPointInEcalCurvedX[numberOfRecoMuons] = externalPointCurved.x () ;
//        myTreeVariables_.muonTkExternalPointInEcalCurvedY[numberOfRecoMuons] = externalPointCurved.y () ;
//        myTreeVariables_.muonTkExternalPointInEcalCurvedZ[numberOfRecoMuons] = externalPointCurved.z () ;
//      
//        myTreeVariables_.muonTkInternalPointInEcalCurvedPx[numberOfRecoMuons] = internalMomentumCurved.x () ;
//        myTreeVariables_.muonTkInternalPointInEcalCurvedPy[numberOfRecoMuons] = internalMomentumCurved.y () ;
//        myTreeVariables_.muonTkInternalPointInEcalCurvedPz[numberOfRecoMuons] = internalMomentumCurved.z () ;
//
//        myTreeVariables_.muonTkExternalPointInEcalCurvedPx[numberOfRecoMuons] = externalMomentumCurved.x () ;
//        myTreeVariables_.muonTkExternalPointInEcalCurvedPy[numberOfRecoMuons] = externalMomentumCurved.y () ;
//        myTreeVariables_.muonTkExternalPointInEcalCurvedPz[numberOfRecoMuons] = externalMomentumCurved.z () ;
//      }
//
//      else
//        myTreeVariables_.muonTkLengthInEcalDetailCurved[numberOfRecoMuons] = -1. ;
//      
//
//
//      //if(neckLace_high.size() > 1)
//      //{
//      //  EcalTrackLengthMeasurement::getDetailedTrackLengthInXtals (XtalMapCurved_high,
//      //                                                             totalLengthCurved_high,
//      //                                                             internalPointCurved_high,
//      //                                                             externalPointCurved_high,
//      //                                                             internalMomentumCurved_high,
//      //                                                             externalMomentumCurved_high,
//      //                                                             & (*theGeometry), 
//      //                                                             & (*theCaloTopology), 
//      //                                                             neckLace_high) ;
//      //  myTreeVariables_.muonTkLengthInEcalDetailCurved_high[numberOfRecoMuons] = totalLengthCurved_high ;
//      //}
//
//      //else
//      //  myTreeVariables_.muonTkLengthInEcalDetailCurved_high[numberOfRecoMuons] = -1. ;
//      
//
//
//      //if(neckLace_low.size() > 1)
//      //{
//      //  EcalTrackLengthMeasurement::getDetailedTrackLengthInXtals (XtalMapCurved_low,
//      //                                                             totalLengthCurved_low,
//      //                                                             internalPointCurved_low,
//      //                                                             externalPointCurved_low,
//      //                                                             internalMomentumCurved_low,
//      //                                                             externalMomentumCurved_low,
//      //                                                             & (*theGeometry), 
//      //                                                             & (*theCaloTopology), 
//      //                                                             neckLace_low) ;
//      //  myTreeVariables_.muonTkLengthInEcalDetailCurved_low[numberOfRecoMuons] = totalLengthCurved_low ;
//      //}
//
//      //else
//      //  myTreeVariables_.muonTkLengthInEcalDetailCurved_low[numberOfRecoMuons] = -1. ;
//      
//      
//      
//
//      int xtalCounter = 0;
//      myTreeVariables_.nMuonCrossedXtals[numberOfRecoMuons] = muonCrossedXtalMap.size();
//      for(std::map<int, float>::const_iterator mapIt = muonCrossedXtalMap.begin();
//          mapIt != muonCrossedXtalMap.end(); ++mapIt)
//      {
//        EBDetId dummyDetId(mapIt -> first);
//
//        myTreeVariables_.muonCrossedXtalHashedIndex[numberOfRecoMuons][xtalCounter] = dummyDetId.hashedIndex();
//        myTreeVariables_.muonCrossedXtalTkLength[numberOfRecoMuons][xtalCounter] = mapIt -> second;
//
//        ++xtalCounter;
//      }
//
//      int xtalCounterCurved = 0;
//      myTreeVariables_.nMuonCrossedXtalsCurved[numberOfRecoMuons] = muonCrossedXtalCurvedMap.size();
//      for(std::map<int, float>::const_iterator mapIt = muonCrossedXtalCurvedMap.begin();
//          mapIt != muonCrossedXtalCurvedMap.end(); ++mapIt)
//      {
//        EBDetId dummyDetId(mapIt -> first);       
//
//	myTreeVariables_.muonCrossedXtalHashedIndexCurved[numberOfRecoMuons][xtalCounterCurved] = dummyDetId.hashedIndex();
//        myTreeVariables_.muonCrossedXtalTkLengthCurved[numberOfRecoMuons][xtalCounterCurved] = mapIt -> second;
//
//        ++xtalCounterCurved;
//      }
//      
//      
//      
//      
//
//          
//      ++numberOfRecoMuons ;
//    } // begin loop on muons
//
//  myTreeVariables_.nRecoMuons = numberOfRecoMuons;
//  
//  return ;
//  
//} // dumpMUinfo
//


// GFdoc GT information, at present anf +-1 bx
void
EcalTimePi0Tree::dump3Ginfo (const edm::Event& iEvent, 
                             const edm::EventSetup& iSetup,
                             EcalTimePi0TreeContent & myTreeVariables_)
{
  std::vector<bool> l1Triggers = determineTriggers (iEvent, iSetup);      
  myTreeVariables_.isECALL1 = l1Triggers[4];
  myTreeVariables_.isHCALL1 = l1Triggers[3];
  myTreeVariables_.isRPCL1 = l1Triggers[2];
  myTreeVariables_.isCSCL1 = l1Triggers[1];
  myTreeVariables_.isDTL1 = l1Triggers[0];
  
  int nbunches = 1;
  for(int bx = -nbunches; bx <=nbunches; bx++){
    std::vector<bool> l1Triggers = determineTriggers (iEvent, iSetup, bx);      
    myTreeVariables_.isECALL1Bx[ bx + nbunches] = l1Triggers[4];
    myTreeVariables_.isHCALL1Bx[ bx + nbunches] = l1Triggers[3];
    myTreeVariables_.isRPCL1Bx [ bx + nbunches] = l1Triggers[2];
    myTreeVariables_.isCSCL1Bx [ bx + nbunches] = l1Triggers[1];
    myTreeVariables_.isDTL1Bx  [ bx + nbunches] = l1Triggers[0];
  }
  //cout << "isECALL1Bx: "
  //     << myTreeVariables_.isECALL1Bx[0]
  //     << myTreeVariables_.isECALL1Bx[1]
  //     << myTreeVariables_.isECALL1Bx[2] << endl;
  //cout << "isHCALL1Bx: " 
  //     << myTreeVariables_.isHCALL1Bx[0]
  //     << myTreeVariables_.isHCALL1Bx[1]
  //     << myTreeVariables_.isHCALL1Bx[2] << endl;
  //cout << "isRPCL1Bx: "  
  //     << myTreeVariables_.isRPCL1Bx [0]
  //     << myTreeVariables_.isRPCL1Bx [1]
  //     << myTreeVariables_.isRPCL1Bx [2] << endl;
  //cout << "isCSCL1Bx: "  
  //     << myTreeVariables_.isCSCL1Bx [0]
  //     << myTreeVariables_.isCSCL1Bx [1]
  //     << myTreeVariables_.isCSCL1Bx [2] << endl;
  //cout << "isDTL1Bx: "  
  //     << myTreeVariables_.isDTL1Bx  [0]
  //     << myTreeVariables_.isDTL1Bx  [1]
  //     << myTreeVariables_.isDTL1Bx  [2] << endl;
 
  return ;
} //PG dump3Ginfo  

void EcalTimePi0Tree::dumpHcalInfo(const CaloGeometry * theGeometry, 
				   const edm::Event& iEvent, 
				   const edm::EventSetup& iSetup,
				   EcalTimePi0TreeContent & myTreeVariables_) 
{
  // Code from F. Santanastasio et al RecoHcal/HcalProm/src/HcalPromBremss.cc
  
  // HCAL rechits
  Handle < HBHERecHitCollection > hbhe;
  iEvent.getByLabel("hbhereco", hbhe);
  const HBHERecHitCollection Hithbhe = *(hbhe.product());
  const HBHERecHitCollection *HBHERecHits = 0;
  iEvent.getByLabel("hbhereco", hbhe);
  HBHERecHits = hbhe.product();   // get a ptr to the product

  // Handle < HORecHitCollection > horh;
  // iEvent.getByLabel("horeco", horh);
  // const HORecHitCollection Hitho = *(horh.product());

  Handle < CaloTowerCollection > calo;
  iEvent.getByLabel("towerMaker", calo);
  const CaloTowerCollection *calohbhe = calo.product();
  
  // -------------------------------------------------------------------------------------------------------------
  
  // %%%%%%%%%%%%%%%% maybe can be removed  
  // Timing plots for different triggers
  myTreeVariables_.hbNRecHits = 0;
  for (HBHERecHitCollection::const_iterator hhit = Hithbhe.begin(); 
       myTreeVariables_.hbNRecHits < MAXHCALRECHITS && hhit != Hithbhe.end(); hhit++) {

    if(hhit->id().subdet()!=HcalBarrel)continue; // %%%%%%%%% ONLY BARREL

    if ( hhit->energy() < hbTreshold_ ) continue; // %%%%%%%%% THRESHOLD ON RECHIT ENERGY SET BY DEFAULT AT 1
    
    GlobalPoint pos = theGeometry->getPosition((*hhit).detid());
    
    myTreeVariables_.hbRecHitDetId [ myTreeVariables_.hbNRecHits ] = (*hhit).detid();
    myTreeVariables_.hbRecHitEta   [ myTreeVariables_.hbNRecHits ] = pos.eta();
    myTreeVariables_.hbRecHitPhi   [ myTreeVariables_.hbNRecHits ] = pos.phi();
    myTreeVariables_.hbRecHitE     [ myTreeVariables_.hbNRecHits ] = (*hhit).energy();
    myTreeVariables_.hbRecHitTime  [ myTreeVariables_.hbNRecHits ] = (*hhit).time();
    myTreeVariables_.hbNRecHits++;
    
  }
  //%%%%%%%%%%%
    
  // for (HORecHitCollection::const_iterator hhit = Hitho.begin(); hhit != Hitho.end(); hhit++) {
  //   if (hhit->energy() > 0.6) {
  //     h_ho_rechit_energy->Fill(hhit->energy());
  // 	h_ho_eta_phi->Fill((hhit->id()).ieta(), (hhit->id()).iphi());
  // 	if (hhit->energy() > maxhorec) {
  // 	  maxhorec = hhit->energy();
  // 	  maxhorec_ETA = hhit->id().ieta();
  // 	  maxhorec_PHI = hhit->id().iphi();
  // 	}
  // 	if (hhit->energy() > next_to_maxhorec && hhit->energy() < maxhorec) {
  // 	  next_to_maxhorec = hhit->energy();
  // 	  next_to_maxhorec_ETA = hhit->id().ieta();
  // 	  next_to_maxhorec_PHI = hhit->id().iphi();
  // 	}
  //   }
  // }

  myTreeVariables_.nCaloTowers=0;
  for (CaloTowerCollection::const_iterator kal = calohbhe->begin();
       myTreeVariables_.nCaloTowers < MAXCALOTOWERS && kal != calohbhe->end(); kal++) 
    {
      if( kal->emEnergy() + kal->hadEnergy() + kal->outerEnergy() <= 0) continue;
      myTreeVariables_.caloTowerEmEnergy    [myTreeVariables_.nCaloTowers] = kal->emEnergy();
      myTreeVariables_.caloTowerHadEnergy   [myTreeVariables_.nCaloTowers] = kal->hadEnergy();
      myTreeVariables_.caloTowerOuterEnergy [myTreeVariables_.nCaloTowers] = kal->outerEnergy();
      myTreeVariables_.caloTowerEmEta       [myTreeVariables_.nCaloTowers] = kal->emPosition().eta();
      myTreeVariables_.caloTowerEmPhi       [myTreeVariables_.nCaloTowers] = kal->emPosition().phi();
      myTreeVariables_.caloTowerHadEta      [myTreeVariables_.nCaloTowers] = kal->hadPosition().eta();
      myTreeVariables_.caloTowerHadPhi      [myTreeVariables_.nCaloTowers] = kal->hadPosition().phi();
      myTreeVariables_.nCaloTowers++;
      if(myTreeVariables_.nCaloTowers>MAXCALOTOWERS) break;
    }
  
} // dumpHcalInfo



// -------------------------------------------------------------------------------------------------------------

//PG stolen from Toyoko's EcalCosmicHists, PGP changed it (added Bx information)
std::vector<bool> 
EcalTimePi0Tree::determineTriggers (const edm::Event& iEvent, const edm::EventSetup& eventSetup, int Bx)
{
  std::vector<bool> l1Triggers; //DT,CSC,RPC,HCAL,ECAL
                                //0 , 1 , 2 , 3  , 4
  for(int i=0;i<5;i++)
    l1Triggers.push_back(false);
  
  // get the GMTReadoutCollection
  Handle<L1MuGMTReadoutCollection> gmtrc_handle; 
  iEvent.getByLabel(l1GMTReadoutRecTag_,gmtrc_handle);
  L1MuGMTReadoutCollection const* gmtrc = gmtrc_handle.product();
  if (!(gmtrc_handle.isValid())) 
    {
      LogWarning("EcalCosmicsHists") << "l1MuGMTReadoutCollection" << " not available";
      return l1Triggers;
    }
  // get hold of L1GlobalReadoutRecord
  Handle<L1GlobalTriggerReadoutRecord> L1GTRR;
  iEvent.getByLabel(l1GMTReadoutRecTag_,L1GTRR);
  
  //Ecal
  edm::ESHandle<L1GtTriggerMenu> menuRcd;
  eventSetup.get<L1GtTriggerMenuRcd>().get(menuRcd) ;
  const L1GtTriggerMenu* menu = menuRcd.product();
  edm::Handle< L1GlobalTriggerReadoutRecord > gtRecord;
  iEvent.getByLabel( edm::InputTag("gtDigis"), gtRecord);
  // Get dWord after masking disabled bits
  const DecisionWord dWord = gtRecord->decisionWord(Bx);

  //cout << "Bx: " << Bx << " TriggerWord: " << endl;
  //copy (dWord.begin(), dWord.end(), std::ostream_iterator<bool>(std::cout, ""));
  //cout << endl;
  
  bool l1SingleEG1 = menu->gtAlgorithmResult("L1_SingleEG1", dWord);
  bool l1SingleEG5 = menu->gtAlgorithmResult("L1_SingleEG5", dWord);
  bool l1SingleEG8 = menu->gtAlgorithmResult("L1_SingleEG8", dWord);
  bool l1SingleEG10 = menu->gtAlgorithmResult("L1_SingleEG10", dWord);
  bool l1SingleEG12 = menu->gtAlgorithmResult("L1_SingleEG12", dWord);
  bool l1SingleEG15 = menu->gtAlgorithmResult("L1_SingleEG15", dWord);
  bool l1SingleEG20 = menu->gtAlgorithmResult("L1_SingleEG20", dWord);
  bool l1SingleEG25 = menu->gtAlgorithmResult("L1_SingleEG25", dWord);
  bool l1DoubleNoIsoEGBTBtight = menu->gtAlgorithmResult("L1_DoubleNoIsoEG_BTB_tight", dWord);
  bool l1DoubleNoIsoEGBTBloose = menu->gtAlgorithmResult("L1_DoubleNoIsoEG_BTB_loose ", dWord);
  bool l1DoubleNoIsoEGTopBottom = menu->gtAlgorithmResult("L1_DoubleNoIsoEGTopBottom", dWord);
  bool l1DoubleNoIsoEGTopBottomCen  = menu->gtAlgorithmResult("L1_DoubleNoIsoEGTopBottomCen", dWord);
  bool l1DoubleNoIsoEGTopBottomCen2  = menu->gtAlgorithmResult("L1_DoubleNoIsoEGTopBottomCen2", dWord);
  bool l1DoubleNoIsoEGTopBottomCenVert  = menu->gtAlgorithmResult("L1_DoubleNoIsoEGTopBottomCenVert", dWord);
  
  l1Triggers[4] = l1SingleEG1 || l1SingleEG5 || l1SingleEG8 || l1SingleEG10 || l1SingleEG12 || l1SingleEG15
    || l1SingleEG20 || l1SingleEG25 || l1DoubleNoIsoEGBTBtight || l1DoubleNoIsoEGBTBloose
    || l1DoubleNoIsoEGTopBottom || l1DoubleNoIsoEGTopBottomCen || l1DoubleNoIsoEGTopBottomCen2
    || l1DoubleNoIsoEGTopBottomCenVert;
  //cout << "l1Triggers[4]" << l1Triggers[4] << endl;


  std::vector<L1MuGMTReadoutRecord> gmt_records = gmtrc->getRecords();
  std::vector<L1MuGMTReadoutRecord>::const_iterator igmtrr;

  for(igmtrr=gmt_records.begin(); igmtrr!=gmt_records.end(); igmtrr++) {

    std::vector<L1MuRegionalCand>::const_iterator iter1;
    std::vector<L1MuRegionalCand> rmc;
    
    //DT triggers
    int idt = 0;
    rmc = igmtrr->getDTBXCands();
    for(iter1=rmc.begin(); iter1!=rmc.end(); iter1++) {
      if ( !(*iter1).empty() ) {
        idt++;
      }
    }
    //if(idt>0) std::cout << "Found " << idt << " valid DT candidates in bx wrt. L1A = " 
    //  << igmtrr->getBxInEvent() << std::endl;
    if(igmtrr->getBxInEvent()==Bx && idt>0) l1Triggers[0] = true;
    //cout << "l1Triggers[0]" << l1Triggers[0] << endl;
    //RPC triggers
    int irpcb = 0;
    rmc = igmtrr->getBrlRPCCands();
    for(iter1=rmc.begin(); iter1!=rmc.end(); iter1++) {
      if ( !(*iter1).empty() ) {
        irpcb++;
      }
    }
    //if(irpcb>0) std::cout << "Found " << irpcb << " valid RPC candidates in bx wrt. L1A = " 
    //  << igmtrr->getBxInEvent() << std::endl;
    if(igmtrr->getBxInEvent()==Bx && irpcb>0) l1Triggers[2] = true;
    //cout << "l1Triggers[2]" << l1Triggers[2] << endl;

    //CSC Triggers
    int icsc = 0;
    rmc = igmtrr->getCSCCands();
    for(iter1=rmc.begin(); iter1!=rmc.end(); iter1++) {
      if ( !(*iter1).empty() ) {
        icsc++;
      }
    }
    //if(icsc>0) std::cout << "Found " << icsc << " valid CSC candidates in bx wrt. L1A = " 
    //  << igmtrr->getBxInEvent() << std::endl;
    if(igmtrr->getBxInEvent()==Bx && icsc>0) l1Triggers[1] = true;
  }
  //cout << "l1Triggers[1]" << l1Triggers[1] << endl;

  L1GlobalTriggerReadoutRecord const* gtrr = L1GTRR.product();
  
  for(int ibx=-1; ibx<=1; ibx++) {
    bool hcal_top = false;
    bool hcal_bot = false;
    const L1GtPsbWord psb = gtrr->gtPsbWord(0xbb0d,ibx);
    std::vector<int> valid_phi;
    if((psb.aData(4)&0x3f) >= 1) {valid_phi.push_back( (psb.aData(4)>>10)&0x1f ); }
    if((psb.bData(4)&0x3f) >= 1) {valid_phi.push_back( (psb.bData(4)>>10)&0x1f ); }
    if((psb.aData(5)&0x3f) >= 1) {valid_phi.push_back( (psb.aData(5)>>10)&0x1f ); }
    if((psb.bData(5)&0x3f) >= 1) {valid_phi.push_back( (psb.bData(5)>>10)&0x1f ); }
    std::vector<int>::const_iterator iphi;
    for(iphi=valid_phi.begin(); iphi!=valid_phi.end(); iphi++) {
      //std::cout << "Found HCAL mip with phi=" << *iphi << " in bx wrt. L1A = " << ibx << std::endl;
      if(*iphi<9) hcal_top=true;
      if(*iphi>8) hcal_bot=true;
    }
    if(ibx==0 && hcal_top && hcal_bot) l1Triggers[3]=true;
  }
  //cout << "l1Triggers[3]" << l1Triggers[3] << endl;

  /*edm::LogInfo("EcalCosmicsHists") << "**** Trigger SourceSource ****";
  if(l1Triggers[0]) edm::LogInfo("EcalCosmicsHists") << "DT";
  if(l1Triggers[2]) edm::LogInfo("EcalCosmicsHists") << "RPC";
  if(l1Triggers[1]) edm::LogInfo("EcalCosmicsHists") << "CSC";
  if(l1Triggers[3]) edm::LogInfo("EcalCosmicsHists") << "HCAL";
  if(l1Triggers[4]) edm::LogInfo("EcalCosmicsHists") << "ECAL";
  edm::LogInfo("EcalCosmicsHists") << "************************";
  */
  //cout << "l1Triggers: "; 
  //for(int i=0;i<5;i++)
  //  cout << l1Triggers[i];
  //cout << endl;
  return l1Triggers;
}
// -------------------------------------------------------------------------------------------------------------


/* for the future?
void
EcalTimePi0Tree::dumpDCCHeaders ()
{
  Handle<EcalRawDataCollection> DCCHeaders;
  iEvent.getByLabel("ecalEBunpacker", DCCHeaders);
  if(!DCCHeaders.isValid())
    LogWarning("EcalCosmicsHists") << "DCC headers not available";

  //make the bx histos right here
  //TODO: Right now we are filling histos for errors...
  int orbit = -100;
  int bx = -100;
  int runType = -100;

  for(EcalRawDataCollection::const_iterator headerItr= DCCHeaders->begin();headerItr != DCCHeaders->end(); 
      ++headerItr) {
    EcalDCCHeaderBlock::EcalDCCEventSettings settings = headerItr->getEventSettings();
    int myorbit = headerItr->getOrbit();
    int mybx = headerItr->getBX();
    int myRunType = headerItr->getRunType();
    int FEDid = headerItr->fedId();
    TH2F* dccRuntypeHist = FEDsAndDCCRuntypeVsBxHists_[FEDid];
    if(dccRuntypeHist==0)
    {
      initHists(FEDid);
      dccRuntypeHist = FEDsAndDCCRuntypeVsBxHists_[FEDid];
    }
    dccRuntypeHist->Fill(mybx,myRunType);
    
    if (bx == -100)
    {
      bx = mybx;
    }
    else if (bx != mybx)
    {
      LogWarning("EcalCosmicsHists") << "This header has a conflicting bx OTHERS were " << bx << " here " << mybx;
      dccBXErrorByFEDHist_->Fill(headerItr->fedId());
      if(bx != -100)
      {
        dccErrorVsBxHist_->Fill(bx,0);
      }
    }

    if (runType == -100)
    {
      runType = myRunType;
    }
    else if (runType != myRunType)
    {
      LogWarning("EcalCosmicsHists") << "This header has a conflicting runType OTHERS were " << bx << " here " << mybx;
      dccRuntypeErrorByFEDHist_->Fill(headerItr->fedId());
      if(bx != -100)
        dccErrorVsBxHist_->Fill(bx,2);
    }
    
    
    if (orbit == -100)
    {
      orbit = myorbit;
    }
    else if (orbit != myorbit)
    {
      LogWarning("EcalCosmicsHists") << "This header has a conflicting orbit; OTHERS were " << orbit << " here " << myorbit;
      dccOrbitErrorByFEDHist_->Fill(headerItr->fedId());
      if(bx != -100)
        dccErrorVsBxHist_->Fill(bx,1);
    }
  }
  
} //PG dumpDCCHeaders
*/

