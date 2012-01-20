// -*- C++ -*-
//
// Package:   EcalTimePhyTreeMaker
// Class:     EcalTimePhyTreeMaker
//
/**\class EcalTimePhyTreeMaker EcalTimePhyTreeMaker.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
 */
//
// Authors:                   Shih-Chuan Kao, Giovanni Franzoni (UMN)
//         Created:  Mo Jul 14 5:46:22 CEST 2008
// $Id: EcalTimePhyTreeMaker.cc,v 1.10 2012/01/17 05:18:09 sckao Exp $
//
//

 
#include "ECALTime/EcalTimePi0/plugins/EcalTimePhyTreeMaker.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalTools.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbRecord.h"

#include "CondFormats/EcalObjects/interface/EcalADCToGeVConstant.h"
#include "CondFormats/DataRecord/interface/EcalADCToGeVConstantRcd.h"

// vertex stuff
#include <DataFormats/VertexReco/interface/VertexFwd.h>
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"


#include <vector>
#include <iostream>
#include <iterator>

using namespace cms ;
using namespace edm ;
using namespace std ;

#define FILL_SHAPE_VARS 0

EcalTimePhyTreeMaker::EcalTimePhyTreeMaker (const edm::ParameterSet& iConfig) :
  barrelEcalRecHitCollection_              (iConfig.getParameter<edm::InputTag> ("barrelEcalRecHitCollection")),
  endcapEcalRecHitCollection_              (iConfig.getParameter<edm::InputTag> ("endcapEcalRecHitCollection")),
  barrelBasicClusterCollection_            (iConfig.getParameter<edm::InputTag> ("barrelBasicClusterCollection")),//GF
  endcapBasicClusterCollection_            (iConfig.getParameter<edm::InputTag> ("endcapBasicClusterCollection")),
  barrelSuperClusterCollection_            (iConfig.getParameter<edm::InputTag> ("barrelSuperClusterCollection")),
  endcapSuperClusterCollection_            (iConfig.getParameter<edm::InputTag> ("endcapSuperClusterCollection")),
  muonCollection_                          (iConfig.getParameter<edm::InputTag> ("muonCollection")),
  JetSource_                            (iConfig.getParameter<edm::InputTag> ("JetSource")),
  METSource_                            (iConfig.getParameter<edm::InputTag> ("METSource")),
  MuonSource_                           (iConfig.getParameter<edm::InputTag> ("MuonSource")),
  ElectronSource_                       (iConfig.getParameter<edm::InputTag> ("ElectronSource")),
  PhotonSource_                         (iConfig.getParameter<edm::InputTag> ("PhotonSource")),
  triggerSource_                        (iConfig.getParameter<edm::InputTag> ("triggerSource")),
  vertexCollection_                     (iConfig.getParameter<edm::InputTag> ("vertexCollection")),
  l1GMTReadoutRecTag_   (iConfig.getParameter<std::string> ("L1GlobalReadoutRecord")),
  gtRecordCollectionTag_ (iConfig.getParameter<std::string> ("GTRecordCollection")),
  runNum_               (iConfig.getParameter<int> ("runNum")),
  minEtEB_              (iConfig.getParameter<double> ("minEtEB")),
  minEtEE_              (iConfig.getParameter<double> ("minEtEE")),
  jetCuts_              (iConfig.getParameter<std::vector<double> >("jetCuts")),
  photonCuts_           (iConfig.getParameter<std::vector<double> >("photonCuts")),
  photonIso_            (iConfig.getParameter<std::vector<double> >("photonIso")),
  electronCuts_         (iConfig.getParameter<std::vector<double> >("electronCuts")),
  muonCuts_             (iConfig.getParameter<std::vector<double> >("muonCuts")),
  fileName_             (iConfig.getParameter<std::string> ("fileName")),
  triggerName_          (iConfig.getParameter<std::string> ("triggerName")),
  doTimeVSAmpliCorrection_(iConfig.getParameter<bool> ("doTimeVSAmpliCorrection")),
  naiveId_ (0)              
{

  // Create File
  fileName_ += "_"+intToString (runNum_)+".root" ;
  file_ = new TFile (fileName_.c_str () , "RECREATE") ;
  file_->cd () ;

  // Initialize Tree
  tree_ = new TTree ( "EcalTimeAnalysis","EcalTimeAnalysis" ) ;
  setBranches (tree_, myTreeVariables_) ;
  
  // initialize the time corrector
  theTimeCorrector_.initEB("EB");
  theTimeCorrector_.initEE("EElow");

}


// -----------------------------------------------------------------------------------------


EcalTimePhyTreeMaker::~EcalTimePhyTreeMaker ()
{
}


// -----------------------------------------------------------------------------------------


void EcalTimePhyTreeMaker::analyze (const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  ++naiveId_ ;
    

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
  const EcalRecHitCollection* theBarrelEcalRecHits = 0;
  if( iEvent.getByLabel (barrelEcalRecHitCollection_, pBarrelEcalRecHits) && pBarrelEcalRecHits.isValid ())
    {
      theBarrelEcalRecHits = pBarrelEcalRecHits.product () ;   
    }
  if (! (pBarrelEcalRecHits.isValid ()) )
    {
      LogWarning ("EcalTimePhyTreeMaker") << barrelEcalRecHitCollection_ 
				     << " not available" ;
      return ;
    }
  
  
  // Ecal endcap RecHits
  edm::Handle<EcalRecHitCollection> pEndcapEcalRecHits ;
  const EcalRecHitCollection* theEndcapEcalRecHits = 0;
  if( iEvent.getByLabel (endcapEcalRecHitCollection_, pEndcapEcalRecHits) && pEndcapEcalRecHits.isValid ())
    {
      theEndcapEcalRecHits = pEndcapEcalRecHits.product () ;   
    }
  
  if (! (pEndcapEcalRecHits.isValid ()))
    {
      LogWarning ("EcalTimePhyTreeMaker") << endcapEcalRecHitCollection_ 
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
      LogWarning ("EcalTimePhyTreeMaker") << barrelSuperClusterCollection_ 
                                     << " not available" ;
      return ;
    } else {
    //std::cout << "EB superclusters: " << barrelSuperClusterCollection_ << std::endl; 
  }
  
  // Endcap SuperClusters
  edm::Handle<reco::SuperClusterCollection> pEndcapSuperClusters ;
  iEvent.getByLabel (endcapSuperClusterCollection_, pEndcapSuperClusters) ;
  const reco::SuperClusterCollection* theEndcapSuperClusters = pEndcapSuperClusters.product () ;
  
  if (! (pEndcapSuperClusters.isValid ()) )
    {
      LogWarning ("EcalTimePhyTreeMaker") << endcapSuperClusterCollection_ 
                                     << " not available" ;
      return ;
    }
  
  
  // Barrel BasicClusters
  edm::Handle<reco::BasicClusterCollection> pBarrelBasicClusters ;
  iEvent.getByLabel (barrelBasicClusterCollection_, pBarrelBasicClusters) ;
  const reco::BasicClusterCollection* theBarrelBasicClusters = pBarrelBasicClusters.product () ;//GF

  if (! (pBarrelBasicClusters.isValid ()) )
    {
      LogWarning ("EcalTimePhyTreeMaker") << barrelBasicClusterCollection_ 
                                     << " not available" ;
      return ;
    }else {
    //std::cout << "EB basiclusters: " << barrelBasicClusterCollection_ << std::endl; 
  }

  // Endcap BasicClusters
  edm::Handle<reco::BasicClusterCollection> pEndcapBasicClusters ;
  iEvent.getByLabel (endcapBasicClusterCollection_, pEndcapBasicClusters) ;
  const reco::BasicClusterCollection* theEndcapBasicClusters = pEndcapBasicClusters.product () ;
  
  if (! (pEndcapBasicClusters.isValid ()) )
    {
      LogWarning ("EcalTimePhyTreeMaker") << endcapBasicClusterCollection_ 
                                     << " not available" ;
      return ;
    }
 
  Handle<reco::VertexCollection> recVtxs;
  iEvent.getByLabel(vertexCollection_, recVtxs);
  const reco::VertexCollection * theRecVtxs = recVtxs.product();

  // ClusterShapes
  EcalClusterLazyTools* lazyTools = new EcalClusterLazyTools(iEvent, iSetup, barrelEcalRecHitCollection_, endcapEcalRecHitCollection_);

  // Xtal - TkLength map
  std::map<int,float> XtalMap ;
  std::map<int,float> XtalMapCurved ;
  std::map<int,float> XtalMapCurved_high ;
  std::map<int,float> XtalMapCurved_low ;

  // Trigger Selection
  //std::cout<<" fk1 "<< std::endl ;
  //int HLTCut = HLTSelection( iEvent ) ; 
    int HLTCut = HLTSelection( iEvent, 75, "HLT_Photon", "_CaloIdVL_IsoL" ) ;


  if ( HLTCut >= 75 ) {

     // GFdoc initialize variables to 0/false
     initializeBranches(tree_, myTreeVariables_);
  
     // Initialize Counter for SCluster and Cluster
     numberOfSuperClusters = 0;
     numberOfClusters = 0;
     sclist.clear() ;

     bool passed = dumpEvtObjectInfo( iEvent ) ;
     if ( passed ) {

        myTreeVariables_.bx          = iEvent.bunchCrossing();
	myTreeVariables_.lumiSection = iEvent.id().luminosityBlock();
	myTreeVariables_.unixTime    = iEvent.eventAuxiliary().time().unixTime();
	myTreeVariables_.orbit       = iEvent.orbitNumber();

	myTreeVariables_.runId         = iEvent.id ().run () ;
	myTreeVariables_.eventId       = iEvent.id ().event () ;
	myTreeVariables_.eventNaiveId  = naiveId_ ;
	myTreeVariables_.timeStampLow  = ( 0xFFFFFFFF & iEvent.time ().value () ) ;
	myTreeVariables_.timeStampHigh = ( iEvent.time ().value () >> 32 ) ;
	myTreeVariables_.trgCut        = HLTCut ;

        dump3Ginfo(iEvent, iSetup, myTreeVariables_) ;
  
	dumpBarrelClusterInfo(iEvent, theGeometry, theCaloTopology,
			theBarrelEcalRecHits, 
			theBarrelBasicClusters, theBarrelSuperClusters, lazyTools, XtalMap, XtalMapCurved, myTreeVariables_) ;
	dumpEndcapClusterInfo(iEvent, theGeometry, theCaloTopology,
			theEndcapEcalRecHits, 
			theEndcapBasicClusters, theEndcapSuperClusters, lazyTools, XtalMap, XtalMapCurved, myTreeVariables_) ;

	dumpJetBarrelClusterInfo(iEvent, theGeometry, theCaloTopology, theBarrelEcalRecHits, 
			theBarrelBasicClusters, lazyTools, XtalMap, XtalMapCurved ) ;
	dumpJetEndcapClusterInfo(iEvent, theGeometry, theCaloTopology, theEndcapEcalRecHits, 
                        theEndcapBasicClusters, lazyTools, XtalMap, XtalMapCurved ) ;
  
	dumpVertexInfo(theRecVtxs, myTreeVariables_);

        tree_ -> Fill();
     }
  }
  delete lazyTools;

}






// -----------------------------------------------------------------------------------------

void EcalTimePhyTreeMaker::endJob ()
{
  using namespace std ;
  file_->cd () ;
  tree_->Write () ;
  file_->Close () ;

}


// -----------------------------------------------------------------------------------------

void EcalTimePhyTreeMaker::beginRun(edm::Run const &, edm::EventSetup const & eventSetup)
{
  // IC's
  eventSetup.get<EcalIntercalibConstantsRcd>().get(ical);
  // ADCtoGeV
  eventSetup.get<EcalADCToGeVConstantRcd>().get(agc);
  // transp corrections
  eventSetup.get<EcalLaserDbRecord>().get(laser);
}



// -----------------------------------------------------------------------------------------

std::string EcalTimePhyTreeMaker::intToString (int num)
{
    using namespace std ;
    ostringstream myStream ;
    myStream << num << flush ;
    return (myStream.str ()) ; //returns the string form of the stringstream object
}


// -----------------------------------------------------------------------------------------
/*
bool EcalTimePhyTreeMaker::HLTSelection( const edm::Event& iEvent ) {

   edm::Handle<edm::TriggerResults> triggers;
   iEvent.getByLabel( triggerSource_, triggers );
   
   const edm::TriggerNames& trgNames = iEvent.triggerNames( *triggers );

   int trgIndex  = trgNames.triggerIndex( triggerName_ );
   int trgResult = 0;
   if ( trgIndex == (int)(trgNames.size()) ) {
          std::cout<<" NO Matched Trigger -- Change Trigger !! "<<std::endl;
   } else {
          trgResult = triggers->accept(trgIndex);
   }

   bool pass =  ( trgResult == 1 ) ? true : false ;
   return pass ;
}
*/
int EcalTimePhyTreeMaker::HLTSelection( const edm::Event& iEvent, int cutVal, string str_head, string str_body ) {

   Handle<edm::TriggerResults> triggers;
   iEvent.getByLabel( triggerSource_, triggers );

   //cout<<" ** Trigger size = "<< triggers->size() <<endl;
   const edm::TriggerNames& trgNames = iEvent.triggerNames( *triggers );
   
   int cutSize = ( cutVal > 99 ) ? 3 : 2 ;
   int trgWordSize = str_head.size() + str_body.size() + cutSize + 3 ; //  3 -> _vX

   int trgResult = 0 ;
   for ( size_t i =0 ; i < trgNames.size(); i++ ) {
       string tName  = trgNames.triggerName( i );
       if  ( (int) tName.size() < trgWordSize ) continue ;
       string strhead = tName.substr(0, str_head.size() ) ;
       string strPt   = tName.substr( str_head.size() , cutSize ) ;
       int cutPt      =  atoi ( strPt.c_str() ) ;
       string strbody = tName.substr( str_head.size() + cutSize , str_body.size() ) ;
       //string strend = tName.substr( trgWordSize-3 ,3) ;
       if ( strhead == str_head && strbody == str_body && cutPt >= cutVal) {
          int trgIndex  = trgNames.triggerIndex(tName);
          int accept    = triggers->accept(trgIndex);
          if ( accept == 1 ) trgResult += cutPt ;
          //if ( accept == 0 ) trgResult = 0 ;
          //cout<<" -> "<< strhead << cutPt << strbody << strend <<" = "<< trgResult <<endl;
       }
   }
   return trgResult ;

}

int EcalTimePhyTreeMaker::HLTSelection( const edm::Event& iEvent ) {

   edm::Handle<edm::TriggerResults> triggers;
   iEvent.getByLabel( triggerSource_, triggers );
   
   const edm::TriggerNames& trgNames = iEvent.triggerNames( *triggers );

   int TrgValue = -1 ;
   int trgIndex1  = trgNames.triggerIndex( "HLT_Photon75_CaloIdVL_IsoL_v8" );
   int trgIndex2  = trgNames.triggerIndex( "HLT_Photon90_CaloIdVL_IsoL_v4" );
   int trgResult = 0;
   if ( trgIndex1 < (int)(trgNames.size()) ) {
          trgResult = triggers->accept(trgIndex1);
          if ( trgResult == 1 ) TrgValue = 75 ;
   }
   if ( trgIndex2 < (int)(trgNames.size()) ) {
          trgResult = triggers->accept(trgIndex2);
          if ( trgResult == 1 ) TrgValue = 90 ;
   } 
   if ( trgIndex1 < (int)(trgNames.size()) && trgIndex2 < (int)(trgNames.size()) ) {
          int trgResult1 = triggers->accept(trgIndex1);
          int trgResult2 = triggers->accept(trgIndex2);
          if ( trgResult1 == 1 && trgResult2 == 1 ) TrgValue =  0 ;
   } 
   if ( trgIndex1 == (int)(trgNames.size()) && trgIndex2 == (int)(trgNames.size()) ) {
          TrgValue = -99 ;
   }

   return TrgValue ;
}


bool EcalTimePhyTreeMaker::dumpEvtObjectInfo (const edm::Event& iEvent )
{
  // applies selection to physics objects and store them (for possible later usage)
  // with thos objects, applies event-based selection
  // if event passes selection, store objects into ntuple

   edm::Handle<reco::PhotonCollection> photons; 
   edm::Handle<reco::GsfElectronCollection> electrons; 
   edm::Handle<reco::MuonCollection> muons; 
   edm::Handle<reco::PFJetCollection> jets; 
   edm::Handle<reco::PFMETCollection> mets; 

   iEvent.getByLabel(PhotonSource_, photons);
   iEvent.getByLabel(ElectronSource_, electrons);
   iEvent.getByLabel(MuonSource_, muons);
   iEvent.getByLabel(JetSource_, jets);
   iEvent.getByLabel(METSource_, mets);
   
   // Ecal barrel RecHits 
   edm::Handle<EcalRecHitCollection> pBarrelEcalRecHits ;
   iEvent.getByLabel (barrelEcalRecHitCollection_, pBarrelEcalRecHits) ;
  
   // Ecal endcap RecHits
   edm::Handle<EcalRecHitCollection> pEndcapEcalRecHits ;
   iEvent.getByLabel (endcapEcalRecHitCollection_, pEndcapEcalRecHits) ;

  
   for(reco::PFMETCollection::const_iterator it = mets->begin(); it != mets->end(); it++) {
       if ( it->pt() < jetCuts_[4]  ) continue ;
       myTreeVariables_.met = it->et() ;
       myTreeVariables_.metPx = it->px() ;
       myTreeVariables_.metPy = it->py() ;
   }

   selectedJets.clear() ;
   for(reco::PFJetCollection::const_iterator it = jets->begin(); it != jets->end(); it++) {
       if ( it->pt() < jetCuts_[0] || fabs( it->eta() ) > jetCuts_[1] ) continue ;
       // Jet ID cuts
       /*
       if ( it->numberOfDaughters() < 2 )               continue ;
       if ( it->chargedEmEnergyFraction() >= 0.99 )     continue ;
       if ( it->neutralHadronEnergyFraction() >= 0.99 ) continue ;
       if ( it->neutralEmEnergyFraction() >= 0.99 )     continue ;
       if ( fabs( it->eta() ) < 2.4 && it->chargedHadronEnergyFraction() <=0 ) continue ;
       if ( fabs( it->eta() ) < 2.4 && it->chargedMultiplicity() <=0 ) continue ;
       */
       selectedJets.push_back( &(*it) ) ;
   }

   selectedElectrons.clear() ;
   float eidx = 11. ;
   for(reco::GsfElectronCollection::const_iterator it = electrons->begin(); it != electrons->end(); it++) {
       if ( it->pt() < electronCuts_[0] || fabs( it->eta() ) > electronCuts_[1] ) continue ;
       //double relIso =  ( it->chargedHadronIso()+ it->neutralHadronIso() + it->photonIso () ) / it->pt();
       // Isolation Cuts
       float ecalSumEt = ( it->isEB() ) ? max(0., it->dr03EcalRecHitSumEt() - 1. ) : it->dr03EcalRecHitSumEt();
       float hcalSumEt = it->dr03HcalTowerSumEt();
       float trkSumPt  = it->dr03TkSumPt();
       double relIso   = (ecalSumEt + hcalSumEt + trkSumPt) / it->pt() ;
       if ( relIso > electronCuts_[2] ) continue ;

       double dR = 999. ;
       for (size_t j=0; j < selectedJets.size(); j++ ) {
           double dR_ =  ROOT::Math::VectorUtil::DeltaR( it->p4(), selectedJets[j]->p4() ) ;
           if ( dR_ < dR ) dR = dR_ ;
       }
       if ( dR <= electronCuts_[3] ) continue ;

       double nLost = it->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() ;
       if ( nLost >= 2 ) continue ;
       eidx += 0.1 ;
       if ( !it->superCluster().isNull() ) sclist.push_back( make_pair(it->superCluster(), eidx ) );
       selectedElectrons.push_back( &(*it) ) ;
   }

   selectedMuons.clear() ;
   float midx = 13.0 ;
   for(reco::MuonCollection::const_iterator it = muons->begin(); it != muons->end(); it++) {
       if ( it->pt() < muonCuts_[0] || fabs( it->eta() ) > muonCuts_[1] ) continue ;
       // Isolation for PAT muon
       //double relIso =  ( it->chargedHadronIso()+ it->neutralHadronIso() + it->photonIso () ) / it->pt();
       // Isolation for RECO muon
       double relIso =0. ;
       if ( it->isIsolationValid() ) {
	 relIso = ( it->isolationR05().emEt + it->isolationR05().hadEt + it->isolationR05().sumPt ) / it->pt();
       }
       if ( relIso > muonCuts_[2] ) continue ;
       double dR = 999. ;
       for (size_t j=0; j < selectedJets.size(); j++ ) {
           double dR_ =  ROOT::Math::VectorUtil::DeltaR( it->p4(), selectedJets[j]->p4() ) ; 
           if ( dR_ < dR ) dR = dR_ ;
       }
       if ( dR <= muonCuts_[3] ) continue ;
       midx += 0.1 ;
       if ( !it->superCluster().isNull() ) sclist.push_back( make_pair(it->superCluster(), midx ) );
       selectedMuons.push_back( &(*it) ) ;
   }

   selectedPhotons.clear() ;
   float gidx = 22.0 ;
   for(reco::PhotonCollection::const_iterator it = photons->begin(); it != photons->end(); it++) {
       if ( it->pt() < photonCuts_[0] || fabs( it->eta() ) > photonCuts_[1] ) continue ;

       //float hcalIsoRatio = it->hcalTowerSumEtConeDR04() / it->pt() ;
       //if  ( ( hcalIsoRatio + it->hadronicOverEm() )*it->energy() >=  photonCuts_[2] ) continue ;

       // S_Minor Cuts from the seed cluster
       reco::CaloClusterPtr SCseed = it->superCluster()->seed() ;
       const EcalRecHitCollection* rechits = ( it->isEB()) ? pBarrelEcalRecHits.product () : pEndcapEcalRecHits.product() ;
       Cluster2ndMoments moments = EcalClusterTools::cluster2ndMoments(*SCseed, *rechits);
       float sMin =  moments.sMin  ;
       if ( sMin <= photonCuts_[5] || sMin >= photonCuts_[6] ) continue ;

       // Isolation cuts 
       /* 
       float ecalSumEt = it->ecalRecHitSumEtConeDR04();
       float hcalSumEt = it->hcalTowerSumEtConeDR04();
       float trkSumPt  = it->trkSumPtSolidConeDR04();
       
       bool trkIso  = ( ( trkSumPt / it->pt())     < photonIso_[0] ) ; 
       bool ecalIso = ( (ecalSumEt / it->energy()) < photonIso_[2] && ecalSumEt < photonIso_[1] ) ; 
       bool hcalIso = ( (hcalSumEt / it->energy()) < photonIso_[4] && hcalSumEt < photonIso_[3] ) ; 
       if ( !trkIso || !ecalIso || !hcalIso ) continue ;
       double dR = 999 ;
       for (size_t j=0; j < selectedJets.size(); j++ ) {
           double dR_ =  ROOT::Math::VectorUtil::DeltaR( it->p4(), selectedJets[j]->p4() ) ;
           if ( dR_ < dR ) dR = dR_ ;
       }
       if ( dR <= photonCuts_[3] ) continue ;
       */ 
       gidx += 0.1 ;
       if ( !it->superCluster().isNull() ) sclist.push_back( make_pair(it->superCluster(), gidx ) );
       selectedPhotons.push_back( &(*it) ) ;
   }

   int nJet   = static_cast<int> ( selectedJets.size() ) ;
   int nEle   = static_cast<int> ( selectedElectrons.size() ) ;
   int nPho   = static_cast<int> ( selectedPhotons.size() ) ;
   int nMuon  = static_cast<int> ( selectedMuons.size() ) ;

   // apply event-based selection
   bool pass = true ;
   if ( nPho< photonCuts_[4] ) pass = false ;
   if ( nJet < jetCuts_[2] )   pass = false ;

   // fill out the information for:
   //%%%%%%% Jets 
   if ( pass ) {
      for ( size_t k=0; k< selectedJets.size(); k++ ) {
          if ( k >= 10 ) break ;
	  myTreeVariables_.jetPx[k] = selectedJets[k]->p4().Px() ;
	  myTreeVariables_.jetPy[k] = selectedJets[k]->p4().Py() ;
	  myTreeVariables_.jetPz[k] = selectedJets[k]->p4().Pz() ;
	  myTreeVariables_.jetE[k]  = selectedJets[k]->p4().E()  ;
	  myTreeVariables_.jetNDau[k] = selectedJets[k]->numberOfDaughters() ;
	  myTreeVariables_.jetCM[k]   = selectedJets[k]->chargedMultiplicity() ;
	  myTreeVariables_.jetCEF[k]  = selectedJets[k]->chargedEmEnergyFraction() ;
	  myTreeVariables_.jetNHF[k]  = selectedJets[k]->neutralHadronEnergyFraction() ;
	  myTreeVariables_.jetNEF[k]  = selectedJets[k]->neutralEmEnergyFraction() ;
      }
      //cout<<" sc collection = "<< sclist.size() ; 
      //cout<<" N_e: "<< selectedElectrons.size() <<" N g: "<< selectedPhotons.size() ;
      //cout<<" N_m: "<< selectedMuons.size() <<" N_J:" << selectedJets.size() << endl;
      myTreeVariables_.nJets = nJet ;

      //%%%%%%% Electrons 
      for ( size_t k=0; k< selectedElectrons.size(); k++ ) {
          if ( k >= 10 ) break ;
	  myTreeVariables_.elePx[k] = selectedElectrons[k]->p4().Px() ;
	  myTreeVariables_.elePy[k] = selectedElectrons[k]->p4().Py() ;
	  myTreeVariables_.elePz[k] = selectedElectrons[k]->p4().Pz() ;
	  myTreeVariables_.eleE[k] = selectedElectrons[k]->p4().E() ;
          myTreeVariables_.eleEcalIso[k] = ( selectedElectrons[k]->isEB() ) ? 
                                        max(0., selectedElectrons[k]->dr03EcalRecHitSumEt() - 1. ) : selectedElectrons[k]->dr03EcalRecHitSumEt();
          myTreeVariables_.eleHcalIso[k] = selectedElectrons[k]->dr03HcalTowerSumEt() ;
          myTreeVariables_.eleTrkIso[k]  = selectedElectrons[k]->dr03TkSumPt() ;
          myTreeVariables_.eleNLostHits[k]  = selectedElectrons[k]->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() ;
      }
      myTreeVariables_.nElectrons = nEle ;

      //%%%%%%% Photons 
      for ( size_t k=0; k< selectedPhotons.size(); k++ ) {
          if ( k >= 10 ) break ;
	  myTreeVariables_.phoPx[k] = selectedPhotons[k]->p4().Px() ;
	  myTreeVariables_.phoPy[k] = selectedPhotons[k]->p4().Py() ;
	  myTreeVariables_.phoPz[k] = selectedPhotons[k]->p4().Pz() ;
	  myTreeVariables_.phoE[k] = selectedPhotons[k]->p4().E() ;
          myTreeVariables_.phoEcalIso[k] = selectedPhotons[k]->ecalRecHitSumEtConeDR04();
          myTreeVariables_.phoHcalIso[k] = selectedPhotons[k]->hcalTowerSumEtConeDR04();
          myTreeVariables_.phoTrkIso[k]  = selectedPhotons[k]->trkSumPtSolidConeDR04();
          myTreeVariables_.phoHovE[k]  = selectedPhotons[k]->hadronicOverEm();

          // S_Minor Cuts from the seed cluster
	  reco::CaloClusterPtr SCseed = selectedPhotons[k]->superCluster()->seed() ;
	  const EcalRecHitCollection* rechits = ( selectedPhotons[k]->isEB()) ? pBarrelEcalRecHits.product () : pEndcapEcalRecHits.product() ;
	  Cluster2ndMoments moments = EcalClusterTools::cluster2ndMoments(*SCseed, *rechits);
	  float sMin =  moments.sMin  ;
	  float sMaj =  moments.sMaj  ;
	  // seed Time 
	  pair<DetId, float> maxRH = EcalClusterTools::getMaximum( *SCseed, rechits );
	  DetId seedCrystalId = maxRH.first;
	  EcalRecHitCollection::const_iterator seedRH = rechits->find(seedCrystalId);
	  float seedTime = (float)seedRH->time();

          myTreeVariables_.phoSmin[k] = sMin ;
          myTreeVariables_.phoSmaj[k] = sMaj ;
          myTreeVariables_.phoTime[k] = seedTime ;

      }
      myTreeVariables_.nPhotons = nPho ;

      //%%%%%%% Muons 
      for ( size_t k=0; k< selectedMuons.size(); k++ ) {
          if ( k >= 10 ) break ;
	  myTreeVariables_.muPx[k] = selectedMuons[k]->p4().Px() ;
	  myTreeVariables_.muPy[k] = selectedMuons[k]->p4().Py() ;
	  myTreeVariables_.muPz[k] = selectedMuons[k]->p4().Pz() ;
	  myTreeVariables_.muE[k] = selectedMuons[k]->p4().E() ;
          myTreeVariables_.muEcalIso[k] = selectedMuons[k]->isolationR05().emEt ;
          myTreeVariables_.muHcalIso[k] = selectedMuons[k]->isolationR05().hadEt ;
          myTreeVariables_.muTrkIso[k]  = selectedMuons[k]->isolationR05().sumPt ;
      }
      myTreeVariables_.nMuons = nMuon ;

   }
   return pass ;

} 

void EcalTimePhyTreeMaker::dumpBarrelClusterInfo (const edm::Event& iEvent,
     				             const CaloGeometry * theGeometry,
					     const CaloTopology * theCaloTopology,
					     const EcalRecHitCollection* theBarrelEcalRecHits,
					     const reco::BasicClusterCollection* theBarrelBasicClusters,
					     const reco::SuperClusterCollection* theBarrelSuperClusters,
					     EcalClusterLazyTools* lazyTools,
					     const std::map<int,float> & XtalMap, //GFdoc unclear
					     const std::map<int,float> & XtalMapCurved, //GFdoc unclear
					     EcalTimePhyTreeContent & myTreeVarnumberOfSuperClustersiables_)
{            
  // get number of of objects already present in the tree (none if dumpBarrelClusterInfo is called first)
  // this will be the index looping over the SC collection;
  //numberOfSuperClusters = myTreeVariables_.nSuperClusters;
  //this will be the index looping over the BC collection;
  //numberOfClusters      = myTreeVariables_.nClusters;

  const EcalIntercalibConstantMap& icalMap = ical->getMap();
  float adcToGeV = float(agc->getEBValue());
  
  //number of superClusters in event (collection = vector!)
  myTreeVariables_.nBarrelSuperClusters = theBarrelSuperClusters->size () ;
  
  /////////////////////////////////////////////////////////////////////////////////////////
  //loop on all superclusters in event
  for (reco::SuperClusterCollection::const_iterator sclus = theBarrelSuperClusters->begin () ; 
       sclus != theBarrelSuperClusters->end () ;  ++sclus)  {//start Loop on SC's

       // matching SC from e, photon, mu
       bool getSC = false ;
       float objMatchId = -1 ;
       // sclist holds list of SC's matched to any of the selected objects
       for ( size_t i=0; i< sclist.size(); i++) {
           if ( sclist[i].first.isNull() )     continue; 
           if ( sclist[i].first.get()->eta() == sclus->eta() && sclist[i].first.get()->phi() == sclus->phi() ) {
              //cout<<" E = "<< sclist[i].first.get()->energy() <<" from "<< sclist[i].second <<endl ; 
              getSC = true ;
              objMatchId = sclist[i].second ;
              break ;
           } // loop over SC's matched to object
       }// loop over all SC's
       if ( !getSC ) continue ;

       // sMin or sMaj cuts
       reco::CaloClusterPtr SCseed = sclus->seed() ;
       Cluster2ndMoments moments = EcalClusterTools::cluster2ndMoments(*SCseed, *theBarrelEcalRecHits );
       float sMin =  moments.sMin  ;
       float sMaj =  moments.sMaj  ;
       myTreeVariables_.sMin[numberOfSuperClusters] = sMin ;
       myTreeVariables_.sMaj[numberOfSuperClusters] = sMaj ;

       //    int numberOfXtalsInSuperCluster = 0 ;//counter for all xtals in supercluster 
       myTreeVariables_.SCPIdx[numberOfSuperClusters] = objMatchId ;
       //cout<<" EB SC "<<numberOfSuperClusters<<"  objId: "<<objMatchId <<" ->"<< myTreeVariables_.SCPIdx[numberOfSuperClusters]<<endl;
       myTreeVariables_.nClustersInSuperCluster[numberOfSuperClusters] = sclus -> clustersSize () ;
       //myTreeVariables_.superClusterType[numberOfSuperClusters] = 0 ;
      
       myTreeVariables_.superClusterEta[numberOfSuperClusters] = sclus -> position ().eta () ;
       myTreeVariables_.superClusterPhi[numberOfSuperClusters] = sclus -> position ().phi () ;
       myTreeVariables_.superClusterX[numberOfSuperClusters] = sclus -> position ().x () ;
       myTreeVariables_.superClusterY[numberOfSuperClusters] = sclus -> position ().y () ;
       myTreeVariables_.superClusterZ[numberOfSuperClusters] = sclus -> position ().z () ;
       myTreeVariables_.superClusterRawEnergy[numberOfSuperClusters] = sclus -> rawEnergy () ;
       myTreeVariables_.superClusterPhiWidth[numberOfSuperClusters] = sclus -> phiWidth () ;
       myTreeVariables_.superClusterEtaWidth[numberOfSuperClusters] = sclus -> etaWidth () ;

       
       ///////////////////////////////////////////////////////////////////////////////////////
       // loop on barrel basic clusters from SC

      for ( reco::CaloCluster_iterator  clus =  sclus->clustersBegin()  ;
            clus != sclus->clustersEnd()  && numberOfClusters<MAXC;   ++clus) // loop on barrel Sclusters
          {

          double energy = (*clus)->energy () ;
	  double phi    = (*clus)->phi () ;
	  double eta    = (*clus)->eta () ;
	  if  ( energy/cosh(eta) <minEtEB_ ) continue;
	  double sinTheta         = fabs( sin( 2 *atan( exp(-1*(*clus)->eta()) ) ) );
	  double transverseEnergy = (*clus)->energy () * sinTheta;
	  double time = -1000.0 ; // gfdoc: work on this to provide a combination of crystals?
	  double ampli = 0. ;
	  double secondMin = 0. ;
	  double secondTime = -1000. ;
	  int numberOfXtalsInCluster=0 ;//counter for all xtals in cluster
	  int numXtalsinCluster = 0 ;   //xtals in cluster above 3sigma noise  
	  EBDetId maxDet ;
          EBDetId secDet ;
       
          // GFdoc clusterDetIds holds crystals that participate to this basic cluster 
          std::vector<std::pair<DetId, float> > clusterDetIds = (*clus)->hitsAndFractions() ; //get these from the cluster
       
          //////////////////////////////////////////////////////////////////////////////////////
          //loop on xtals in cluster
          for (std::vector<std::pair<DetId, float> >::const_iterator detitr = clusterDetIds.begin () ; 
	      detitr != clusterDetIds.end ()  && numberOfXtalsInCluster<MAXXTALINC; // && numberOfXtals<MAXXTAL ; 
	      ++detitr)// loop on rechics of barrel basic clusters
	   {
	      // Here I use the "find" on a recHit collection... I have been warned...   (GFdoc: ??)
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
	     //The checking above should no longer be needed... 
	     //as only those in the cluster would already have rechits..
	   
	     // GFdoc this is one crystal in the basic cluster
	     EcalRecHit myhit = (*thishit) ;
	   
             // SIC Feb 14 2011 -- Add check on RecHit flags (takes care of spike cleaning in 42X)
             if( !( myhit.checkFlag(EcalRecHit::kGood) ||
                    myhit.checkFlag(EcalRecHit::kOutOfTime) ||
                    myhit.checkFlag(EcalRecHit::kPoorCalib) 
		    ) 
		 )
                 continue;

             // thisamp is the EB amplitude of the current rechit
	     double thisamp  = myhit.energy () ;
	     double thisChi2 = myhit.chi2 ();
	     double thisOutOfTimeChi2 = myhit.outOfTimeChi2 ();
	   
	     EcalIntercalibConstantMap::const_iterator icalit = icalMap.find(detitr->first);
	     EcalIntercalibConstant icalconst = 1;
	     if( icalit!=icalMap.end() ) {
	       icalconst = (*icalit);
	     } else {
	       edm::LogError("EcalTimePhyTreeMaker") << "No intercalib const found for xtal "
						<< (detitr->first).rawId();
   	     }
	   
	     // get laser coefficient
	     float lasercalib = 1.;
	     lasercalib = laser->getLaserCorrection( detitr->first, iEvent.time());
	     // std::cout << "GF EB debug: " << lasercalib << std::endl;

	     // discard rechits with A/sigma < 12
	     if ( thisamp/(icalconst*lasercalib*adcToGeV) < (1.1*12) ) continue;

	     if (thisamp > 0.027) //cut on energy->number of crystals in cluster above 3sigma noise; gf: desirable?
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
	     
	     GlobalPoint pos = theGeometry->getPosition((myhit).detid());
	     myTreeVariables_.xtalInBCEta[numberOfClusters][numberOfXtalsInCluster]=      pos.eta();
	     myTreeVariables_.xtalInBCPhi[numberOfClusters][numberOfXtalsInCluster]=      pos.phi();
	     double thistime = myhit.time();
	     if(doTimeVSAmpliCorrection_)  {
	       thistime += theTimeCorrector_.getCorrection((float) thisamp/(icalconst*lasercalib*adcToGeV), pos.eta()  );
	       //std::cout << "EB: " << thisamp << "\t" << ((float) thisamp/(icalconst*lasercalib*adcToGeV)) << "\t" << pos.eta() 
	       //<< "\t\t" << theTimeCorrector_.getCorrection((float) thisamp/(icalconst*lasercalib*adcToGeV), pos.eta()  ) 
	       //<< std::endl;
	     } else { /* std::cout << "not doing correction " << std::endl;*/ }
	   
	    if(myhit.isTimeErrorValid())
              myTreeVariables_.xtalInBCTimeErr[numberOfClusters][numberOfXtalsInCluster]= myhit.timeError();
            else
              myTreeVariables_.xtalInBCTimeErr[numberOfClusters][numberOfXtalsInCluster]= -999999;

	    // xtal variables inside a barrel basic cluster 
	    myTreeVariables_.xtalInBCEnergy[numberOfClusters][numberOfXtalsInCluster]=       (float) thisamp;
	    myTreeVariables_.xtalInBCTime[numberOfClusters][numberOfXtalsInCluster]=         (float) thistime; 
	    myTreeVariables_.xtalInBCHashedIndex[numberOfClusters][numberOfXtalsInCluster]=  EBDetId (detitr -> first).hashedIndex () ; 
	    myTreeVariables_.xtalInBCIEta[numberOfClusters][numberOfXtalsInCluster]=         EBDetId((detitr -> first)).ieta();
	    myTreeVariables_.xtalInBCIPhi[numberOfClusters][numberOfXtalsInCluster]=         EBDetId((detitr -> first)).iphi();
	    myTreeVariables_.xtalInBCIx[numberOfClusters][numberOfXtalsInCluster]=           -999999; 
	    myTreeVariables_.xtalInBCIy[numberOfClusters][numberOfXtalsInCluster]=           -999999; 
	    myTreeVariables_.xtalInBCFlag[numberOfClusters][numberOfXtalsInCluster]=         myhit.recoFlag(); 
	    myTreeVariables_.xtalInBCAmplitudeADC[numberOfClusters][numberOfXtalsInCluster]= (float) thisamp/(icalconst*lasercalib*adcToGeV);
	    myTreeVariables_.xtalInBCChi2[numberOfClusters][numberOfXtalsInCluster]=         thisChi2;
	    myTreeVariables_.xtalInBCOutOfTimeChi2[numberOfClusters][numberOfXtalsInCluster]=thisOutOfTimeChi2;
	    // note: SwissCross = 1 - E4/E1   
            myTreeVariables_.xtalInBCSwissCross[numberOfClusters][numberOfXtalsInCluster] =
             EcalTools::swissCross(detitr->first,*theBarrelEcalRecHits,0.5);

	   
	    numberOfXtalsInCluster++ ; // increment number of crystals in basic cluster
	   
	 } //end loop on rechits within barrel basic clusters
       //////////////////////////////////////////////////////
       
       float E2 = (float)ampli + (float)secondMin ;
       
       //Set some more values for this cluster
       myTreeVariables_.clusterMom[numberOfClusters] = numberOfSuperClusters;
       myTreeVariables_.CPIdx[numberOfClusters]      = objMatchId;
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
       if(FILL_SHAPE_VARS) {
	 myTreeVariables_.clusterE2x2[numberOfClusters] = lazyTools -> e2x2(*(*clus));
	 myTreeVariables_.clusterE3x2[numberOfClusters] = lazyTools -> e3x2(*(*clus));
	 myTreeVariables_.clusterE3x3[numberOfClusters] = lazyTools -> e3x3(*(*clus));
	 myTreeVariables_.clusterE4x4[numberOfClusters] = lazyTools -> e4x4(*(*clus));
	 myTreeVariables_.clusterE5x5[numberOfClusters] = lazyTools -> e5x5(*(*clus));
	 
	 myTreeVariables_.clusterE2x5Right[numberOfClusters]  = lazyTools -> e2x5Right(*(*clus));
	 myTreeVariables_.clusterE2x5Left[numberOfClusters]   = lazyTools -> e2x5Left(*(*clus));
	 myTreeVariables_.clusterE2x5Top[numberOfClusters]    = lazyTools -> e2x5Top(*(*clus));
	 myTreeVariables_.clusterE2x5Bottom[numberOfClusters] = lazyTools -> e2x5Bottom(*(*clus));
	 
	 myTreeVariables_.clusterCovEtaEta[numberOfClusters] = (lazyTools -> covariances(*(*clus)))[0];
	 myTreeVariables_.clusterCovPhiPhi[numberOfClusters] = (lazyTools -> covariances(*(*clus)))[2];
	 myTreeVariables_.clusterCovEtaPhi[numberOfClusters] = (lazyTools -> covariances(*(*clus)))[1];
	 
	 myTreeVariables_.clusterLat[numberOfClusters] = (lazyTools -> lat(*(*clus)))[2];
	 myTreeVariables_.clusterEtaLat[numberOfClusters] = (lazyTools -> lat(*(*clus)))[0];
	 myTreeVariables_.clusterPhiLat[numberOfClusters] = (lazyTools -> lat(*(*clus)))[1];
	 
	 myTreeVariables_.clusterZernike20[numberOfClusters] = lazyTools -> zernike20(*(*clus));
	 myTreeVariables_.clusterZernike42[numberOfClusters] = lazyTools -> zernike42(*(*clus));
       }
       
       numberOfClusters++ ;
     } //end (BASIC)cluster loop
     numberOfSuperClusters++ ;
    } //end supercluster loop
   
   
   //number of Bclusters already stored in the ntuple
   myTreeVariables_.nSuperClusters = numberOfSuperClusters;
   myTreeVariables_.nClusters = numberOfClusters;
   //myTreeVariables_.nXtals    = numberOfXtals ; - obsolete
   return ;
} // end dumpBarrelClusterInfo  


// -----------------------------------------------------------------------------------------

void EcalTimePhyTreeMaker::dumpEndcapClusterInfo (const edm::Event& iEvent,
                                             const CaloGeometry * theGeometry,
					     const CaloTopology * theCaloTopology,
					     const EcalRecHitCollection* theEndcapEcalRecHits,
					     const reco::BasicClusterCollection* theEndcapBasicClusters,
					     const reco::SuperClusterCollection* theEndcapSuperClusters,
					     EcalClusterLazyTools* lazyTools,
					     const std::map<int,float> & XtalMap,
					     const std::map<int,float> & XtalMapCurved,
					     EcalTimePhyTreeContent & myTreeVariables_)
{
  // counters come from the ntuple are to account for what was added in dumpBarrelClusterInf
  // numberOfSuperClusters = myTreeVariables_.nSuperClusters;
  // numberOfClusters      = myTreeVariables_.nClusters;
  // int numberOfXtals         = myTreeVariables_.nXtals ; // this is number of crystals associated to any cluster

  const EcalIntercalibConstantMap& icalMap = ical->getMap();
  float adcToGeV = float(agc->getEEValue());
  
  
  //number of superClusters in event (collection = vector!)
  myTreeVariables_.nEndcapSuperClusters = theEndcapSuperClusters->size () ;
  
  //loop on all endcap superclusters in event
  for (reco::SuperClusterCollection::const_iterator sclus = theEndcapSuperClusters->begin () ; 
       sclus != theEndcapSuperClusters->end () ;  ++sclus) 
    {//loop on SC's
      
       // matching SC from objects
       bool getSC = false ;
       float objMatchId = -1 ;
       for ( size_t i=0; i< sclist.size(); i++) {
           if ( sclist[i].first.isNull() ) continue; 
           if ( sclist[i].first.get()->eta() == sclus->eta() && sclist[i].first.get()->phi() == sclus->phi() ) {
              //cout<<" E = "<< sclist[i].first.get()->energy() <<" from "<< sclist[i].second <<endl ; 
              getSC = true ;
              objMatchId = sclist[i].second ;
              break ;
           }
       }

      if ( !getSC ) continue ;

      // sMin or sMaj cuts
      reco::CaloClusterPtr SCseed = sclus->seed() ;
      Cluster2ndMoments moments = EcalClusterTools::cluster2ndMoments(*SCseed, *theEndcapEcalRecHits );
      float sMin =  moments.sMin  ;
      float sMaj =  moments.sMaj  ;
      myTreeVariables_.sMin[numberOfSuperClusters] = sMin ;
      myTreeVariables_.sMaj[numberOfSuperClusters] = sMaj ;

      myTreeVariables_.SCPIdx[numberOfSuperClusters] = objMatchId ;
      //cout<<" EE SC "<<numberOfSuperClusters<<"  objId: "<<objMatchId <<" ->"<< myTreeVariables_.SCPIdx[numberOfSuperClusters]<<endl;
      myTreeVariables_.nClustersInSuperCluster[numberOfSuperClusters] = sclus -> clustersSize () ;
     
      myTreeVariables_.superClusterEta[numberOfSuperClusters] = sclus -> position ().eta () ;
      myTreeVariables_.superClusterPhi[numberOfSuperClusters] = sclus -> position ().phi () ;
      myTreeVariables_.superClusterX[numberOfSuperClusters] = sclus -> position ().x () ;
      myTreeVariables_.superClusterY[numberOfSuperClusters] = sclus -> position ().y () ;
      myTreeVariables_.superClusterZ[numberOfSuperClusters] = sclus -> position ().z () ;
      myTreeVariables_.superClusterRawEnergy[numberOfSuperClusters] = sclus -> energy () ;
      myTreeVariables_.superClusterPhiWidth[numberOfSuperClusters] = sclus -> phiWidth () ;
      myTreeVariables_.superClusterEtaWidth[numberOfSuperClusters] = sclus -> etaWidth () ;

  

  ///////////////////////////////////////////////////////////////////////////////////////
   for ( reco::CaloCluster_iterator  clus =  sclus->clustersBegin()  ;
         clus != sclus->clustersEnd()  && numberOfClusters<MAXC;   ++clus) // loop on barrel Sclusters
       {

         double energy = (*clus)->energy () ;
         double phi    = (*clus)->phi () ;
         double eta    = (*clus)->eta () ;
	 if  ( energy/cosh(eta) <minEtEE_ ) continue;
	 double sinTheta         = fabs( sin( 2 *atan( exp(-1*(*clus)->eta()) ) ) );
	 double transverseEnergy = (*clus)->energy () * sinTheta;
         double time = -1000.0 ;  // gfdoc: work on this to provide a combination of crystals?
         double ampli = 0. ;
         double secondMin = 0. ;
         double secondTime = -1000. ;
         int numberOfXtalsInCluster=0 ;//counter for all xtals in cluster
         int numXtalsinCluster = 0 ;//xtals in cluster above 3sigma noise  
         EEDetId maxDet ;
         EEDetId secDet ;
         
         std::vector<std::pair<DetId, float> > clusterDetIds = (*clus)->hitsAndFractions() ; //get these from the cluster

        //loop on xtals in cluster
         for (std::vector<std::pair<DetId, float> >::const_iterator detitr = clusterDetIds.begin () ; 
              detitr != clusterDetIds.end () && numberOfXtalsInCluster<MAXXTALINC;// && numberOfXtals<MAXXTAL ; 
              ++detitr)// loop on rechits of endcap basic clusters
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
             
             // SIC Feb 14 2011 -- Add check on RecHit flags (takes care of spike cleaning in 42X)
             if( !( myhit.checkFlag(EcalRecHit::kGood) ||
                    myhit.checkFlag(EcalRecHit::kOutOfTime) ||
                    myhit.checkFlag(EcalRecHit::kPoorCalib) 
		    ) 
		 )
                 continue;

	     // thisamp is the EE amplitude of the current rechit
	     double thisamp  = myhit.energy () ;
	     double thisChi2 = myhit.chi2 ();
	     double thisOutOfTimeChi2 = myhit.outOfTimeChi2 ();

	     // get laser coefficient
	     float lasercalib = 1.;
	     lasercalib = laser->getLaserCorrection( detitr->first, iEvent.time());
	     //std::cout << "GF EE debug: " << lasercalib << std::endl;
	     
             EcalIntercalibConstantMap::const_iterator icalit = icalMap.find(detitr->first);
             EcalIntercalibConstant icalconst = 1;
             if( icalit!=icalMap.end() ) {
               icalconst = (*icalit);
             } else {
               edm::LogError("EcalTimePhyTreeMaker") << "No intercalib const found for xtal "
						  << (detitr->first).rawId();
             }
	     // don't store rechits with A/sigma < 12
	     if ( thisamp/(icalconst*lasercalib*adcToGeV) < (2.2*12) ) continue;
	     
	     
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

	     
	     
	     GlobalPoint pos = theGeometry->getPosition((myhit).detid());
	     myTreeVariables_.xtalInBCEta[numberOfClusters][numberOfXtalsInCluster]=      pos.eta();
	     myTreeVariables_.xtalInBCPhi[numberOfClusters][numberOfXtalsInCluster]=      pos.phi();
	     
	     double thistime = myhit.time ();
	     if(doTimeVSAmpliCorrection_)  {
	       thistime += theTimeCorrector_.getCorrection((float) thisamp/(icalconst*lasercalib*adcToGeV), pos.eta()  );
	       //std::cout << "EE: " << thisamp << "\t" << ((float) thisamp/(icalconst*lasercalib*adcToGeV)) << "\t" << pos.eta() 
	       // << "\t\t" << theTimeCorrector_.getCorrection((float) thisamp/(icalconst*lasercalib*adcToGeV), pos.eta()  ) 
	       //<< std::endl;
	     } else { /*std::cout << "not doing correction " << std::endl; */}

	     // xtal variables inside an endcap basic cluster 
	      myTreeVariables_.xtalInBCEnergy[numberOfClusters][numberOfXtalsInCluster]=      (float) thisamp;
	      myTreeVariables_.xtalInBCTime[numberOfClusters][numberOfXtalsInCluster]=        (float) thistime;
	      myTreeVariables_.xtalInBCHashedIndex[numberOfClusters][numberOfXtalsInCluster]= EEDetId (detitr -> first).hashedIndex () ; 
	      myTreeVariables_.xtalInBCIEta[numberOfClusters][numberOfXtalsInCluster]=-999999;
	      myTreeVariables_.xtalInBCIPhi[numberOfClusters][numberOfXtalsInCluster]=-999999;
	      myTreeVariables_.xtalInBCIx[numberOfClusters][numberOfXtalsInCluster]=          EEDetId((detitr -> first)).ix();
	      myTreeVariables_.xtalInBCIy[numberOfClusters][numberOfXtalsInCluster]=          EEDetId((detitr -> first)).iy();
	      myTreeVariables_.xtalInBCFlag[numberOfClusters][numberOfXtalsInCluster]=         myhit.recoFlag(); 
	      myTreeVariables_.xtalInBCAmplitudeADC[numberOfClusters][numberOfXtalsInCluster]=      (float) thisamp/(icalconst*lasercalib*adcToGeV);
              myTreeVariables_.xtalInBCChi2[numberOfClusters][numberOfXtalsInCluster]=         thisChi2;
              myTreeVariables_.xtalInBCOutOfTimeChi2[numberOfClusters][numberOfXtalsInCluster]=thisOutOfTimeChi2;
              myTreeVariables_.xtalInBCSwissCross[numberOfClusters][numberOfXtalsInCluster] =
                EcalTools::swissCross(detitr->first,*theEndcapEcalRecHits,0.5);

              numberOfXtalsInCluster++ ; // increment number of crystals in basic cluster
	     
           } //end loop on rechics within endcap basic clusters
         //////////////////////////////////////////////////////
	 
         float E2 = (float)ampli + (float)secondMin ;
         
         //Set some more values for this cluster
         myTreeVariables_.clusterMom[numberOfClusters] = numberOfSuperClusters;
         myTreeVariables_.CPIdx[numberOfClusters]      = objMatchId;
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

	 // (basic) cluster shapes for endcap
         myTreeVariables_.clusterE2x2[numberOfClusters] = lazyTools -> e2x2(*(*clus));

         myTreeVariables_.clusterE3x2[numberOfClusters] = lazyTools -> e3x2(*(*clus));
         myTreeVariables_.clusterE3x3[numberOfClusters] = lazyTools -> e3x3(*(*clus));
         myTreeVariables_.clusterE4x4[numberOfClusters] = lazyTools -> e4x4(*(*clus));
         myTreeVariables_.clusterE5x5[numberOfClusters] = lazyTools -> e5x5(*(*clus));

         myTreeVariables_.clusterE2x5Right[numberOfClusters] = lazyTools -> e2x5Right(*(*clus));
         myTreeVariables_.clusterE2x5Left[numberOfClusters] = lazyTools -> e2x5Left(*(*clus));
         myTreeVariables_.clusterE2x5Top[numberOfClusters] = lazyTools -> e2x5Top(*(*clus));
         myTreeVariables_.clusterE2x5Bottom[numberOfClusters] = lazyTools -> e2x5Bottom(*(*clus));

         myTreeVariables_.clusterCovEtaEta[numberOfClusters] = (lazyTools -> covariances(*(*clus)))[0];
         myTreeVariables_.clusterCovPhiPhi[numberOfClusters] = (lazyTools -> covariances(*(*clus)))[2];
         myTreeVariables_.clusterCovEtaPhi[numberOfClusters] = (lazyTools -> covariances(*(*clus)))[1];

         myTreeVariables_.clusterLat[numberOfClusters] = (lazyTools -> lat(*(*clus)))[2];
         myTreeVariables_.clusterEtaLat[numberOfClusters] = (lazyTools -> lat(*(*clus)))[0];
         myTreeVariables_.clusterPhiLat[numberOfClusters] = (lazyTools -> lat(*(*clus)))[1];

         myTreeVariables_.clusterZernike20[numberOfClusters] = lazyTools -> zernike20(*(*clus));
         myTreeVariables_.clusterZernike42[numberOfClusters] = lazyTools -> zernike42(*(*clus));
    
         numberOfClusters++ ;
    }//end endcap basic cluster loop

    numberOfSuperClusters++ ;
   } //end endcap supercluster loopsuperClusterEnergySum
     
  //number of B clusters already stored in the ntuple
  myTreeVariables_.nSuperClusters  = numberOfSuperClusters ;
  myTreeVariables_.nClusters  = numberOfClusters ;

  return ;
} // end dumpEndcapClusterInfo  


void EcalTimePhyTreeMaker::dumpJetBarrelClusterInfo (const edm::Event& iEvent,
     				             const CaloGeometry * theGeometry,
					     const CaloTopology * theCaloTopology,
					     const EcalRecHitCollection* theBarrelEcalRecHits,
					     const reco::BasicClusterCollection* theBarrelBasicClusters,
					     EcalClusterLazyTools* lazyTools,
					     const std::map<int,float> & XtalMap, 
					     const std::map<int,float> & XtalMapCurved  )
{            
  // get number of of objects already present in the tree (none if dumpBarrelClusterInfo is called first)
  // this will be the index looping over the SC collection;
  //numberOfSuperClusters = myTreeVariables_.nSuperClusters;
  //numberOfClusters      = myTreeVariables_.nClusters;

  const EcalIntercalibConstantMap& icalMap = ical->getMap();
  float adcToGeV = float(agc->getEBValue());
  
  ///////////////////////////////////////////////////////////////////////////////////////
  // loop on barrel basic clusters 
   bool getC = false;
   for (reco::BasicClusterCollection::const_iterator clus = theBarrelBasicClusters->begin () ;
        clus != theBarrelBasicClusters->end ()  && numberOfClusters<MAXC;  ++clus) {   

          double energy = clus->energy () ;
	  double phi    = clus->phi () ;
	  double eta    = clus->eta () ;
	  if  ( energy/cosh(eta) <minEtEB_ ) continue;
	  double sinTheta         = fabs( sin( 2 *atan( exp(-1*clus->eta()) ) ) );
	  double transverseEnergy = clus->energy () * sinTheta;
	  double time = -1000.0 ; // gfdoc: work on this to provide a combination of crystals?
	  double ampli = 0. ;
	  double secondMin = 0. ;
	  double secondTime = -1000. ;
	  int numberOfXtalsInCluster=0 ;//counter for all xtals in cluster
	  int numXtalsinCluster = 0 ;   //xtals in cluster above 3sigma noise  
	  EBDetId maxDet ;
          EBDetId secDet ;
       
          // matching SC from objects
          float jidx = 100. ;
          getC = false;
          for ( size_t i=0; i< selectedJets.size(); i++) {
               bool matched = matching( clus->energy(), clus->position(), selectedJets[i]->vertex(), selectedJets[i]->p4() ) ;
               jidx += 0.1 ;
               if ( matched ) { 
                  getC = true ;
                  //cout <<" ** from jet "<< i <<" / "<< selectedJets.size();
                  break ;
               }
          }
          if ( !getC ) continue ;


          //number of superClusters in event (collection = vector!)
          myTreeVariables_.SCPIdx[numberOfSuperClusters] = jidx ;
          //cout<<" BJ Cl "<<numberOfClusters<<"  objId: "<< jidx <<" ->"<< myTreeVariables_.SCPIdx[numberOfSuperClusters]<<endl;
	  myTreeVariables_.superClusterEta[numberOfSuperClusters] = -99999 ;
	  myTreeVariables_.superClusterPhi[numberOfSuperClusters] = -99999 ;
	  myTreeVariables_.superClusterX[numberOfSuperClusters] = -99999 ;
	  myTreeVariables_.superClusterY[numberOfSuperClusters] = -99999 ;
	  myTreeVariables_.superClusterZ[numberOfSuperClusters] = -99999 ;
	  myTreeVariables_.superClusterRawEnergy[numberOfSuperClusters] = -999999 ;
	  myTreeVariables_.superClusterPhiWidth[numberOfSuperClusters]  = -999999 ;
	  myTreeVariables_.superClusterEtaWidth[numberOfSuperClusters]  = -999999 ; 
          myTreeVariables_.sMin[numberOfSuperClusters] = -1 ;
          myTreeVariables_.sMaj[numberOfSuperClusters] = -1 ;

          // GFdoc clusterDetIds holds crystals that participate to this basic cluster 
          std::vector<std::pair<DetId, float> > clusterDetIds = clus->hitsAndFractions() ; //get these from the cluster
       
          //////////////////////////////////////////////////////////////////////////////////////
          //loop on xtals in cluster
          for (std::vector<std::pair<DetId, float> >::const_iterator detitr = clusterDetIds.begin () ; 
	      detitr != clusterDetIds.end ()  && numberOfXtalsInCluster<MAXXTALINC; // && numberOfXtals<MAXXTAL ; 
	      ++detitr)// loop on rechics of barrel basic clusters
	   {

	      // Here I use the "find" on a digi collection... I have been warned...   (GFdoc: ??)
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
	     //The checking above should no longer be needed... 
	     //as only those in the cluster would already have rechits..
	   
	     // GFdoc this is one crystal in the basic cluster
	     EcalRecHit myhit = (*thishit) ;
	   
             // SIC Feb 14 2011 -- Add check on RecHit flags (takes care of spike cleaning in 42X)
             if( !( myhit.checkFlag(EcalRecHit::kGood) ||
                    myhit.checkFlag(EcalRecHit::kOutOfTime) ||
                    myhit.checkFlag(EcalRecHit::kPoorCalib) 
		    ) 
		 )
                 continue;
	   
	     // thisamp is the EB amplitude of the current rechit
	     double thisamp  = myhit.energy () ;
	     double thistime = myhit.time ();
	     double thisChi2 = myhit.chi2 ();
	     double thisOutOfTimeChi2 = myhit.outOfTimeChi2 ();

	   
	     EcalIntercalibConstantMap::const_iterator icalit = icalMap.find(detitr->first);
	     EcalIntercalibConstant icalconst = 1;
	     if( icalit!=icalMap.end() ) {
	       icalconst = (*icalit);
	     } else {
	       edm::LogError("EcalTimePhyTreeMaker") << "No intercalib const found for xtal "
						<< (detitr->first).rawId();
   	     }
	   
	     // get laser coefficient
	     float lasercalib = 1.;
	     lasercalib = laser->getLaserCorrection( detitr->first, iEvent.time());
	     // std::cout << "GF EB debug: " << lasercalib << std::endl;

	     // discard rechits with A/sigma < 12
	     if ( thisamp/(icalconst*lasercalib*adcToGeV) < (1.1*12) ) continue;

	     if (thisamp > 0.027) //cut on energy->number of crystals in cluster above 3sigma noise; gf: desirable?
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

	   
	    if(myhit.isTimeErrorValid())
              myTreeVariables_.xtalInBCTimeErr[numberOfClusters][numberOfXtalsInCluster]= myhit.timeError();
            else
              myTreeVariables_.xtalInBCTimeErr[numberOfClusters][numberOfXtalsInCluster]= -999999;

	    // xtal variables inside a barrel basic cluster 
	    myTreeVariables_.xtalInBCEnergy[numberOfClusters][numberOfXtalsInCluster]=       (float) thisamp;
	    myTreeVariables_.xtalInBCTime[numberOfClusters][numberOfXtalsInCluster]=         (float) thistime; 
	    myTreeVariables_.xtalInBCHashedIndex[numberOfClusters][numberOfXtalsInCluster]=  EBDetId (detitr -> first).hashedIndex () ; 
	    myTreeVariables_.xtalInBCIEta[numberOfClusters][numberOfXtalsInCluster]=         EBDetId((detitr -> first)).ieta();
	    myTreeVariables_.xtalInBCIPhi[numberOfClusters][numberOfXtalsInCluster]=         EBDetId((detitr -> first)).iphi();
	    myTreeVariables_.xtalInBCIx[numberOfClusters][numberOfXtalsInCluster]=           -999999; 
	    myTreeVariables_.xtalInBCIy[numberOfClusters][numberOfXtalsInCluster]=           -999999; 
	    myTreeVariables_.xtalInBCFlag[numberOfClusters][numberOfXtalsInCluster]=         myhit.recoFlag(); 
	    myTreeVariables_.xtalInBCAmplitudeADC[numberOfClusters][numberOfXtalsInCluster]= (float) thisamp/(icalconst*lasercalib*adcToGeV);
	    myTreeVariables_.xtalInBCChi2[numberOfClusters][numberOfXtalsInCluster]=         thisChi2;
	    myTreeVariables_.xtalInBCOutOfTimeChi2[numberOfClusters][numberOfXtalsInCluster]=thisOutOfTimeChi2;
	    // note: SwissCross = 1 - E4/E1   
            myTreeVariables_.xtalInBCSwissCross[numberOfClusters][numberOfXtalsInCluster] =
             EcalTools::swissCross(detitr->first,*theBarrelEcalRecHits,0.5);

  
	    GlobalPoint pos = theGeometry->getPosition((myhit).detid());
	    myTreeVariables_.xtalInBCEta[numberOfClusters][numberOfXtalsInCluster]=      pos.eta();
	    myTreeVariables_.xtalInBCPhi[numberOfClusters][numberOfXtalsInCluster]=      pos.phi();
	   
	    numberOfXtalsInCluster++ ; // increment number of crystals in basic cluster
	   
	 } //end loop on rechits within barrel basic clusters
       //////////////////////////////////////////////////////
       
       float E2 = (float)ampli + (float)secondMin ;
       
       //Set some more values for this cluster
       myTreeVariables_.clusterMom[numberOfClusters] = numberOfSuperClusters ;
       myTreeVariables_.CPIdx[numberOfClusters] = jidx ;
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
       if(FILL_SHAPE_VARS) {
	 myTreeVariables_.clusterE2x2[numberOfClusters] = lazyTools -> e2x2(*clus);
	 myTreeVariables_.clusterE3x2[numberOfClusters] = lazyTools -> e3x2(*clus);
	 myTreeVariables_.clusterE3x3[numberOfClusters] = lazyTools -> e3x3(*clus);
	 myTreeVariables_.clusterE4x4[numberOfClusters] = lazyTools -> e4x4(*clus);
	 myTreeVariables_.clusterE5x5[numberOfClusters] = lazyTools -> e5x5(*clus);
	 
	 myTreeVariables_.clusterE2x5Right[numberOfClusters]  = lazyTools -> e2x5Right(*clus);
	 myTreeVariables_.clusterE2x5Left[numberOfClusters]   = lazyTools -> e2x5Left(*clus);
	 myTreeVariables_.clusterE2x5Top[numberOfClusters]    = lazyTools -> e2x5Top(*clus);
	 myTreeVariables_.clusterE2x5Bottom[numberOfClusters] = lazyTools -> e2x5Bottom(*clus);
	 
	 myTreeVariables_.clusterCovEtaEta[numberOfClusters] = (lazyTools -> covariances(*clus))[0];
	 myTreeVariables_.clusterCovPhiPhi[numberOfClusters] = (lazyTools -> covariances(*clus))[2];
	 myTreeVariables_.clusterCovEtaPhi[numberOfClusters] = (lazyTools -> covariances(*clus))[1];
	 
	 myTreeVariables_.clusterLat[numberOfClusters] = (lazyTools -> lat( *clus ))[2];
	 myTreeVariables_.clusterEtaLat[numberOfClusters] = (lazyTools -> lat( *clus ))[0];
	 myTreeVariables_.clusterPhiLat[numberOfClusters] = (lazyTools -> lat( *clus ))[1];
	 
	 myTreeVariables_.clusterZernike20[numberOfClusters] = lazyTools -> zernike20(*clus);
	 myTreeVariables_.clusterZernike42[numberOfClusters] = lazyTools -> zernike42(*clus);
       }
       numberOfClusters++ ;
       
     } //end (BASIC)cluster loop

     if ( getC ) {
        //cout<<"    SC "<< numberOfSuperClusters <<endl; 
        myTreeVariables_.nClustersInSuperCluster[numberOfSuperClusters] = numberOfClusters - myTreeVariables_.nClusters ;
        numberOfSuperClusters++ ;      
        myTreeVariables_.nSuperClusters = numberOfSuperClusters;
        myTreeVariables_.nClusters = numberOfClusters;
     }
     return ;
} // end dumpBarrelClusterInfo  


void EcalTimePhyTreeMaker::dumpJetEndcapClusterInfo (const edm::Event& iEvent,
     				             const CaloGeometry * theGeometry,
					     const CaloTopology * theCaloTopology,
                                             const EcalRecHitCollection* theEndcapEcalRecHits,
                                             const reco::BasicClusterCollection* theEndcapBasicClusters,
					     EcalClusterLazyTools* lazyTools,
					     const std::map<int,float> & XtalMap, 
					     const std::map<int,float> & XtalMapCurved  )
{            

  const EcalIntercalibConstantMap& icalMap = ical->getMap();
  float adcToGeV = float(agc->getEEValue());
  
  //numberOfSuperClusters = myTreeVariables_.nSuperClusters;
  //numberOfClusters      = myTreeVariables_.nClusters;

  /////////////////////////////////////////////////////////////////////////////////////////


      ///////////////////////////////////////////////////////////////////////////////////////
      // loop on endcap basic clusters 

   bool getC = false;
   for (reco::BasicClusterCollection::const_iterator clus = theEndcapBasicClusters->begin () ;
        clus != theEndcapBasicClusters->end () && numberOfClusters<MAXC;  ++clus)  {

          double energy = clus->energy () ;
	  double phi    = clus->phi () ;
	  double eta    = clus->eta () ;
	  if  ( energy/cosh(eta) <minEtEE_ ) continue;
	  double sinTheta         = fabs( sin( 2 *atan( exp(-1*clus->eta()) ) ) );
	  double transverseEnergy = clus->energy () * sinTheta;
	  double time = -1000.0 ; // gfdoc: work on this to provide a combination of crystals?
	  double ampli = 0. ;
	  double secondMin = 0. ;
	  double secondTime = -1000. ;
	  int numberOfXtalsInCluster=0 ;//counter for all xtals in cluster
	  int numXtalsinCluster = 0 ;   //xtals in cluster above 3sigma noise  
	  EEDetId maxDet ;
          EEDetId secDet ;
       
          // matching SC from objects
          float jidx = 100. ;
          getC = false ;
          for ( size_t i=0; i< selectedJets.size(); i++) {
               bool matched = matching( clus->energy(), clus->position(), selectedJets[i]->vertex(), selectedJets[i]->p4() ) ;
               jidx += 0.1 ;
               if ( matched ) { 
                  getC = true ;
                  //cout <<" ** from jet "<<i <<" / "<< selectedJets.size() ;
                  break ;
               }
          }
          if ( !getC ) continue ;

          //number of superClusters in event (collection = vector!)
	  //myTreeVariables_.nSuperClusters += theEndcapSuperClusters->size () ;
	  //myTreeVariables_.nEndcapSuperClusters = theEndcapSuperClusters->size () ;
	  myTreeVariables_.SCPIdx[numberOfSuperClusters] = jidx ;
          //cout<<" EJ Cl "<<numberOfClusters<<"  objId: "<< jidx <<" ->"<< myTreeVariables_.SCPIdx[numberOfSuperClusters]<<endl;
	  myTreeVariables_.superClusterEta[numberOfSuperClusters] = -99999 ;
	  myTreeVariables_.superClusterPhi[numberOfSuperClusters] = -99999 ;
	  myTreeVariables_.superClusterX[numberOfSuperClusters] = -99999 ;
	  myTreeVariables_.superClusterY[numberOfSuperClusters] = -99999 ;
	  myTreeVariables_.superClusterZ[numberOfSuperClusters] = -99999 ;
	  myTreeVariables_.superClusterRawEnergy[numberOfSuperClusters] = -999999 ;
	  myTreeVariables_.superClusterPhiWidth[numberOfSuperClusters]  = -999999 ;
	  myTreeVariables_.superClusterEtaWidth[numberOfSuperClusters]  = -999999 ; 
          myTreeVariables_.sMin[numberOfSuperClusters] = -1 ;
          myTreeVariables_.sMaj[numberOfSuperClusters] = -1 ;
          // GFdoc clusterDetIds holds crystals that participate to this basic cluster 
          std::vector<std::pair<DetId, float> > clusterDetIds = clus->hitsAndFractions() ; //get these from the cluster
       
          //////////////////////////////////////////////////////////////////////////////////////
          //loop on xtals in cluster
          for (std::vector<std::pair<DetId, float> >::const_iterator detitr = clusterDetIds.begin () ; 
	      detitr != clusterDetIds.end ()  && numberOfXtalsInCluster<MAXXTALINC; // && numberOfXtals<MAXXTAL ; 
	      ++detitr)// loop on rechics of barrel basic clusters
	   {

	      // Here I use the "find" on a digi collection... I have been warned...   (GFdoc: ??)
   	      // GFdoc: check if DetId belongs to ECAL; if so, find it among those if this basic cluster
    	     if ( (detitr -> first).det () != DetId::Ecal) 
	        { 
   	          std::cout << " det is " << (detitr -> first).det () << " (and not DetId::Ecal)" << std::endl ;
	          continue ;
	        }
	     if ( (detitr -> first).subdetId () != EcalEndcap) 
	        {
	          std::cout << " subdet is " << (detitr -> first).subdetId () << " (and not EcalEndcap)" << std::endl ; 
	          continue ; 
	        }
	   
	     // GFdoc now find it!
             EcalRecHitCollection::const_iterator thishit = theEndcapEcalRecHits->find ( (detitr -> first) ) ;
             if (thishit == theEndcapEcalRecHits->end ()) 
	        {
	          continue ;
	        }
	     //The checking above should no longer be needed... 
	     //as only those in the cluster would already have rechits..
	   
	     // GFdoc this is one crystal in the basic cluster
	     EcalRecHit myhit = (*thishit) ;
	   
             // SIC Feb 14 2011 -- Add check on RecHit flags (takes care of spike cleaning in 42X)
             if( !( myhit.checkFlag(EcalRecHit::kGood) ||
                    myhit.checkFlag(EcalRecHit::kOutOfTime) ||
                    myhit.checkFlag(EcalRecHit::kPoorCalib) 
		    ) 
		 )
                 continue;

	     // thisamp is the EE amplitude of the current rechit
	     double thisamp  = myhit.energy () ;
	     double thistime = myhit.time ();
	     double thisChi2 = myhit.chi2 ();
	     double thisOutOfTimeChi2 = myhit.outOfTimeChi2 ();

	   
	     EcalIntercalibConstantMap::const_iterator icalit = icalMap.find(detitr->first);
	     EcalIntercalibConstant icalconst = 1;
	     if( icalit!=icalMap.end() ) {
	       icalconst = (*icalit);
	     } else {
	       edm::LogError("EcalTimePhyTreeMaker") << "No intercalib const found for xtal "
						<< (detitr->first).rawId();
   	     }
	   
	     // get laser coefficient
	     float lasercalib = 1.;
	     lasercalib = laser->getLaserCorrection( detitr->first, iEvent.time());
	     // std::cout << "GF EE debug: " << lasercalib << std::endl;

	     // discard rechits with A/sigma < 12
	     if ( thisamp/(icalconst*lasercalib*adcToGeV) < (2.2*12) ) continue;

	     if (thisamp > 0.027) //cut on energy->number of crystals in cluster above 3sigma noise; gf: desirable?
             {  
	        numXtalsinCluster++ ; 
             }
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

	   
	    if(myhit.isTimeErrorValid())
              myTreeVariables_.xtalInBCTimeErr[numberOfClusters][numberOfXtalsInCluster]= myhit.timeError();
            else
              myTreeVariables_.xtalInBCTimeErr[numberOfClusters][numberOfXtalsInCluster]= -999999;

	    // xtal variables inside a barrel basic cluster 
	    myTreeVariables_.xtalInBCEnergy[numberOfClusters][numberOfXtalsInCluster]=       (float) thisamp;
	    myTreeVariables_.xtalInBCTime[numberOfClusters][numberOfXtalsInCluster]=         (float) thistime; 
	    myTreeVariables_.xtalInBCHashedIndex[numberOfClusters][numberOfXtalsInCluster]=  EEDetId (detitr -> first).hashedIndex () ; 
	    myTreeVariables_.xtalInBCIEta[numberOfClusters][numberOfXtalsInCluster]=           -999999; 
	    myTreeVariables_.xtalInBCIPhi[numberOfClusters][numberOfXtalsInCluster]=           -999999; 
	    myTreeVariables_.xtalInBCIx[numberOfClusters][numberOfXtalsInCluster]=         EEDetId((detitr -> first)).ix();
	    myTreeVariables_.xtalInBCIy[numberOfClusters][numberOfXtalsInCluster]=         EEDetId((detitr -> first)).iy();
	    myTreeVariables_.xtalInBCFlag[numberOfClusters][numberOfXtalsInCluster]=         myhit.recoFlag(); 
	    myTreeVariables_.xtalInBCAmplitudeADC[numberOfClusters][numberOfXtalsInCluster]= (float) thisamp/(icalconst*lasercalib*adcToGeV);
	    myTreeVariables_.xtalInBCChi2[numberOfClusters][numberOfXtalsInCluster]=         thisChi2;
	    myTreeVariables_.xtalInBCOutOfTimeChi2[numberOfClusters][numberOfXtalsInCluster]=thisOutOfTimeChi2;
	    // note: SwissCross = 1 - E4/E1   
            myTreeVariables_.xtalInBCSwissCross[numberOfClusters][numberOfXtalsInCluster] =
             EcalTools::swissCross(detitr->first,*theEndcapEcalRecHits,0.5);

            GlobalPoint pos = theGeometry->getPosition((myhit).detid());
            myTreeVariables_.xtalInBCEta[numberOfClusters][numberOfXtalsInCluster]=      pos.eta();
            myTreeVariables_.xtalInBCPhi[numberOfClusters][numberOfXtalsInCluster]=      pos.phi();
  

	    numberOfXtalsInCluster++ ; // increment number of crystals in basic cluster
	   
	 } //end loop on rechits within barrel basic clusters
       //////////////////////////////////////////////////////
       
       float E2 = (float)ampli + (float)secondMin ;
       
       //Set some more values for this cluster
       myTreeVariables_.clusterMom[numberOfClusters] = numberOfSuperClusters ;
       myTreeVariables_.CPIdx[numberOfClusters] = jidx;
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
       if(FILL_SHAPE_VARS) {
	 myTreeVariables_.clusterE2x2[numberOfClusters] = lazyTools -> e2x2(*clus);
	 myTreeVariables_.clusterE3x2[numberOfClusters] = lazyTools -> e3x2(*clus);
	 myTreeVariables_.clusterE3x3[numberOfClusters] = lazyTools -> e3x3(*clus);
	 myTreeVariables_.clusterE4x4[numberOfClusters] = lazyTools -> e4x4(*clus);
	 myTreeVariables_.clusterE5x5[numberOfClusters] = lazyTools -> e5x5(*clus);
	 
	 myTreeVariables_.clusterE2x5Right[numberOfClusters]  = lazyTools -> e2x5Right(*clus);
	 myTreeVariables_.clusterE2x5Left[numberOfClusters]   = lazyTools -> e2x5Left(*clus);
	 myTreeVariables_.clusterE2x5Top[numberOfClusters]    = lazyTools -> e2x5Top(*clus);
	 myTreeVariables_.clusterE2x5Bottom[numberOfClusters] = lazyTools -> e2x5Bottom(*clus);
	 
	 myTreeVariables_.clusterCovEtaEta[numberOfClusters] = (lazyTools -> covariances(*clus))[0];
	 myTreeVariables_.clusterCovPhiPhi[numberOfClusters] = (lazyTools -> covariances(*clus))[2];
	 myTreeVariables_.clusterCovEtaPhi[numberOfClusters] = (lazyTools -> covariances(*clus))[1];
	 
	 myTreeVariables_.clusterLat[numberOfClusters] = (lazyTools -> lat( *clus ))[2];
	 myTreeVariables_.clusterEtaLat[numberOfClusters] = (lazyTools -> lat( *clus ))[0];
	 myTreeVariables_.clusterPhiLat[numberOfClusters] = (lazyTools -> lat( *clus ))[1];
	 
	 myTreeVariables_.clusterZernike20[numberOfClusters] = lazyTools -> zernike20(*clus);
	 myTreeVariables_.clusterZernike42[numberOfClusters] = lazyTools -> zernike42(*clus);
       }
       
       numberOfClusters++ ;
     } //end (BASIC)cluster loop

     if ( getC ) {
        //cout<<"    SC "<< numberOfSuperClusters <<endl; 
        myTreeVariables_.nClustersInSuperCluster[numberOfSuperClusters] = numberOfClusters - myTreeVariables_.nClusters ;
        numberOfSuperClusters++ ;      
        myTreeVariables_.nSuperClusters = numberOfSuperClusters;
        myTreeVariables_.nClusters = numberOfClusters;
     }
     return ;
} // end dumpEndcapClusterInfo  


void
EcalTimePhyTreeMaker::dumpVertexInfo(const reco::VertexCollection* recVtxs, EcalTimePhyTreeContent & myTreeVariables_){

  int thisVertex=0;
  myTreeVariables_.nVertices= recVtxs->size();
  
  // make sure you don't exceed with number of vertices the available slots in the ntuple format (thisVertex) 
  for(reco::VertexCollection::const_iterator v=recVtxs->begin(); 
      v!=recVtxs->end() && thisVertex<MAXVTX; 
      ++v){
    myTreeVariables_.vtxNTracks[thisVertex]=v->tracksSize();
    myTreeVariables_.vtxIsFake[thisVertex]= v->isFake();
    myTreeVariables_.vtxChi2[thisVertex]=   v->chi2();
    myTreeVariables_.vtxNdof[thisVertex]=   v->ndof();
    myTreeVariables_.vtxX[thisVertex]=      v->x();
    myTreeVariables_.vtxDx[thisVertex]=     v->xError();
    myTreeVariables_.vtxY[thisVertex]=      v->y();
    myTreeVariables_.vtxDy[thisVertex]=     v->yError();
    myTreeVariables_.vtxZ[thisVertex]=      v->z();
    myTreeVariables_.vtxDz[thisVertex]=     v->zError();
    thisVertex++;
  }// end loop on vertices
  
}


// GFdoc GT information, at present anf +-1 bx
void
EcalTimePhyTreeMaker::dump3Ginfo (const edm::Event& iEvent, 
                             const edm::EventSetup& iSetup,
                             EcalTimePhyTreeContent & myTreeVariables_)
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

  //SIC - July 5 2010
  // Fill the trigger bit arrays, from EcalTimingAnalysis
  myTreeVariables_.l1NActiveTriggers = 0 ;
  myTreeVariables_.l1NActiveTechTriggers = 0;
  //I (Jason) Decided ONLY to look at the L1 triggers that took part in the decision, not just the ACTIVE triggers
  // HOPEFULLY this wasn't a bad decision
  edm::Handle< L1GlobalTriggerReadoutRecord > gtRecord;
  iEvent.getByLabel(gtRecordCollectionTag_, gtRecord);
  DecisionWord dWord = gtRecord->decisionWord();   // this will get the decision word *before* masking disabled bits
  int iBit = -1;
  for (std::vector<bool>::iterator itBit = dWord.begin(); itBit != dWord.end(); ++itBit)
  {
    iBit++;
    if (*itBit)
    {
      myTreeVariables_.l1ActiveTriggers[myTreeVariables_.l1NActiveTriggers] = iBit ;
      myTreeVariables_.l1NActiveTriggers++;
    }
  }

  TechnicalTriggerWord tw = gtRecord->technicalTriggerWord();
  if ( ! tw.empty() )
  {
    // loop over dec. bit to get total rate (no overlap)
    for ( int itechbit = 0; itechbit < 64; ++itechbit )
    {
      myTreeVariables_.l1ActiveTechTriggers[myTreeVariables_.l1NActiveTechTriggers] = 0; // ADD THIS 

      if ( tw[itechbit] )
      {
        myTreeVariables_.l1ActiveTechTriggers[myTreeVariables_.l1NActiveTechTriggers] = itechbit;
        myTreeVariables_.l1NActiveTechTriggers++;
      }

    }
  }
 
  return ;
} //PG dump3Ginfo  

bool EcalTimePhyTreeMaker::matching ( float sc_E, math::XYZPoint sc_pos, math::XYZPoint obj_v3, math::XYZTLorentzVector obj_p4 ) 
{
    math::XYZVector sc_v3( sc_pos.X() - obj_v3.X(), sc_pos.Y() - obj_v3.Y(), sc_pos.Z() - obj_v3.Z() );
    double dR =  ROOT::Math::VectorUtil::DeltaR( sc_v3, obj_p4 ) ;
    double dE =  obj_p4.E() - sc_E ;

    bool matched = ( dR < 0.5 && dE > 0 ) ? true : false ;

    return matched ;
}


// -------------------------------------------------------------------------------------------------------------

//PG stolen from Toyoko's EcalCosmicHists, PGP changed it (added Bx information)
std::vector<bool> 
EcalTimePhyTreeMaker::determineTriggers (const edm::Event& iEvent, const edm::EventSetup& eventSetup, int Bx)
{
  std::vector<bool> l1Triggers; //DT,CSC,RPC,HCAL,ECAL
                                //0 , 1 , 2 , 3  , 4
  for(int i=0;i<5;i++)
    l1Triggers.push_back(false);

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

  //cout << "l1Triggers: "; 
  //for(int i=0;i<5;i++)
  //  cout << l1Triggers[i];
  //cout << endl;
  return l1Triggers;
}
// -------------------------------------------------------------------------------------------------------------
