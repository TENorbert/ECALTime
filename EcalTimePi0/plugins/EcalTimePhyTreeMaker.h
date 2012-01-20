#ifndef EcalTimePhyTreeMaker_h
#define EcalTimePhyTreeMaker_h

// -*- C++ -*-
//
// Package:   EcalTimePhyTreeMaker
// Class:     EcalTimePhyTreeMaker
//
/**\class EcalTimePhyTreeMaker EcalTimePhyTreeMaker.h

Description: <one line class summary>

Implementation:
<Notes on implementation>
 */
//
// Authors:                   Shih-Chuan Kao, Giovanni Franzoni (UMN)
//         Created:  Mo Jul 14 5:46:22 CEST 2008
// $Id: EcalTimePhyTreeMaker.h,v 1.5 2012/01/15 09:35:08 franzoni Exp $
//
//

// system include files
#include <memory>
#include <vector>
#include <map>
#include <set>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuRegionalCand.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTReadoutCollection.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GtPsbWord.h"

#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"

#include "CondFormats/EcalObjects/interface/EcalIntercalibConstants.h"
#include "CondFormats/EcalObjects/interface/EcalADCToGeVConstant.h"
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbService.h"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TGraph.h"
#include "TTree.h"

#include <vector>

// *** for TrackAssociation
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
// ***

//for track length
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "RecoCaloTools/Navigation/interface/CaloNavigator.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElement.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"

#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"

#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"

// containers for vertices
#include <DataFormats/VertexReco/interface/VertexFwd.h>

#include <TMath.h>
#include <Math/VectorUtil.h>

#include "ECALTime/EcalTimePi0/interface/EcalTimePhyTreeContent.h"
#include "ECALTime/EcalTimePi0/interface/timeVsAmpliCorrector.h"


typedef std::pair<reco::SuperClusterRef, float> ParticleSC  ;

class EcalTimePhyTreeMaker : public edm::EDAnalyzer 
{
   public:

      explicit EcalTimePhyTreeMaker (const edm::ParameterSet&) ;
      ~EcalTimePhyTreeMaker () ;

   protected:

      virtual void beginJob () {} ;
      virtual void analyze (const edm::Event&, const edm::EventSetup&) ;
      virtual void endJob () ;
      virtual void beginRun(edm::Run const &, edm::EventSetup const &) ;

  private:

      std::string intToString (int num) ;
      void initHists (int) ;

      //bool HLTSelection( const edm::Event& iEvent ) ;
      int  HLTSelection( const edm::Event& iEvent, int cutVal,  std::string str_head = "HLT_Photon", std::string str_body = "_CaloIdVL_IsoL" ) ;
 
      int  HLTSelection( const edm::Event& iEvent ) ;

      bool dumpEvtObjectInfo( const edm::Event& iEvent ) ;

      //! dump Cluster information
      //! has to run after dumpMUinfo, to have the XtalMap already filled
      void dumpBarrelClusterInfo(const edm::Event& iEvent,
				 const CaloGeometry * theGeometry,
				 const CaloTopology * theCaloTopology,
				 const EcalRecHitCollection* theBarrelEcalRecHits,
				 const reco::BasicClusterCollection* theBarrelBasicClusters,
				 const reco::SuperClusterCollection* theBarrelSuperClusters,
				 EcalClusterLazyTools* lazyTools,
				 const std::map<int,float> & XtalMap,
				 const std::map<int,float> & XtalMapCurved,
				 EcalTimePhyTreeContent & myTreeVariables_) ;

      void dumpEndcapClusterInfo(const edm::Event& iEvent,
				 const CaloGeometry * theGeometry,
				 const CaloTopology * theCaloTopology,
				 const EcalRecHitCollection* theEndcapEcalRecHits,
				 const reco::BasicClusterCollection* theEndcapBasicClusters,
				 const reco::SuperClusterCollection* theEndcapSuperClusters,
				 EcalClusterLazyTools* lazyTools,
				 const std::map<int,float> & XtalMap,
				 const std::map<int,float> & XtalMapCurved,
				 EcalTimePhyTreeContent & myTreeVariables_) ;

      void dumpJetBarrelClusterInfo (const edm::Event& iEvent,
                                 const CaloGeometry * theGeometry,
				 const CaloTopology * theCaloTopology,
				 const EcalRecHitCollection* theBarrelEcalRecHits,
				 const reco::BasicClusterCollection* theBarrelBasicClusters,
				 EcalClusterLazyTools* lazyTools,
				 const std::map<int,float> & XtalMap,
				 const std::map<int,float> & XtalMapCurved  ) ;

      void dumpJetEndcapClusterInfo (const edm::Event& iEvent,
                                 const CaloGeometry * theGeometry,
				 const CaloTopology * theCaloTopology,
				 const EcalRecHitCollection* theEndcapEcalRecHits,
				 const reco::BasicClusterCollection* theEndcapBasicClusters,
				 EcalClusterLazyTools* lazyTools,
				 const std::map<int,float> & XtalMap,
				 const std::map<int,float> & XtalMapCurved  ) ;

      
      void dumpVertexInfo(const reco::VertexCollection* recVtxs, EcalTimePhyTreeContent & myTreeVariables_);
 
      //! dump trigger information
      void dump3Ginfo (const edm::Event& iEvent, const edm::EventSetup& eventSetup,
                       EcalTimePhyTreeContent & myTreeVariables_) ;
  
      bool matching ( float sc_Energy, math::XYZPoint sc_pos, math::XYZPoint obj_v3, math::XYZTLorentzVector obj_p4 ) ;

        //! collect trigger information to be dumped
      std::vector<bool> determineTriggers (const edm::Event& iEvent, 
                                           const edm::EventSetup& eventSetup, int Bx=0) ;



    // ----------member data ---------------------------

   protected:

      edm::InputTag barrelEcalRecHitCollection_ ;
      edm::InputTag endcapEcalRecHitCollection_ ;
      edm::InputTag barrelBasicClusterCollection_ ;
      edm::InputTag endcapBasicClusterCollection_ ;
      edm::InputTag barrelSuperClusterCollection_ ;
      edm::InputTag endcapSuperClusterCollection_ ;
      edm::InputTag muonCollection_ ;
      edm::InputTag JetSource_ ;
      edm::InputTag METSource_ ;
      edm::InputTag MuonSource_ ;
      edm::InputTag ElectronSource_ ;
      edm::InputTag PhotonSource_ ;
      edm::InputTag triggerSource_ ;
      edm::InputTag vertexCollection_ ;
      edm::InputTag l1GMTReadoutRecTag_ ;
      edm::InputTag gtRecordCollectionTag_ ;

      int runNum_ ;
      double minEtEB_;
      double minEtEE_;
      std::vector<double> jetCuts_ ;
      std::vector<double> metCuts_ ;
      std::vector<double> photonCuts_ ;
      std::vector<double> photonIso_ ;
      std::vector<double> electronCuts_ ;
      std::vector<double> muonCuts_ ;
      std::string fileName_ ;
      std::string triggerName_ ;

      EcalTimePhyTreeContent myTreeVariables_ ;

      TFile* file_ ;
      TTree* tree_ ;
      
      std::vector<int> l1Accepts_ ;
      std::vector<std::string> l1Names_ ;

      double hbTreshold_;
      edm::ESHandle<EcalIntercalibConstants> ical;
      edm::ESHandle<EcalADCToGeVConstant> agc;
      edm::ESHandle<EcalLaserDbService> laser;

      std::vector<const reco::PFJet*> selectedJets ;
      //std::vector<const reco::Electron*> selectedElectrons ;
      std::vector<const reco::GsfElectron*> selectedElectrons ;
      std::vector<const reco::Muon*> selectedMuons ;
      std::vector<const reco::Photon*> selectedPhotons ;
      std::vector<ParticleSC> sclist ;

      int numberOfSuperClusters ;
      int numberOfClusters ;

      bool          doTimeVSAmpliCorrection_;
      timeCorrector theTimeCorrector_;
      int  naiveId_ ; 

} ;

#endif
