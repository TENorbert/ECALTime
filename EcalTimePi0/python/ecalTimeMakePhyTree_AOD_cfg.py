import FWCore.ParameterSet.Config as cms

process = cms.Process("displacedPhotonsNtuplizer")

# to get clustering 
process.load('Configuration/StandardSequences/GeometryExtended_cff')

# Geometry
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cff")
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("Geometry.EcalMapping.EcalMapping_cfi")
process.load("Geometry.EcalMapping.EcalMappingRecord_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi") # gfwork: need this?


# Global Tag
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_noesprefer_cff")
process.GlobalTag.globaltag = 'GR_R_42_V19::All'


# Trigger
process.load("L1TriggerConfig.L1ScalesProducers.L1MuTriggerScalesConfig_cff")
process.load("L1TriggerConfig.L1ScalesProducers.L1MuTriggerPtScaleConfig_cff")
process.load("L1TriggerConfig.L1GtConfigProducers.L1GtBoardMapsConfig_cff")
process.load("L1TriggerConfig.L1GtConfigProducers.L1GtConfig_cff")
process.load("L1TriggerConfig.L1GtConfigProducers.Luminosity.startup.L1Menu_startup2_v2_Unprescaled_cff")
import FWCore.Modules.printContent_cfi
process.dumpEv = FWCore.Modules.printContent_cfi.printContent.clone()

import EventFilter.L1GlobalTriggerRawToDigi.l1GtUnpack_cfi
process.gtDigis = EventFilter.L1GlobalTriggerRawToDigi.l1GtUnpack_cfi.l1GtUnpack.clone()


# this is the ntuple producer
process.load("ECALTime.EcalTimePi0.ecalTimePhyTree_cfi")
process.ecalTimePhyTree.fileName = 'EcalTimeTree'
process.ecalTimePhyTree.barrelEcalRecHitCollection = cms.InputTag("reducedEcalRecHitsEB","")
process.ecalTimePhyTree.endcapEcalRecHitCollection = cms.InputTag("reducedEcalRecHitsEE","")
process.ecalTimePhyTree.barrelBasicClusterCollection = cms.InputTag("hybridSuperClusters","hybridBarrelBasicClusters")
process.ecalTimePhyTree.endcapBasicClusterCollection = cms.InputTag("multi5x5BasicClusters","multi5x5EndcapBasicClusters")
process.ecalTimePhyTree.barrelSuperClusterCollection = cms.InputTag("correctedHybridSuperClusters","")
process.ecalTimePhyTree.endcapSuperClusterCollection = cms.InputTag("correctedMulti5x5SuperClustersWithPreshower","")
process.ecalTimePhyTree.PhotonSource = cms.InputTag("myphotons")
#process.ecalTimePhyTree.PhotonSource = cms.InputTag("photons")
process.ecalTimePhyTree.muonCollection = cms.InputTag("muons")
process.ecalTimePhyTree.runNum = 999999



# Set up cuts for physics objects
# jet cuts                                           pt    eta  nJets
process.ecalTimePhyTree.jetCuts       = cms.vdouble( 25. , 2.4, 3 )
process.ecalTimePhyTree.metCuts       = cms.vdouble( 20  )
# photon cuts                                        pt  eta  dR   nPhoton
process.ecalTimePhyTree.photonCuts    = cms.vdouble( 30, 2.4, 0.3, 1 )
process.ecalTimePhyTree.electronCuts  = cms.vdouble( 25, 2.4, 0.15, 0.3 )
process.ecalTimePhyTree.muonCuts      = cms.vdouble( 25, 2.1, 0.2, 0.3 )

process.dumpEvContent = cms.EDAnalyzer("EventContentAnalyzer")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10000))




#################### ala romans ####################################
process.load("RecoEgamma.EgammaPhotonProducers.photonSequence_cff")
process.load("RecoEgamma.PhotonIdentification.photonId_cff")
process.load("RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi")
###########  USE UNCLEANED SUPERCLUSTERS  ######################### MS

process.load("RecoEcal.EgammaClusterProducers.uncleanSCRecovery_cfi") 
process.uncleanSCRecovered.cleanScCollection=cms.InputTag ("correctedHybridSuperClusters")	
process.photonCore.scHybridBarrelProducer=cms.InputTag ("uncleanSCRecovered:uncleanHybridSuperClusters")	

process.photons.barrelEcalHits=cms.InputTag("reducedEcalRecHitsEB")	
process.photons.endcapEcalHits=cms.InputTag("reducedEcalRecHitsEE")	

   
from RecoEgamma.PhotonIdentification.isolationCalculator_cfi import*	
newisolationSumsCalculator = isolationSumsCalculator.clone()	  
newisolationSumsCalculator.barrelEcalRecHitCollection = cms.InputTag('reducedEcalRecHitsEB')	
newisolationSumsCalculator.endcapEcalRecHitCollection = cms.InputTag('reducedEcalRecHitsEE')	
process.photons.isolationSumsCalculatorSet=newisolationSumsCalculator	




################################################################################# gf
import RecoEgamma.EgammaPhotonProducers.photonCore_cfi
import RecoEgamma.EgammaPhotonProducers.photons_cfi

process.myphotons=RecoEgamma.EgammaPhotonProducers.photons_cfi.photons.clone()
process.myphotons.barrelEcalHits=cms.InputTag("reducedEcalRecHitsEB")	
process.myphotons.endcapEcalHits=cms.InputTag("reducedEcalRecHitsEE")
process.myphotons.isolationSumsCalculatorSet=newisolationSumsCalculator	

process.myphotonCore=RecoEgamma.EgammaPhotonProducers.photonCore_cfi.photonCore.clone()
process.myphotonCore.scHybridBarrelProducer=cms.InputTag ("uncleanSCRecovered:uncleanHybridSuperClusters")

process.myPhotonSequence = cms.Sequence(process.myphotonCore+process.myphotons)

from RecoEgamma.PhotonIdentification.photonId_cfi import *
# photonID sequence
process.myPhotonIDSequence = cms.Sequence(PhotonIDProd)


###########  USE UNCLEANED SUPERCLUSTERS  ################ MS
process.uncleanPhotons = cms.Sequence(
               process.uncleanSCRecovered*
               #process.photonSequence *      # romans
               process.myPhotonSequence *     # gf
               #process.photonIDSequence *
               process.myPhotonIDSequence
               )


process.p = cms.Path(
    process.uncleanPhotons * 
    #process.dumpEvContent  *
    process.ecalTimePhyTree
    )


process.options   = cms.untracked.PSet(
                    wantSummary = cms.untracked.bool(True),
                    SkipEvent = cms.untracked.vstring('ProductNotFound')
)


process.MessageLogger = cms.Service("MessageLogger",
    cout = cms.untracked.PSet(
        threshold = cms.untracked.string('WARNING')
    ),
    categories = cms.untracked.vstring('ecalTimePhyTree'),
    destinations = cms.untracked.vstring('cout')
)
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)


# dbs search --query "find file where dataset=/ExpressPhysics/BeamCommissioning09-Express-v2/FEVT and run=124020" | grep store | awk '{printf "\"%s\",\n", $1}'
process.source = cms.Source(
    "PoolSource",
    skipEvents = cms.untracked.uint32(0),
    
    # a few files from:    /MinimumBias/Commissioning10-GR_R_35X_V7A_SD_EG-v2/RECO
    fileNames = (cms.untracked.vstring(
    'file:/data/franzoni/data/Run2011B-DoubleElectron-AOD-PromptReco-v1-000-179-889-F87DC321-BB01-E111-AD50-001D09F24DA8.root'
    )
                 ),
    # explicitly drop photons resident in AOD/RECO, to make sure only those locally re-made (uncleaned photons) are used
    inputCommands = cms.untracked.vstring('keep *'
                                          #,'drop  *_photonCore_*_RECO' # drop hfRecoEcalCandidate as remade in this process
                                          #, 'drop *_photons_*_RECO' # drop photons as remade in this process
                                          )
)
