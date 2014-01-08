import FWCore.ParameterSet.Config as cms

process = cms.Process("displacedPhotonsNtuplizer")

# Trigger
process.load("L1TriggerConfig.L1ScalesProducers.L1MuTriggerScalesConfig_cff")
process.load("L1TriggerConfig.L1ScalesProducers.L1MuTriggerPtScaleConfig_cff")
process.load("L1TriggerConfig.L1GtConfigProducers.L1GtBoardMapsConfig_cff")
process.load("L1TriggerConfig.L1GtConfigProducers.L1GtConfig_cff")
process.load("L1TriggerConfig.L1GtConfigProducers.Luminosity.startup.L1Menu_startup2_v2_Unprescaled_cff")

import EventFilter.L1GlobalTriggerRawToDigi.l1GtUnpack_cfi
process.gtDigis = EventFilter.L1GlobalTriggerRawToDigi.l1GtUnpack_cfi.l1GtUnpack.clone()

######################  NTUPLE PRODUCER  ####################################################################
# this is the ntuple producer
process.load("ECALTime.EcalTimePi0.ecalTimePhyTree_cfi")
process.ecalTimePhyTree.fileName = 'EcalTimeTree'
process.ecalTimePhyTree.barrelEcalRecHitCollection = cms.InputTag("reducedEcalRecHitsEB","")
process.ecalTimePhyTree.endcapEcalRecHitCollection = cms.InputTag("reducedEcalRecHitsEE","")
process.ecalTimePhyTree.barrelBasicClusterCollection = cms.InputTag("uncleanSCRecovered","uncleanHybridBarrelBasicClusters")
#process.ecalTimePhyTree.barrelBasicClusterCollection = cms.InputTag("hybridSuperClusters","hybridBarrelBasicClusters")
#process.ecalTimePhyTree.endcapBasicClusterCollection = cms.InputTag("multi5x5BasicClusters","multi5x5EndcapBasicClusters")
process.ecalTimePhyTree.endcapBasicClusterCollection = cms.InputTag("multi5x5SuperClusters","multi5x5EndcapBasicClusters")
process.ecalTimePhyTree.barrelSuperClusterCollection = cms.InputTag("uncleanSCRecovered","uncleanHybridSuperClusters")
#process.ecalTimePhyTree.barrelSuperClusterCollection = cms.InputTag("correctedHybridSuperClusters","")
process.ecalTimePhyTree.endcapSuperClusterCollection = cms.InputTag("correctedMulti5x5SuperClustersWithPreshower","")
process.ecalTimePhyTree.PhotonSource = cms.InputTag("myphotons")
#process.ecalTimePhyTree.PhotonSource = cms.InputTag("photons")
process.ecalTimePhyTree.muonCollection = cms.InputTag("muons")
# switch on or off Tambe's analysis level corrections
process.ecalTimePhyTree.doTimeVSAmpliCorrection = cms.bool(False) # True if in CMSSW4XY
process.ecalTimePhyTree.runNum = 999999
#process.ecalTimePhyTree.triggerHeader  = cms.untracked.string('HLT_Photon')
#process.ecalTimePhyTree.triggerBody    = cms.untracked.string('EBOnly_CaloIdVL_IsoL_TriPFJet25')
process.ecalTimePhyTree.triggerHeader  = cms.untracked.string('HLT_DisplacedPhoton65')
process.ecalTimePhyTree.triggerBody    = cms.untracked.string('_CaloIdVL_IsoL_PFMET25')
#process.ecalTimePhyTree.triggerBody    = cms.untracked.string('HLT_Photon50_CaloIdVL_IsoL','HLT_DisplacedPhoton65_CaloIdVL_IsoL_PFMET25'),
#process.ecalTimePhyTree.triggerBody    = cms.untracked.string('_CaloIdVL_IsoL')
process.ecalTimePhyTree.trigSource     = cms.InputTag("TriggerResults","","HLT")
process.ecalTimePhyTree.L1GlobalReadoutRecord = cms.string('gtDigis')


#################### PRE-SELECTION-CUTS ##################################

# Set up cuts for physics objects;
# nJets and nPhoton constitute event-based selections 

###################### signal preselection cuts ##########################
## jet cuts                                           pt   |eta|  nJets
#process.ecalTimePhyTree.jetCuts       = cms.vdouble( 25. , 2.4, 3 )
## met cuts                                           Et
#process.ecalTimePhyTree.metCuts       = cms.vdouble( 20  )
## photon cuts                                        pt |eta|  dR   nPhoton
#process.ecalTimePhyTree.photonCuts    = cms.vdouble( 30, 2.4, 0.3, 1 )
# photon Isolation                                  trk,  ecalEt, ecalR, hcalEt, hcalR
#process.ecalTimePhyTree.photonIso     = cms.vdouble( 0.2,    4.5,   0.1,    4.0,   0.1 )
## electron cuts                                      pt |eta| relIso dR
#process.ecalTimePhyTree.electronCuts  = cms.vdouble( 25, 2.4, 0.15, 0.3 )
## muon cuts                                          pt |eta| relIso dR
#process.ecalTimePhyTree.muonCuts      = cms.vdouble( 25, 2.1, 0.2, 0.3 )

##################### loose selection for studies ##########################
# jet cuts                                           pt  |eta|  NJet  MaxNJet MET  
process.ecalTimePhyTree.jetCuts       = cms.vdouble( 25,  2.4,    0,      99,    0 )
# photon cuts                                        Pt  eta  hcal   dR  nPho sMinMin sMinMax                                 
process.ecalTimePhyTree.photonCuts    = cms.vdouble( 60, 9.9,   99,   0.,   1,    0.,    1.0 )
# photon Isolation                                  trk,  ecalEt, ecalR, hcalEt, hcalR
process.ecalTimePhyTree.photonIso     = cms.vdouble( 0.5,    99.,   1.,    99.,     1. )
# electron cuts                                      pt |eta| relIso dR
process.ecalTimePhyTree.electronCuts  = cms.vdouble( 25, 2.4, 0.5,  0.5 )
# muon cuts                                          pt |eta| relIso dR
process.ecalTimePhyTree.muonCuts      = cms.vdouble( 25, 2.1, 0.5,  0.5 )


###########  USE UNCLEANED SUPERCLUSTERS  #############################################

################# RECO ON THE FLY FROM AOD ###################################
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('RecoEgamma.EgammaPhotonProducers.conversionTracks_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("RecoEcal.Configuration.RecoEcal_cff")
from Configuration.StandardSequences.Reconstruction_cff import *
from RecoEcal.Configuration.RecoEcal_cff import *
from RecoEcal.EgammaClusterProducers.hybridSuperClusters_cfi import *



# Global Tag
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_noesprefer_cff")
# global tag for 42x
#process.GlobalTag.globaltag = 'GR_P_V22::All'
# global tag for 44x
#process.GlobalTag.globaltag = 'GR_R_44_V13::All'
#process.GlobalTag.globaltag = 'GR_R_53_V18::All'
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag( process.GlobalTag, 'GR_R_53_V18::All' )


# to get clustering 
#Geometry

process.load("Configuration.Geometry.GeometryIdeal_cff")
#process.load("Configuration.StandardSequences.Geometry_cff")
process.load('Configuration/StandardSequences/GeometryExtended_cff')
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cff")
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("Geometry.EcalMapping.EcalMapping_cfi")
process.load("Geometry.EcalMapping.EcalMappingRecord_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi") # gfwork: need this?
process.CaloTowerConstituentsMapBuilder = cms.ESProducer("CaloTowerConstituentsMapBuilder")


##################### PHOTON PRODUCTION PROCESS ###############################################

process.load("RecoEgamma.PhotonIdentification.photonId_cff")
process.load("RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi")

import RecoEgamma.EgammaPhotonProducers.photonCore_cfi
import RecoEgamma.EgammaPhotonProducers.photons_cfi

process.myphotonCores=RecoEgamma.EgammaPhotonProducers.photonCore_cfi.photonCore.clone()
process.load("RecoEcal.EgammaClusterProducers.uncleanSCRecovery_cfi") 
process.uncleanSCRecovered.cleanScCollection=cms.InputTag ("correctedHybridSuperClusters")	
process.myphotonCores.scHybridBarrelProducer=cms.InputTag ("uncleanSCRecovered:uncleanHybridSuperClusters")

from RecoEgamma.PhotonIdentification.mipVariable_cfi import *
newMipVariable = mipVariable.clone()
newMipVariable.barrelEcalRecHitCollection = cms.InputTag('reducedEcalRecHitsEB')
newMipVariable.endcapEcalRecHitCollection = cms.InputTag('reducedEcalRecHitsEE')

from RecoEgamma.PhotonIdentification.isolationCalculator_cfi import*	
newisolationSumsCalculator = isolationSumsCalculator.clone()
newisolationSumsCalculator.barrelEcalRecHitProducer = cms.string('reducedEcalRecHitsEB')
newisolationSumsCalculator.endcapEcalRecHitProducer = cms.string('reducedEcalRecHitsEE')
newisolationSumsCalculator.barrelEcalRecHitCollection = cms.InputTag('reducedEcalRecHitsEB')	
newisolationSumsCalculator.endcapEcalRecHitCollection = cms.InputTag('reducedEcalRecHitsEE')

process.myphotons=RecoEgamma.EgammaPhotonProducers.photons_cfi.photons.clone()
process.myphotons.barrelEcalHits=cms.InputTag("reducedEcalRecHitsEB")	
process.myphotons.endcapEcalHits=cms.InputTag("reducedEcalRecHitsEE")
process.myphotons.isolationSumsCalculatorSet=newisolationSumsCalculator
process.myphotons.mipVariableSet = newMipVariable
process.myphotons.photonCoreProducer=cms.InputTag("myphotonCores")

process.myPhotonSequence = cms.Sequence(process.myphotonCores+
                                        process.myphotons)

from RecoEgamma.PhotonIdentification.photonId_cfi import *
# photonID sequence
process.myPhotonIDSequence = cms.Sequence(PhotonIDProd)
process.PhotonIDProd.photonProducer=cms.string("myphotons")

##### Unclean Photon Process ############################################
process.uncleanPhotons = cms.Sequence(
                process.uncleanSCRecovered *
                process.myPhotonSequence *
               process.myPhotonIDSequence
               )

process.dumpEvContent = cms.EDAnalyzer("EventContentAnalyzer")
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))

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
    #'file:/data/franzoni/data/Run2011B-PhotonHad-AOD-PromptReco-v1-000-179-558-5CDAF51F-A800-E111-ADD4-BCAEC518FF52.root'
    'file:/hdfs/cms/phedex/store/data/Run2012C/SinglePhoton/RECO/EXODisplacedPhoton-PromptSkim-v3/000/198/941/00000/0EA7C91A-B8CF-E111-9766-002481E150EA.root'
#    'file:/hdfs/cms/phedex/store/data/Run2012A/DoubleElectron/AOD/22Jan2013-v1/30000/4A0C8B83-6470-E211-BBA0-003048679182.root'
    )
                 ),
    # explicitly drop photons resident in AOD/RECO, to make sure only those locally re-made (uncleaned photons) are used
    inputCommands = cms.untracked.vstring('keep *'
                                          #,'drop  *_photonCore_*_RECO' # drop hfRecoEcalCandidate as remade in this process
                                          #, 'drop *_photons_*_RECO' # drop photons as remade in this process
                                          )
)
