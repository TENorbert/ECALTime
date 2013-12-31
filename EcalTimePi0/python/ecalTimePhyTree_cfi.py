import FWCore.ParameterSet.Config as cms

ecalTimePhyTree = cms.EDAnalyzer("EcalTimePhyTreeMaker",
    barrelEcalRecHitCollection = cms.InputTag("reducedEcalRecHitsEB",""),
    endcapEcalRecHitCollection = cms.InputTag("reducedEcalRecHitsEE",""),

    useRaw = cms.bool(False),

    barrelEcalUncalibratedRecHitCollection = cms.InputTag("ecalRatioUncalibRecHit","EcalUncalibRecHitsEB"),
    endcapEcalUncalibratedRecHitCollection = cms.InputTag("ecalRatioUncalibRecHit","EcalUncalibRecHitsEE"),

    # SUPER and BASIC cluster collections come from AOD (no longer in-house made 3x3's)
    barrelBasicClusterCollection = cms.InputTag("hybridSuperClusters","hybridBarrelBasicClusters"),
    endcapBasicClusterCollection = cms.InputTag("multi5x5BasicClusters","multi5x5EndcapBasicClusters"),
    endcapSuperClusterCollection = cms.InputTag("correctedMulti5x5SuperClustersWithPreshower",""),
    barrelSuperClusterCollection = cms.InputTag("correctedHybridSuperClusters",""),
    JetSource = cms.InputTag("kt6PFJets"), # Romans use these jets
    METSource = cms.InputTag("pfMet"),
    MuonSource = cms.InputTag("muons"),
    ElectronSource = cms.InputTag("gsfElectrons"),
    PhotonSource   = cms.InputTag("photons"),
    triggerSource    = cms.InputTag("TriggerResults","","HLT"),
   # triggerHeader  = cms.untracked.string('HLT_Photon'),
   # triggerBody    = cms.untracked.string('_CaloIdVL_IsoL'),

    triggerHeader  = cms.untracked.string('HLT_DisplacedPhoton65'),
    triggerBody    = cms.untracked.string('_CaloIdVL_IsoL_PFMET25'),
    
    # jet/MET cuts                  Pt  eta  NJet  MaxNJet, MET          
    jetCuts          = cms.vdouble( 30, 2.4,    3,      99,  20 ),
    # photon cuts                   Pt  eta  hcal   dR  nPho sMinMin sMinMax                                 
    photonCuts       = cms.vdouble( 60, 2.4,    6, 0.3,   1,    0.1,    0.53 ),
    # photon Isolation              trk,  ecalEt, ecalR, hcalEt, hcalR
    photonIso        = cms.vdouble( 0.2,     4.5,   0.1,    4.0,   0.1 ),
    # electron cuts                 Pt  eta  iso   dR(photon, jets)                                 
    electronCuts     = cms.vdouble( 25, 2.4, 0.15,    0.3 ),
    # muon cuts                     Pt  eta  iso  dR(photon, jets)                                 
    muonCuts         = cms.vdouble( 25, 2.1, 0.2, 0.3 ),
    vertexCollection  = cms.InputTag("offlinePrimaryVertices",""),
    doTimeVSAmpliCorrection = cms.bool(True),
    muonCollection = cms.InputTag("GLBMuons"),
    hbTreshold = cms.double(1.),                               
    L1GlobalReadoutRecord = cms.string('gtDigis'),
    GTRecordCollection = cms.string('gtDigis'),
    runNum = cms.int32(-1),
    minEtEB = cms.double(2.),
    minEtEE = cms.double(2.),
    fileName = cms.string('EcalTimeTree'),
    TrackAssociatorParameters = cms.PSet(
    muonMaxDistanceSigmaX = cms.double(0.0),
    muonMaxDistanceSigmaY = cms.double(0.0),
    CSCSegmentCollectionLabel = cms.InputTag("cscSegments"),
    dRHcal = cms.double(9999.0),
    dREcal = cms.double(9999.0),
    dRPreshowerPreselection = cms.double(9999.0),
    CaloTowerCollectionLabel = cms.InputTag("towerMaker"),
    useEcal = cms.bool(True),
    usePreshower = cms.bool(True),
    dREcalPreselection = cms.double(0.05),
    HORecHitCollectionLabel = cms.InputTag("horeco"),
    dRMuon = cms.double(9999.0),
    crossedEnergyType = cms.string('SinglePointAlongTrajectory'),
    propagateAllDirections = cms.bool(True),
    muonMaxDistanceX = cms.double(5.0),
    muonMaxDistanceY = cms.double(5.0),
    useHO = cms.bool(True),
    trajectoryUncertaintyTolerance = cms.double(-1),
    accountForTrajectoryChangeCalo = cms.bool(False),
    DTRecSegment4DCollectionLabel = cms.InputTag("dt4DSegments"),
    EERecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
    dRHcalPreselection = cms.double(0.2),
    useMuon = cms.bool(True),
    useCalo = cms.bool(False),
    EBRecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
    dRMuonPreselection = cms.double(0.2),
    truthMatch = cms.bool(False),
    HBHERecHitCollectionLabel = cms.InputTag("hbhereco"),
    useHcal = cms.bool(True)
    )
)
