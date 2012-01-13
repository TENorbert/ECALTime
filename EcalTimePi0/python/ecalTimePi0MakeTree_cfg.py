import FWCore.ParameterSet.Config as cms

process = cms.Process("TIMEPI0ANALYSIS")

# gfworks: to get clustering 
process.load('Configuration/StandardSequences/GeometryExtended_cff')

# Geometry
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cff")
# process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi") # gfwork: need this? 
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
# process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi") # gfwork: need this?
process.load("Geometry.EcalMapping.EcalMapping_cfi")
process.load("Geometry.EcalMapping.EcalMappingRecord_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi") # gfwork: need this?


# Global Tag
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_noesprefer_cff")
#process.GlobalTag.globaltag = 'CRAFT_ALL_V12::All'
process.GlobalTag.globaltag = 'GR09_31X_V5P::All' # gfwork: update this? 


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



# general basic- and super- clustering sequences
import RecoEcal.EgammaClusterProducers.multi5x5ClusteringSequence_cff

# 3x3 clustering for barrel
process.multi5x5BasicClustersTimePi0Barrel =  RecoEcal.EgammaClusterProducers.multi5x5BasicClusters_cfi.multi5x5BasicClusters.clone(
    # which regions should be clusterized
    doEndcap = cms.bool(False),
    doBarrel = cms.bool(True),

    # gfwork: this is standard, can go away 
    barrelHitProducer = cms.string('ecalRecHit'),
    barrelHitCollection = cms.string('EcalRecHitsEB'),
    endcapHitProducer = cms.string('ecalRecHit'),
    endcapHitCollection = cms.string('EcalRecHitsEE'),
    
    IslandBarrelSeedThr = cms.double(0.150),   # barrelSeedThreshold
    IslandEndcapSeedThr = cms.double(0.250),   # endcapSeedThreshold

    barrelClusterCollection = cms.string('multi5x5BarrelBasicClusters'),
    endcapClusterCollection = cms.string('multi5x5EndcapBasicClusters'),
    clustershapecollectionEE = cms.string('multi5x5EndcapShape'),
    clustershapecollectionEB = cms.string('multi5x5BarrelShape'),
    barrelShapeAssociation = cms.string('multi5x5BarrelShapeAssoc'),
    endcapShapeAssociation = cms.string('multi5x5EndcapShapeAssoc'),
    )


# 3x3 clustering for endcap
process.multi5x5BasicClustersTimePi0Endcap =  RecoEcal.EgammaClusterProducers.multi5x5BasicClusters_cfi.multi5x5BasicClusters.clone(
    # which regions should be clusterized
    doEndcap = cms.bool(True),
    doBarrel = cms.bool(False),

    barrelHitProducer = cms.string('ecalRecHit'),
    barrelHitCollection = cms.string('EcalRecHitsEB'),
    endcapHitProducer = cms.string('ecalRecHit'),
    endcapHitCollection = cms.string('EcalRecHitsEE'),
    
    IslandBarrelSeedThr = cms.double(0.150),              # endcapSeedThreshold
    IslandEndcapSeedThr = cms.double(0.200),             # barrelSeedThreshold

    barrelClusterCollection = cms.string('multi5x5BarrelBasicClusters'),
    endcapClusterCollection = cms.string('multi5x5EndcapBasicClusters'),
    clustershapecollectionEE = cms.string('multi5x5EndcapShape'),
    clustershapecollectionEB = cms.string('multi5x5BarrelShape'),
    barrelShapeAssociation = cms.string('multi5x5BarrelShapeAssoc'),
    endcapShapeAssociation = cms.string('multi5x5EndcapShapeAssoc'),
    )


# super clustering for the ECAL BARREL, staring from multi5x5 3x3 clusters
process.multi5x5SuperClustersTimePi0Barrel =  RecoEcal.EgammaClusterProducers.multi5x5SuperClusters_cfi.multi5x5SuperClusters.clone(
    doBarrel = cms.bool(True),
    doEndcaps = cms.bool(False),
    barrelClusterProducer = cms.string('multi5x5BasicClustersTimePi0Barrel'),
    barrelClusterCollection = cms.string('multi5x5BarrelBasicClusters'),
    endcapClusterProducer = cms.string('multi5x5BasicClustersTimePi0Endcap'),
    endcapClusterCollection = cms.string('multi5x5EndcapBasicClusters')
 )


# super clustering for the ECAL ENDCAP, staring from multi5x5 3x3 clusters
process.multi5x5SuperClustersTimePi0Endcap =  RecoEcal.EgammaClusterProducers.multi5x5SuperClusters_cfi.multi5x5SuperClusters.clone(
    doBarrel = cms.bool(False),
    doEndcaps = cms.bool(True),
    barrelClusterProducer = cms.string('multi5x5BasicClustersTimePi0Barrel'),
    barrelClusterCollection = cms.string('multi5x5BarrelBasicClusters'),
    endcapClusterProducer = cms.string('multi5x5BasicClustersTimePi0Endcap'),
    endcapClusterCollection = cms.string('multi5x5EndcapBasicClusters')
 )




# this is the ntuple producer
process.load("ECALTime.EcalTimePi0.ecalTimePi0Tree_cfi")
process.ecalTimePi0Tree.fileName = 'EcalTimePi0Tree'
process.ecalTimePi0Tree.muonCollection = cms.InputTag("muons")
process.ecalTimePi0Tree.runNum = 108645
# gfworks: replace these names
process.ecalTimePi0Tree.barrelSuperClusterCollection = cms.InputTag("multi5x5SuperClustersTimePi0Barrel","multi5x5BarrelSuperClusters")
process.ecalTimePi0Tree.endcapSuperClusterCollection = cms.InputTag("multi5x5SuperClustersTimePi0Endcap","multi5x5EndcapSuperClusters")
process.ecalTimePi0Tree.barrelBasicClusterCollection = cms.InputTag("multi5x5BasicClustersTimePi0Barrel","multi5x5BarrelBasicClusters")
process.ecalTimePi0Tree.endcapBasicClusterCollection = cms.InputTag("multi5x5BasicClustersTimePi0Endcap","multi5x5EndcapBasicClusters")
process.ecalTimePi0Tree.barrelClusterShapeAssociationCollection = cms.InputTag("multi5x5BasicClustersTimePi0Barrel","multi5x5BarrelShapeAssoc")
process.ecalTimePi0Tree.endcapClusterShapeAssociationCollection = cms.InputTag("multi5x5BasicClustersTimePi0Endcap","multi5x5EndcapShapeAssoc") 



process.dumpEvContent = cms.EDAnalyzer("EventContentAnalyzer")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.p = cms.Path(
    process.multi5x5BasicClustersTimePi0Barrel *
    process.multi5x5BasicClustersTimePi0Endcap *
    process.multi5x5SuperClustersTimePi0Barrel *
    process.multi5x5SuperClustersTimePi0Endcap *
    process.dumpEvContent  *
    process.ecalTimePi0Tree
    )



process.MessageLogger = cms.Service("MessageLogger",
    cout = cms.untracked.PSet(
        threshold = cms.untracked.string('DEBUG')
    ),
    categories = cms.untracked.vstring('ecalTimePi0Tree'),
    destinations = cms.untracked.vstring('cout')
)




# GF: some legacy reco files to test; replace w/ collision data
# dbs search --query "find file where dataset=/ExpressPhysics/BeamCommissioning09-Express-v2/FEVT and run=124020" | grep store | awk '{printf "\"%s\",\n", $1}'
process.source = cms.Source(
    "PoolSource",
    skipEvents = cms.untracked.uint32(0),

     fileNames = (cms.untracked.vstring("file:/data/franzoni/data/pi0TimeAnalysis/FE37D88F-D6E6-DE11-92D0-000423D9989E.root"
    
#     fileNames = (cms.untracked.vstring(
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/FE37D88F-D6E6-DE11-92D0-000423D9989E.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/FE2CF646-CBE6-DE11-A65B-001D09F251FE.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/F2514441-D3E6-DE11-901C-000423D944F0.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/F2487F56-CDE6-DE11-8D0F-001D09F2A690.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/F0E91277-CAE6-DE11-876B-000423D98DD4.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/F09D067F-D4E6-DE11-8FEF-001D09F29619.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/EE1CB2FD-CDE6-DE11-84DF-00304879FBB2.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/EC331195-D6E6-DE11-8DF5-000423D94E1C.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/EC2C4C60-D7E6-DE11-B94E-000423D998BA.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/E604EFF6-D5E6-DE11-8B51-000423D6CA42.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/E4430D00-DBE6-DE11-9984-000423D6CA42.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/E29F29FB-DCE6-DE11-ACB2-0030487A1990.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/DE2C6387-CFE6-DE11-BD52-000423D99394.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/DCE92EFD-CDE6-DE11-8D8D-0030487A3232.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/DC682406-D8E6-DE11-9CC9-000423D9939C.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/DC2842CD-D8E6-DE11-934E-000423D951D4.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/D8CDDAC0-D3E6-DE11-90EF-000423D991D4.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/D6F6B846-CBE6-DE11-BFE9-001D09F2305C.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/D67260CF-D8E6-DE11-9EAF-001D09F2516D.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/D66FA010-D3E6-DE11-85F6-001617DBD224.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/D65FB5F7-DCE6-DE11-B894-0030487A322E.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/D63619DD-CBE6-DE11-9F80-000423D986C4.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/D078A8E9-D0E6-DE11-A2D5-001D09F28E80.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/C88068CC-D8E6-DE11-9875-000423D94E70.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/C23A4931-DAE6-DE11-B360-000423D6B48C.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/BC5A20B7-DBE6-DE11-B62F-000423D991D4.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/B2BC96BB-DBE6-DE11-BF92-003048D2C108.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/B297AB7D-D4E6-DE11-9222-001D09F2426D.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/B0455557-CDE6-DE11-AB59-001D09F2527B.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/B019982E-D0E6-DE11-BAC2-000423D952C0.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/AE50E4A6-CCE6-DE11-B556-001D09F2A690.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/AE180245-CBE6-DE11-B940-000423D9A212.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/ACEE1931-DAE6-DE11-88BD-0019DB29C5FC.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/AC0C9231-DAE6-DE11-908D-001617E30D52.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/A8AF7552-D2E6-DE11-AD3B-000423D6006E.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/A80E386B-CFE6-DE11-913C-001617C3B5E4.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/A2FB0C31-D5E6-DE11-BFD4-001D09F252DA.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/A276C352-D2E6-DE11-952A-000423D94908.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/A0E651FB-CDE6-DE11-8F3E-001617C3B66C.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/9A130BA9-CCE6-DE11-BC12-001D09F24F65.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/98F8482A-D0E6-DE11-A7A8-000423D998BA.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/98D260D0-CEE6-DE11-B2DF-0030487D1BCC.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/98BF4DFE-D5E6-DE11-BE59-000423D6BA18.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/9485FD04-DBE6-DE11-A7E3-001D09F2AD4D.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/8CC62FEA-D0E6-DE11-837C-001D09F242EA.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/8C392E6B-D9E6-DE11-82D1-000423D9A212.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/80596BF8-DCE6-DE11-8C4D-0030487A3232.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/7ED370A4-D1E6-DE11-9C2B-000423D8F63C.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/7684C2C0-D3E6-DE11-A9C1-000423D8FA38.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/74820E10-D3E6-DE11-A69E-000423D98B08.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/722883F9-D5E6-DE11-B8A2-000423DD2F34.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/6ED30A05-D8E6-DE11-9D01-001617E30D12.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/68F75F77-CAE6-DE11-ABB4-000423D98C20.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/68A686B6-DDE6-DE11-B607-000423D99614.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/66928B2A-D0E6-DE11-BA0A-000423D944F8.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/5EF23848-CBE6-DE11-88EB-001D09F2423B.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/5E4A3072-DEE6-DE11-BBAA-000423D98F98.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/5AC0A861-D7E6-DE11-A1EF-000423D99394.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/5A478D2A-D0E6-DE11-803A-000423D94A20.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/5A080CBC-DBE6-DE11-985C-001D09F26509.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/588472EB-D0E6-DE11-9C23-001D09F24763.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/54B2CB6E-D9E6-DE11-BF52-000423D99658.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/5455E7FE-D5E6-DE11-B934-000423D6CA02.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/540D0DBC-DBE6-DE11-BB58-001D09F24D4E.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/52F962CF-CEE6-DE11-A1B7-000423D6CA72.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/50F3AFB2-DDE6-DE11-8795-000423D99B3E.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/4EFD6331-D5E6-DE11-9576-001D09F291D7.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/4EE80411-D3E6-DE11-8B81-000423D987FC.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/4A7B3501-D8E6-DE11-AE34-000423D99F1E.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/4A22DE52-D2E6-DE11-B8F1-000423D987E0.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/48F2B0DB-CBE6-DE11-9D14-000423D98750.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/48A381DC-CBE6-DE11-912E-000423D98E6C.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/467E8532-DAE6-DE11-AD1D-001617DBCF6A.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/465A6D7E-D4E6-DE11-82D9-001D09F252DA.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/40C27308-D8E6-DE11-991F-001D09F24763.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/384860FF-D5E6-DE11-A781-000423D98804.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/348B7BFD-CDE6-DE11-9AE8-00304879FA4A.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/2EB7537D-DCE6-DE11-948F-001D09F24EE3.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/2C6DE034-D5E6-DE11-900F-001D09F26C5C.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/2A8767A8-CCE6-DE11-BE44-001D09F29533.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/24A8B3D0-CEE6-DE11-8601-000423D999CA.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/2403EA33-DAE6-DE11-B905-003048D37538.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/20D5DC52-D2E6-DE11-A574-000423D998BA.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/207224E8-D0E6-DE11-ADC1-001D09F24EE3.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/1C871B8A-D1E6-DE11-8B5F-000423D98E6C.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/0EABC47E-DCE6-DE11-BB6F-000423D94524.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/0E554646-CBE6-DE11-B28E-0030487A18F2.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/0CA53501-CEE6-DE11-AF74-0030487A322E.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/0A1AA97F-DCE6-DE11-9A83-001D09F24489.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/08CE77B2-DDE6-DE11-86E6-000423D991F0.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/08B873FA-CDE6-DE11-B52A-000423D991D4.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/06C6BC2A-D0E6-DE11-A897-000423D94534.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/060388DB-CBE6-DE11-AD80-000423D944F8.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/046ECC7F-D4E6-DE11-845D-001D09F291D7.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/04389BCE-CEE6-DE11-A077-000423D6B42C.root",
#              "/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/124/020/02EB23BF-DBE6-DE11-9F9B-000423D991F0.root"
                                        )
                  )
    )



