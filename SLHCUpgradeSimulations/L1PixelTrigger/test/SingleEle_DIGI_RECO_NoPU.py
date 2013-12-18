# Auto generated configuration file
# using: 
# Revision: 1.14 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: DIGI -s DIGI,L1,DIGI2RAW,RAW2DIGI,RECO:pixeltrackerlocalreco --magField 38T_PostLS1 --geometry ExtendedPhase2TkBE5D --conditions POSTLS261_V2::All --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1,SLHCUpgradeSimulations/Configuration/HCalCustoms.customise_HcalPhase0,SLHCUpgradeSimulations/Configuration/phase2TkCustomsBE5D.customise --eventcontent FEVTDEBUG --datatier GEN-SIM-DIGI-RAW-RECO --filein /store/group/comm_trigger/L1TrackTrigger/BE5D_612_SLHC6_patch1/singleEle/SingleEle_GEN_SIM.root --python_filename SingleEle_DIGI_RECO_NoPU.py -n 1 --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('RECO')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring('/store/group/comm_trigger/L1TrackTrigger/BE5D_612_SLHC6_patch1/singleEle/SingleEle_GEN_SIM.root')
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.14 $'),
    annotation = cms.untracked.string('DIGI nevts:1'),
    name = cms.untracked.string('Applications')
)

# Output definition

process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    fileName = cms.untracked.string('DIGI_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW-RECO')
    )
)

# Additional output definition
#################################################################################################
# Beam Spot
#################################################################################################
process.BeamSpotFromSim = cms.EDProducer("BeamSpotFromSimProducer")
process.BeamSpot = cms.Path(process.BeamSpotFromSim)

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'POSTLS261_V2::All', '')

# Path and EndPath definitions
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.raw2digi_step = cms.Path(process.RawToDigi)
process.reconstruction_step = cms.Path(process.pixeltrackerlocalreco)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

# --- Run the L1Calo and l1extra :
process.calolocalreco_step = cms.Path( process.calolocalreco )
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load("SLHCUpgradeSimulations.L1CaloTrigger.SLHCCaloTrigger_cff")
process.L1CaloTowerProducer.ECALDigis = cms.InputTag("simEcalTriggerPrimitiveDigis")
process.L1CaloTowerProducer.HCALDigis = cms.InputTag("simHcalTriggerPrimitiveDigis")
process.pL1Calo = cms.Path( process.SLHCCaloTrigger )
#################################################################################################
# Single Crystal
#################################################################################################
process.L1EGammaCrystalsProducer = cms.EDProducer("L1EGCrystalClusterProducer",
    DEBUG = cms.bool(False)
)
process.l1ExtraCrystalProducer = cms.EDProducer("L1ExtraCrystalPosition",
    eGammaSrc = cms.InputTag("SLHCL1ExtraParticles","EGamma"),
    eClusterSrc = cms.InputTag("L1EGammaCrystalsProducer","EGCrystalCluster")
)
process.crystal_producer = cms.Path(process.L1EGammaCrystalsProducer)
process.egcrystal_producer = cms.Path(process.l1ExtraCrystalProducer)

#################################################################################################
# Analyzer 
#################################################################################################
process.NtupleMaker = cms.EDAnalyzer('L1PixelTrigger',
    eGamma = cms.InputTag("SLHCL1ExtraParticles","EGamma"),
    eGammaCrystal = cms.InputTag("l1ExtraCrystalProducer","EGammaCrystal"),
    associatePixel = cms.bool(True),
    associateStrip = cms.bool(False),
    associateRecoTracks = cms.bool(False),
    ROUList = cms.vstring(
            'g4SimHitsTrackerHitsPixelBarrelLowTof', 
            'g4SimHitsTrackerHitsPixelBarrelHighTof', 
            'g4SimHitsTrackerHitsPixelEndcapLowTof', 
            'g4SimHitsTrackerHitsPixelEndcapHighTof')
)
process.p = cms.Path(process.NtupleMaker)
process.TFileService = cms.Service("TFileService", fileName = cms.string('SingleEle_NoPU_ntuple.root') )

# Schedule definition
process.schedule = cms.Schedule(process.digitisation_step,process.L1simulation_step,process.digi2raw_step,process.raw2digi_step,process.reconstruction_step,process.calolocalreco_step,process.pL1Calo,process.crystal_producer,process.egcrystal_producer,process.BeamSpot,process.p)

# customisation of the process.

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.phase2TkCustomsBE5D
from SLHCUpgradeSimulations.Configuration.phase2TkCustomsBE5D import customise 

#call to customisation function customise imported from SLHCUpgradeSimulations.Configuration.phase2TkCustomsBE5D
process = customise(process)

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.HCalCustoms
from SLHCUpgradeSimulations.Configuration.HCalCustoms import customise_HcalPhase0 

#call to customisation function customise_HcalPhase0 imported from SLHCUpgradeSimulations.Configuration.HCalCustoms
process = customise_HcalPhase0(process)

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.postLS1Customs
from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1 

#call to customisation function customisePostLS1 imported from SLHCUpgradeSimulations.Configuration.postLS1Customs
process = customisePostLS1(process)

# End of customisation functions
