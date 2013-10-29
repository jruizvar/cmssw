# Auto generated configuration file
# using: 
# Revision: 1.14 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step3 -s RAW2DIGI,RECO:pixeltrackerlocalreco --magField 38T_PostLS1 --geometry ExtendedPhase2TkBE5D --conditions POSTLS261_V2::All --customise SLHCUpgradeSimulations/Configuration/combinedCustoms.cust_phase2_BE5D --eventcontent FEVTDEBUG --datatier RECO --filein /store/group/comm_trigger/L1TrackTrigger/BE5D_612_SLHC6_patch1/singleEle/OldHcal/SingleElectron_BE5D_PU140.root --fileout file:dummy.root --python_filename SingleEle_RERECO_noPU.py -n 1 --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('ntuple')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring('/store/group/comm_trigger/L1TrackTrigger/BE5D_612_SLHC6_patch1/singleEle/OldHcal/SingleElectron_BE5D_PU140.root')
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
process.raw2digi_step = cms.Path(process.RawToDigi)
process.reconstruction_step = cms.Path(process.pixeltrackerlocalreco)
process.calolocalreco_step = cms.Path( process.calolocalreco )

# --- Run the L1Calo and l1extra :
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load("SLHCUpgradeSimulations.L1CaloTrigger.SLHCCaloTrigger_cff")
process.L1CaloTowerProducer.ECALDigis = cms.InputTag("simEcalTriggerPrimitiveDigis")
process.L1CaloTowerProducer.HCALDigis = cms.InputTag("simHcalTriggerPrimitiveDigis")
process.pL1Calo = cms.Path( process.SLHCCaloTrigger )
#################################################################################################
# Analyzer 
#################################################################################################
process.NtupleMaker = cms.EDAnalyzer('L1PixelTrigger',
    eGamma = cms.InputTag("SLHCL1ExtraParticles","EGamma"),
    eGammaCrystal = cms.InputTag("SLHCL1ExtraParticlesNewClustering","EGamma")
)
process.p = cms.Path(process.NtupleMaker)
process.TFileService = cms.Service("TFileService", fileName = cms.string('SingleElectron_PU140_ntuple.root') )

# Schedule definition
process.schedule = cms.Schedule(process.pL1Calo,process.BeamSpot,process.reconstruction_step,process.p)

# customisation of the process.
from SLHCUpgradeSimulations.Configuration.phase2TkCustomsBE5D import customise
process=customise(process)
