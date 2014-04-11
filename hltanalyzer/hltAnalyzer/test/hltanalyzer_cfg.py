import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("L1TriggerConfig.L1GtConfigProducers.L1GtConfig_cff")

from RSGravToZZ_kMpl01_M_1000_PU20bx25 import source
process.source = source
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.genElectrons = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("genParticles"), 
    cut = cms.string( "abs(pdgId())==11 && mother(0).pdgId()==23" ),
    filter = cms.bool(True)
)
process.genEleFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("genElectrons"),
    minNumber = cms.uint32(3)
)
process.genElectron1 = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("genElectrons"),
    cut = cms.string( "charge > 0" ),
)
process.genElectron2 = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("genElectrons"),
    cut = cms.string( "charge < 0" ),
)
process.demo = cms.EDAnalyzer('HltAnalyzer',
    hltLabel = cms.InputTag("TriggerResults","","HLT"),
    procName = cms.string("HLT"),
    hltPaths = cms.vstring("HLT_DoubleEle33_"),
    elesTag  = cms.InputTag('genElectrons'),
    e1Tag    = cms.InputTag('genElectron1'),
    e2Tag    = cms.InputTag('genElectron2')
)
process.p = cms.Path(process.genElectrons + ~process.genEleFilter + process.genElectron1+process.genElectron2+process.demo)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('histoGenElectrons.root')
)
