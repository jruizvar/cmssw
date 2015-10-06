import FWCore.ParameterSet.Config as cms

options = dict()

options['HLTProcessName']      = "HLT"
options['TnPPATHS']            = ["HLT_Ele27_eta2p1_WPLoose_Gsf_v2",]
options['TnPHLTTagFilters']    = ["hltEle27WPLooseGsfTrackIsoFilter",]
options['TnPHLTProbeFilters']  = ["*"]
options['HLTFILTERTOMEASURE']  = "hltEle105CaloIdVTGsfTrkIdTGsfDphiFilter"
options['GLOBALTAG']           = '74X_dataRun2_v2'
options['ELECTRON_CUTS']       = "(abs(superCluster.eta)<2.5) && (ecalEnergy*sin(superClusterPosition.theta)>10.0)"
options['ELECTRON_TAG_CUTS']   = "(abs(superCluster.eta)<=2.5) && !(1.4442<=abs(superCluster.eta)<=1.566) && pt >= 30.0"

###################################################################

process = cms.Process("tnp")

process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
process.hltHighLevel.throw = cms.bool(True)
process.hltHighLevel.HLTPaths = options['TnPPATHS']

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = options['GLOBALTAG']

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# source
import sys
SAMPLE = str(sys.argv[2])
process.load("ExoDiBosonResonances.EDBRCommon.PromptReco.Run2015D.SingleElectron."+SAMPLE)
import FWCore.PythonUtilities.LumiList as LumiList
process.source.lumisToProcess = LumiList.LumiList(
    filename = 'Cert_246908-257599_13TeV_PromptReco_Collisions15_25ns_JSON.txt').getVLuminosityBlockRange()

###################################################################
## ELECTRON MODULES
###################################################################

process.goodElectrons = cms.EDFilter("PATElectronRefSelector",
                                     src = cms.InputTag("slimmedElectrons"),
                                     cut = cms.string(options['ELECTRON_CUTS']) 
                                     )

process.eleVarHelper = cms.EDProducer("ElectronVariableHelper",
                                      probes = cms.InputTag("slimmedElectrons"),
                                      vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices")
                                      )

process.miniIsoElectronsMap = cms.EDProducer("PatElectronMiniIsolationValueMap",
                                             r_iso_min = cms.double(0.05),
                                             r_iso_max = cms.double(0.2),
                                             kt_scale  = cms.double(10.),
                                             charged_only = cms.bool(False),
                                             leptons = cms.InputTag("slimmedElectrons"),
                                             pfCands = cms.InputTag("packedPFCandidates"))

process.miniIsolatedGoodElectrons = cms.EDProducer("MiniIsolationSelector",
                                                   input         = cms.InputTag("goodElectrons"),
                                                   cut           = cms.string(options['ELECTRON_TAG_CUTS']),
                                                   selection     = cms.InputTag("miniIsoElectronsMap"),
                                                   id_cut        = cms.double(0.1),
                                                   isGreaterThan = cms.bool(False)
                                                   )

###################################################################
## ID MODULES
###################################################################

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff',
                 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff']

for idmod in my_id_modules:
    setupAllVIDIdsInModule(process, idmod, setupVIDElectronSelection)

process.goodElectronsHeepMiniIso = cms.EDProducer("PatElectronNm1Selector",
                                            input     = cms.InputTag("miniIsolatedGoodElectrons"),
                                            cut       = cms.string(""),
                                            selection = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"),
                                            cutNamesToMask = cms.vstring("GsfEleTrkPtIsoCut_0", "GsfEleEmHadD1IsoRhoCut_0")
                                            )

process.goodElectronsCutBasedTight = cms.EDProducer("PatElectronSelectorByValueMap",
                                            input     = cms.InputTag("goodElectrons"),
                                            cut       = cms.string(options['ELECTRON_TAG_CUTS']),
                                            selection = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
                                            id_cut    = cms.bool(True),
                                            )

###################################################################
## TRIGGER MATCHING
###################################################################

process.goodElectronsTagHLT = cms.EDProducer("PatElectronTriggerCandProducer",
                                             filterNames = cms.vstring(options['TnPHLTTagFilters']),
                                             inputs      = cms.InputTag("goodElectronsCutBasedTight"),
                                             bits        = cms.InputTag('TriggerResults::HLT'),
                                             objects     = cms.InputTag('selectedPatTrigger'),
                                             dR          = cms.double(0.3),
                                             isAND       = cms.bool(True)
                                             )

process.goodElectronsProbeHLT = cms.EDProducer("PatElectronTriggerCandProducer",
                                               filterNames = cms.vstring(options['TnPHLTProbeFilters']),
                                               inputs      = cms.InputTag("goodElectronsHeepMiniIso"),
                                               bits        = cms.InputTag('TriggerResults::HLT'),
                                               objects     = cms.InputTag('selectedPatTrigger'),
                                               dR          = cms.double(0.3),
                                               isAND       = cms.bool(True)
                                               )

process.goodElectronsProbeMeasureHLT = cms.EDProducer("PatElectronTriggerCandProducer",
                                                      filterNames = cms.vstring(options['HLTFILTERTOMEASURE']),
                                                      inputs      = cms.InputTag("goodElectronsProbeHLT"),
                                                      bits        = cms.InputTag('TriggerResults::HLT'),
                                                      objects     = cms.InputTag('selectedPatTrigger'),
                                                      dR          = cms.double(0.3),
                                                      isAND       = cms.bool(False)
                                                      )
  
###################################################################
## SEQUENCES
###################################################################

process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag("slimmedElectrons")
process.ele_sequence = cms.Sequence(
    process.goodElectrons                 +
    process.miniIsolatedGoodElectrons     +
    process.egmGsfElectronIDs             +
    process.goodElectronsHeepMiniIso      +
    process.goodElectronsCutBasedTight    +
    process.goodElectronsTagHLT           +
    process.goodElectronsProbeHLT         +
    process.goodElectronsProbeMeasureHLT  )
  
###################################################################
## TnP PAIRS
###################################################################

process.tagHeepMiniIsoHLT = cms.EDProducer("CandViewShallowCloneCombiner",
                                         decay = cms.string("goodElectronsTagHLT@+ goodElectronsProbeHLT@-"), 
                                         checkCharge = cms.bool(True),
                                         cut = cms.string("40<mass<1000"),
                                         )

process.allTagsAndProbes = cms.Sequence( process.tagHeepMiniIsoHLT )

##########################################################################
## TREE CONTENT
#########################################################################

ZVariablesToStore = cms.PSet(
    eta    = cms.string("eta"),
    abseta = cms.string("abs(eta)"),
    pt     = cms.string("pt"),
    mass   = cms.string("mass"),
    )   

ProbeVariablesToStore = cms.PSet(
    probe_Ele_eta    = cms.string("eta"),
    probe_Ele_abseta = cms.string("abs(eta)"),
    probe_Ele_pt     = cms.string("pt"),
    probe_Ele_et     = cms.string("et"),
    probe_Ele_e      = cms.string("energy"),
    probe_Ele_q      = cms.string("charge"),

#id based
    probe_Ele_dEtaIn   = cms.string("deltaEtaSuperClusterTrackAtVtx"),
    probe_Ele_dPhiIn   = cms.string("deltaPhiSuperClusterTrackAtVtx"),
    probe_Ele_dEtaSeed = cms.string("deltaEtaSeedClusterTrackAtVtx"),
    probe_Ele_hoe      = cms.string("hadronicOverEm"),
    probe_Ele_ooemoop  = cms.string("(1.0/ecalEnergy - eSuperClusterOverP/ecalEnergy)"),
    probe_Ele_mHits    = cms.InputTag("eleVarHelper:missinghits"),
    probe_Ele_dz       = cms.InputTag("eleVarHelper:dz"),
    probe_Ele_dxy      = cms.InputTag("eleVarHelper:dxy"),

#isolation
    probe_Ele_trkIso   = cms.string("dr03TkSumPt"),
    probe_Ele_miniIso  = cms.InputTag("miniIsoElectronsMap")
)

TagVariablesToStore = cms.PSet(
    Ele_eta    = cms.string("eta"),
    Ele_abseta = cms.string("abs(eta)"),
    Ele_pt     = cms.string("pt"),
    Ele_et     = cms.string("et"),
    Ele_e      = cms.string("energy"),
    Ele_q      = cms.string("charge"),
)

CommonStuffForGsfElectronProbe = cms.PSet(
    variables = cms.PSet(ProbeVariablesToStore),
    ignoreExceptions =  cms.bool (True),
    addRunLumiInfo   =  cms.bool (True),
    addEventVariablesInfo   =  cms.bool(True),
    vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    pairVariables =  cms.PSet(ZVariablesToStore),
    pairFlags     =  cms.PSet(
        mass60to120 = cms.string("60 < mass < 120")
        ),
    tagVariables   =  cms.PSet(TagVariablesToStore),
    tagFlags       =  cms.PSet(),    
)

##########################################################################
## TREE MAKER OPTIONS
##########################################################################

process.GsfElectronToTrigger = cms.EDAnalyzer("TagProbeFitTreeProducer",
                                              CommonStuffForGsfElectronProbe,
                                              isMC          = cms.bool(False),
                                              tagProbePairs = cms.InputTag("tagHeepMiniIsoHLT"),
                                              arbitration   = cms.string("None"),
                                              flags         = cms.PSet( passingHLT = cms.InputTag("goodElectronsProbeMeasureHLT") ),
                                              allProbes     = cms.InputTag("goodElectronsProbeHLT")
                                              )

process.tree_sequence = cms.Sequence ( process.GsfElectronToTrigger )

##########################################################################
##    ____       _   _     
##   |  _ \ __ _| |_| |__  
##   | |_) / _` | __| '_ \ 
##   |  __/ (_| | |_| | | |
##   |_|   \__,_|\__|_| |_|
##
##########################################################################

process.p = cms.Path(
    process.hltHighLevel        +
    process.miniIsoElectronsMap +
    process.ele_sequence        + 
    process.allTagsAndProbes    +
    process.eleVarHelper        +
    process.tree_sequence
    )

process.TFileService = cms.Service(
    "TFileService", fileName = cms.string( "TnPTree_data_"+SAMPLE+".root" )
    )
