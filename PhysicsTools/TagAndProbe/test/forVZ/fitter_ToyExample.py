import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.TnPMeasurement = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    InputFileNames = cms.vstring('TnPTree_data_SingleElectron_Run2015D.root'),
    InputDirectoryName = cms.string('GsfElectronToTrigger'),
    InputTreeName = cms.string('fitter_tree'),
    OutputFileName = cms.string('efficiency-passingHLT.root'),
    Variables = cms.PSet( 
        mass         = cms.vstring("Tag-Probe mass", "40.", "140.", "GeV"), 
        probe_Ele_et = cms.vstring("Probe e_{T}",    "0.",  "500.", "GeV") 
    ),
    Categories = cms.PSet( 
        passingHLT = cms.vstring('passingHLT', 'dummy[pass=1,fail=0]') 
    ),
    PDFs = cms.PSet(
        gaussPlusLinear = cms.vstring(
            "Gaussian::signal(mass, mean[91.19,80.,100.], sigma[10.,5.,20.])",
            "Chebychev::backgroundPass(mass, cPass[0,-1,1])",
            "Chebychev::backgroundFail(mass, cFail[0,-1,1])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        ),
    ),
    Efficiencies = cms.PSet(
        et = cms.PSet(
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet( 
                probe_Ele_et = cms.vdouble(35., 115., 500.) 
            ),
            EfficiencyCategoryAndState = cms.vstring('passingHLT', 'pass'),
            BinToPDFmap = cms.vstring("gaussPlusLinear")
        )
    ),
)

process.fit = cms.Path( process.TnPMeasurement )
