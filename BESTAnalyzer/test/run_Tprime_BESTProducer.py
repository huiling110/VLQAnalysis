import FWCore.ParameterSet.Config as cms
from Configuration.AlCa.GlobalTag import GlobalTag


GT = '94X_mc2017_realistic_v17'

process = cms.Process("run")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("JetMETCorrections.Configuration.JetCorrectionServices_cff")
process.load("JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, GT)



process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/mc/RunIIFall17MiniAODv2/TprimeTprime_M-1400_TuneCP5_13TeV-madgraph-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/30000/626C8F3E-4962-E811-8F8A-008CFAE4504C.root'
        )
                            )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100))

#  process.run = cms.EDAnalyzer('BESTAnalyzer',
process.run = cms.EDProducer('BESTProducer_v2',
#                             graphDefinitions = cms.VPSet(
                             name = cms.string('BESTGraph'),
                             path = cms.FileInPath('VLQAnalysis/BESTAnalyzer/data/constantgraph.pb'),
                             means = cms.FileInPath('VLQAnalysis/BESTAnalyzer/data/ScalerParameters.txt'),
#        ),
                             inputJetColl = cms.string('selectedAK8Jets'),
                             isMC = cms.bool(True),
                             isSignal = cms.bool(True),
                             GT_ = cms.string(GT)
                             )

#  process.TFileService = cms.Service("TFileService", fileName = cms.string("BESToutputs.root") )

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string("BESTProducer_v2_out.root"),
                               #  SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
                               #  outputCommands = cms.untracked.vstring('drop *',
                                                                      #  'keep *_fixedGridRhoAll_*_*',
                                                                      #  'keep *_run_*_*',
                                                                      #  )
                               )
process.outpath = cms.EndPath(process.out)

process.p = cms.Path(process.run)
