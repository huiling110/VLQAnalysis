import FWCore.ParameterSet.Config as cms
from Configuration.AlCa.GlobalTag import GlobalTag


GT = '94X_dataRun2_v11'

process = cms.Process("run")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("JetMETCorrections.Configuration.JetCorrectionServices_cff")
process.load("JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, GT)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/data/Run2017B/JetHT/MINIAOD/31Mar2018-v1/50000/F81C048B-2443-E811-8ADC-0CC47A4D76D0.root'
        )
                            )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


process.run = cms.EDAnalyzer('BESTAnalyzer',
                             name = cms.string('BESTGraph'),
                             path = cms.FileInPath('VLQAnalysis/BESTAnalyzer/data/constantgraph.pb'),
                             means = cms.FileInPath('VLQAnalysis/BESTAnalyzer/data/ScalerParameters.txt'),
                             inputJetColl = cms.string('selectedAK8Jets'),
                             isMC = cms.bool(False),
                             isSignal = cms.bool(False),
                             GT_ = cms.string(GT)
                             )

process.TFileService = cms.Service("TFileService", fileName = cms.string("BESToutputs.root") )

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string("ana_out.root"),
                               SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
                               outputCommands = cms.untracked.vstring('drop *',
								      'keep *_fixedGridRhoAll_*_*',
                                                                      'keep *_run_*_*',
                                                                      #, 'keep *_goodPatJetsCATopTagPF_*_*'
                                                                      #, 'keep recoPFJets_*_*_*'
                                                                      ) 
                               )
process.outpath = cms.EndPath(process.out)

process.p = cms.Path(process.run)
