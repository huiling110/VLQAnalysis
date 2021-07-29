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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/mc/RunIIFall17MiniAODv2/QCD_Pt-15to7000_TuneCP5_Flat2017_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/AA231798-AA43-E811-BCDC-0CC47A7C35A8.root'
        )
                            )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


process.run = cms.EDAnalyzer('BESTAnalyzer',
                             name = cms.string('BESTGraph'),
                             path = cms.FileInPath('VLQAnalysis/BESTAnalyzer/test/constantgraph.pb'),
                             means = cms.FileInPath('VLQAnalysis/BESTAnalyzer/test/ScalerParameters.txt'),
                             inputJetColl = cms.string('selectedAK8Jets'),
                             isMC = cms.bool(True),
                             isSignal = cms.bool(False),
                             GT_ = cms.string(GT)
                             )

process.TFileService = cms.Service("TFileService", fileName = cms.string("BESTInputs.root") )
#???what is the role of a cms.Service ?

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
