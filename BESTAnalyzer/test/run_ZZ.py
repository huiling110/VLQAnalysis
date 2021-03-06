import FWCore.ParameterSet.Config as cms
#from JMEAnalysis.JetToolbox.jetToolbox_cff import jetToolbox
#from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
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
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        '/store/mc/RunIIFall17MiniAODv2/RadionToZZ_narrow_M-5000_TuneCP5_13TeV-madgraph/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/240000/88C32EBF-A689-E911-BE4D-A4BF0112BCD4.root'
#        'file://88C32EBF-A689-E911-BE4D-A4BF0112BCD4.root'
        )
                            )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


process.run = cms.EDAnalyzer('BESTAnalyzer',
#                             graphDefinitions = cms.VPSet(
                             name = cms.string('BESTGraph'),
                             path = cms.FileInPath('VLQAnalysis/BESTAnalyzer/test/constantgraph.pb'),
                             means = cms.FileInPath('VLQAnalysis/BESTAnalyzer/test/ScalerParameters.txt'),
#        ),
                             inputJetColl = cms.string('selectedAK8Jets'),
                             isMC_ = cms.bool(True),
                             GT_ = cms.string(GT)
                             )

process.TFileService = cms.Service("TFileService", fileName = cms.string("BESTInputs.root") )

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
