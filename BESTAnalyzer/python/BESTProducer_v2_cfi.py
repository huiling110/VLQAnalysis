import FWCore.ParameterSet.Config as cms

GT = '94X_mc2017_realistic_v17'
BESTOutput = cms.EDProducer('BESTProducer_v2',
        name = cms.string('BESTGraph'),
        path = cms.FileInPath('VLQAnalysis/BESTAnalyzer/data/constantgraph.pb'),
        means = cms.FileInPath('VLQAnalysis/BESTAnalyzer/data/ScalerParameters.txt'),
        inputJetColl = cms.string('selectedAK8Jets'),
        isMC = cms.bool(True),
        isSignal = cms.bool(True),
        GT_ = cms.string(GT)
        )
