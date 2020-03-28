#=========================================================================================
# run_JetHT_2017.py ----------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Load Modules and Settings --------------------------------------------------------------
#-----------------------------------------------------------------------------------------

import FWCore.ParameterSet.Config as cms
#from JMEAnalysis.JetToolbox.jetToolbox_cff import jetToolbox
#from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
from Configuration.AlCa.GlobalTag import GlobalTag

# Global Tag
GT = '94X_dataRun2_ReReco_EOY17_v6'

process = cms.Process("run")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("JetMETCorrections.Configuration.JetCorrectionServices_cff")
process.load("JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, GT)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))

#-----------------------------------------------------------------------------------------
# Input Parameters -----------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

# Define the input files
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        '/store/data/Run2017F/JetHT/MINIAOD/31Mar2018-v1/30000/06971FF9-9637-E811-AECE-B496910A9828.root'
        )
                            )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# Inputs to the EDAnalyzer
process.run = cms.EDAnalyzer('BESTAnalyzer',
#                             graphDefinitions = cms.VPSet(
                             name = cms.string('BESTGraph'),
                             path = cms.FileInPath('VLQAnalysis/BESTAnalyzer/data/constantgraph.pb'),
                             means = cms.FileInPath('VLQAnalysis/BESTAnalyzer/data/ScalerParameters.txt'),
#        ),
                             inputJetColl = cms.string('selectedAK8Jets'),
                             isMC = cms.bool(False),
                             isSignal = cms.bool(False),
                             GT_ = cms.string(GT)
                             )

# Name the output files
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

# list steps for running
process.p = cms.Path(process.run)
