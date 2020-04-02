from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'JetHT_2017F'
config.General.workArea = 'CrabBEST'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run_Data2017.py'
config.JobType.inputFiles = ['/uscms_data/d3/rband/TT_Analysis/CMSSW_10_2_18/src/VLQAnalysis/BESTAnalyzer/data/constantgraph.pb', '/uscms_data/d3/rband/TT_Analysis/CMSSW_10_2_18/src/VLQAnalysis/BESTAnalyzer/data/ScalerParameters.txt', '/uscms_data/d3/rband/TT_Analysis/CMSSW_10_2_18/src/VLQAnalysis/BESTAnalyzer/data/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt']
config.JobType.outputFiles = ['BESToutputs.root']
#config.JobType.allowUndistributedCMSSW = True

config.section_("Data")
config.Data.inputDataset = '/JetHT/Run2017F-31Mar2018-v1/MINIAOD'
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 20
config.Data.ignoreLocality = True
config.Data.publication = False
# This string is used to construct the output dataset name

config.section_("Site")
config.Site.storageSite = 'T3_US_FNALLPC'
config.Site.whitelist = ['T2_US_*']
