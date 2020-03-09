from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'Tprime_1400GeV'
config.General.workArea = 'CrabBEST'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run_Tprime.py'
config.JobType.inputFiles = ['constantgraph.pb', 'ScalerParameters.txt']
config.JobType.outputFiles = ['BESToutputs.root']
#config.JobType.allowUndistributedCMSSW = True

config.section_("Data")
config.Data.inputDataset = '/TprimeTprime_M-1400_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.ignoreLocality = True
config.Data.publication = False
# This string is used to construct the output dataset name

config.section_("Site")
config.Site.storageSite = 'T3_US_FNALLPC'
config.Site.whitelist = ['T2_US_*']
