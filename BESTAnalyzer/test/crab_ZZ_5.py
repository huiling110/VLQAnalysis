from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'RadionZZ_4TeV_trees_Analyzer'
config.General.workArea = 'CrabBEST'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run_ZZ.py'
#config.JobType.inputFiles = ['TMVARegression_MLP.weights.xml']
config.JobType.outputFiles = ['BESTInputs.root']
#config.JobType.allowUndistributedCMSSW = True

config.section_("Data")
config.Data.inputDataset = '/RadionToZZ_narrow_M-4000_TuneCP5_13TeV-madgraph/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.ignoreLocality = True
config.Data.publication = False
# This string is used to construct the output dataset name

config.section_("Site")
config.Site.storageSite = 'T3_US_FNALLPC'
config.Site.whitelist = ['T2_US_*']
