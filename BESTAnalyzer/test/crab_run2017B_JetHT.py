from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'run2017B_JetHT'
config.General.workArea = 'CrabBEST'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run_JetHT_2017.py'
#config.JobType.inputFiles = ['constantgraph.pb', 'ScalerParameters.txt']
config.JobType.outputFiles = ['BESToutputs.root']
#config.JobType.allowUndistributedCMSSW = True

config.section_("Data")
config.Data.inputDataset = '/JetHT/Run2017B-31Mar2018-v1/MINIAOD'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.ignoreLocality = True
config.Data.publication = False
# This string is used to construct the output dataset name

config.section_("Site")
config.Site.storageSite = 'T3_US_FNALLPC'
config.Site.whitelist = ['T2_US_*']
