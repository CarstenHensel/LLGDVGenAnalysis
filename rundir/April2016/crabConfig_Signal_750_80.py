from WMCore.Configuration import Configuration
config = Configuration()
config.section_("General")
config.General.requestName = 'TruthAnalysis_750_80'
config.General.workArea = 'truthAnalysis'
config.General.transferOutputs = True
config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'truthAnalysis_26042016.py'
config.JobType.outputFiles = ['output.root']
config.section_("Data")
config.Data.inputDataset = '/CRAB_PrivateMC/mhamer-LLG_750_80_Fall15_stepRECO-413386f6eddb08329706f28eff10fb19/USER'
config.Data.inputDBS='phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
NJOBS = 250 # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.publication = False
config.Data.outputDatasetTag = 'TruthAnalysis_Signal_750_80_Apr2016' 
config.section_("Site")
config.Site.storageSite = 'T2_DE_DESY'
