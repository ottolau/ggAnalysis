if __name__ == '__main__':

# Usage : python crabConfig.py (to create jobs)
#         ./multicrab -c status -d <work area> (to check job status)

    from CRABAPI.RawCommand import crabCommand
    from httplib import HTTPException

    from CRABClient.UserUtilities import config
    config = config()
    
    from multiprocessing import Process

    # Common configuration

    config.General.workArea     = 'crab_projects_ntuples'
    config.General.transferLogs = False
    config.JobType.pluginName   = 'Analysis' # PrivateMC
    config.JobType.psetName     = 'run_data_101X_BsToJPsiPhi_aod.py'
    config.JobType.inputFiles   = ['Summer16_23Sep2016BCDV4_DATA_L2Relative_AK8PFchs.txt', 'Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK8PFchs.txt', 'Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK8PFchs.txt', 'Summer16_23Sep2016AllV4_DATA.db']
    config.JobType.sendExternalFolder = True
    config.Data.inputDBS        = 'global'    
    config.Data.splitting       = 'LumiBased' # EventBased, FileBased, LumiBased (1 lumi ~= 300 events)
    #config.Data.splitting       = 'Automatic' # EventBased, FileBased, LumiBased (1 lumi ~= 300 events)

    config.Data.totalUnits      = -1
    config.Data.publication     = False
    #config.Data.lumiMask        = 'testjson.txt'
    config.Data.lumiMask        = 'Cert_314472-317080_13TeV_PromptReco_Collisions18_JSON.txt'

    config.Site.storageSite     = 'T3_US_FNALLPC'
    config.JobType.maxMemoryMB = 4000


    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException, hte:
            print hte.headers

    # dataset dependent configuration

    #config.General.requestName = 'DoubleMuon_Run2016B-03Feb2017_ver2-v2_MINIAOD_18May18'
    #config.General.requestName = 'DoubleMuon_Run2016B-07Aug17_ver2-v1_AOD_26May18_GSF'
    config.General.requestName = 'ParkingBPH6_Run2018A-14May2018-v1_AOD_25Jul18_JPsiPhi_4tracks_muTrgObj'
    #config.General.requestName = 'ParkingBPH1_Run2018B-PromptReco-v1_MINIAOD_18Jun18'

    config.Data.unitsPerJob    = 5
    #config.Data.inputDataset   = '/DoubleMuon/Run2016B-23Sep2016-v3/MINIAOD'
    #config.Data.inputDataset   = '/DoubleMuon/Run2016B-03Feb2017_ver2-v2/MINIAOD'
    #config.Data.inputDataset   = '/DoubleMuon/Run2016D-03Feb2017-v1/MINIAOD'
    #config.Data.inputDataset   = '/DoubleMuon/Run2016B-07Aug17_ver2-v1/AOD'
    #config.Data.inputDataset   = '/ParkingBPH1/Run2018B-PromptReco-v1/MINIAOD'
    config.Data.inputDataset   = '/ParkingBPH6/Run2018A-14May2018-v1/AOD'

    config.Data.outLFNDirBase  = '/store/user/klau/'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()




