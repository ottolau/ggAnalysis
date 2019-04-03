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
    config.JobType.psetName     = 'run_data_105X_lowPtElectrons_aod.py'
    #config.JobType.psetName     = 'run_data_101X_BsToJPsiPhi_mini.py'
    config.JobType.inputFiles   = ['Summer16_23Sep2016BCDV4_DATA_L2Relative_AK8PFchs.txt', 'Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK8PFchs.txt', 'Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK8PFchs.txt', 'Summer16_23Sep2016AllV4_DATA.db']
    config.JobType.sendExternalFolder = True
    config.Data.inputDBS        = 'global'    
    config.Data.splitting       = 'LumiBased' # EventBased, FileBased, LumiBased (1 lumi ~= 300 events)
    #config.Data.splitting       = 'Automatic' # EventBased, FileBased, LumiBased (1 lumi ~= 300 events)

    config.Data.totalUnits      = -1
    config.Data.publication     = False
    #config.Data.lumiMask        = 'testjson.txt'
    #config.Data.lumiMask        = 'Cert_314472-317080_13TeV_PromptReco_Collisions18_JSON.txt'
    config.Data.lumiMask        = 'Cert_314472-323523_13TeV_PromptReco_Collisions18_JSON.txt'

    config.Site.storageSite     = 'T3_US_FNALLPC'
    config.JobType.maxMemoryMB = 4000
    config.Data.allowNonValidInputDataset = True

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException, hte:
            print hte.headers

    # dataset dependent configuration

    config.General.requestName = 'ParkingBPH1_Run2018A-22Mar2019-v1_AOD_02Apr19_BsPhiLL_lowPtElectrons'

    config.Data.unitsPerJob    = 5
    config.Data.inputDataset   = '/ParkingBPH1/Run2018A-22Mar2019-v1/AOD'

    config.Data.outLFNDirBase  = '/store/user/klau/LowPtElectrons'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()




