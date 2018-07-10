import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask

process = cms.Process('ggKit')

patAlgosToolsTask = getPatAlgosToolsTask(process)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

process.load('Configuration.StandardSequences.Services_cff')
#process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_2016SeptRepro_v7')
#process.GlobalTag = GlobalTag(process.GlobalTag, '101X_dataRun2_HLT_v7')

#process.Tracer = cms.Service("Tracer")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#jec from sqlite
#process.load("CondCore.DBCommon.CondDBCommon_cfi")
#from CondCore.DBCommon.CondDBSetup_cfi import *
#process.jec = cms.ESSource("PoolDBESSource",CondDBSetup,
# connect = cms.string('sqlite:Summer16_23Sep2016AllV4_DATA.db'),
# toGet = cms.VPSet(
# cms.PSet(
#  record = cms.string('JetCorrectionsRecord'),
#  tag = cms.string('JetCorrectorParametersCollection_Summer16_23Sep2016AllV4_DATA_AK4PFchs'),
#  label = cms.untracked.string('AK4PFchs')
# ),
# cms.PSet(
#  record = cms.string('JetCorrectionsRecord'),
#  tag = cms.string('JetCorrectorParametersCollection_Summer16_23Sep2016AllV4_DATA_AK8PFchs'),
#  label = cms.untracked.string('AK8PFchs')
#  )))
#process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        #'/store/data/Run2018A/ParkingBPH3/AOD/14May2018-v1/70000/3AAB8D5C-BF59-E811-9475-002590E7DFE4.root'
        #'/store/data/Run2018A/ParkingBPH3/AOD/14May2018-v1/70000/D82B8CFE-D15A-E811-96DA-0CC47AF9B2FE.root'
        #'/store/data/Run2018A/DoubleMuon/AOD/22May2018-v1/80000/C694B48F-595E-E811-9634-7CD30AD0A684.root'
        'file:3AAB8D5C-BF59-E811-9475-002590E7DFE4.root'
        )
                            )

#process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.load('PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff')
patAlgosToolsTask.add(process.patCandidatesTask)
#process.load( "PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff" )
#patAlgosToolsTask.add(process.triggerProducerTask)
process.load('PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff')
patAlgosToolsTask.add(process.selectedPatCandidatesTask)
process.load('PhysicsTools.PatAlgos.cleaningLayer1.cleanPatCandidates_cff')
patAlgosToolsTask.add(process.cleanPatCandidatesTask)

### add trigger information to the configuration
from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger( process, None, None, None, None, '' )

from PhysicsTools.PatAlgos.tools.coreTools import *
runOnData( process,  names=['Photons', 'Electrons','Muons','Taus','Jets'], outputModules = [] )
#runOnData( process, outputModules = [] )
removeMCMatching(process, names=['All'], outputModules=[])

# this loads all available b-taggers
#process.load("RecoBTag.Configuration.RecoBTag_cff")
#process.load("RecoBTag.SecondaryVertex.pfBoostedDoubleSecondaryVertexAK8BJetTags_cfi")
#process.pfImpactParameterTagInfosAK8.primaryVertex = cms.InputTag("offlineSlimmedPrimaryVertices")
#process.pfImpactParameterTagInfosAK8.candidates = cms.InputTag("packedPFCandidates")
#process.pfImpactParameterTagInfosAK8.jets = cms.InputTag("slimmedJetsAK8")
#process.load("RecoBTag.SecondaryVertex.pfInclusiveSecondaryVertexFinderTagInfosAK8_cfi")
#process.pfInclusiveSecondaryVertexFinderTagInfosAK8.extSVCollection = cms.InputTag("slimmedSecondaryVertices")

process.TFileService = cms.Service("TFileService", fileName = cms.string('ggtree_data.root'))


useAOD = True
#####VID framework####################
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
if useAOD == True :
    dataFormat = DataFormat.AOD
    process.load("ggAnalysis.ggNtuplizer.ggNtuplizer_cfi")
    #process.load("ggAnalysis.ggNtuplizer.ggMETFilters_cff")
    process.ggNtuplizer.addFilterInfo=cms.bool(False)
    process.ggNtuplizer.development=cms.bool(True)
    #from JMEAnalysis.JetToolbox.jetToolbox_cff import *
    #jetToolbox( process, 'ak4', 'ak4PFJetsCHS', 'out', runOnMC = False, miniAOD= False, addSoftDrop=True, addSoftDropSubjets=True, addNsub=True, addPUJetID=True, JETCorrPayload='AK4PFchs', JETCorrLevels=['L1FastJet','L2Relative', 'L3Absolute','L2L3Residual'] )
    #jetToolbox( process, 'ak8', 'ak8PFJetsCHS', 'out', runOnMC = False, miniAOD= False, addSoftDrop=True, addSoftDropSubjets=True, addNsub=True, bTagDiscriminators=['pfBoostedDoubleSecondaryVertexAK8BJetTags'] )
    process.ggNtuplizer.dumpSoftDrop= cms.bool(False)

else :
    dataFormat = DataFormat.MiniAOD
    process.load("ggAnalysis.ggNtuplizer.ggNtuplizer_miniAOD_cfi")
    process.load("CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi")
    process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
    process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False) 
    process.HBHENoiseFilterResultProducer.defaultDecision = cms.string("HBHENoiseFilterResultRun2Loose")
    process.ggNtuplizer.dumpSoftDrop= cms.bool(True)
    process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")
    from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetCorrFactorsUpdated
    process.patJetCorrFactorsReapplyJEC = process.patJetCorrFactorsUpdated.clone(
        src = cms.InputTag("slimmedJets"),
        levels = ['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual'],
        payload = 'AK4PFchs' ) # Make sure to choose the appropriate levels and payload here!

    from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetsUpdated
    process.patJetsReapplyJEC = process.patJetsUpdated.clone(
        jetSource = cms.InputTag("slimmedJets"),
        jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
        )

    process.patJetAK8CorrFactorsReapplyJEC = process.patJetCorrFactorsUpdated.clone(
        src = cms.InputTag("slimmedJetsAK8"),
        levels = ['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual'],
        payload = 'AK8PFchs' ) # Make sure to choose the appropriate levels and payload here!

    process.patJetsAK8ReapplyJEC = process.patJetsUpdated.clone(
        jetSource = cms.InputTag("slimmedJetsAK8"),
        jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetAK8CorrFactorsReapplyJEC"))
        )

    process.reapplyJEC = cms.Sequence( process.patJetCorrFactorsReapplyJEC + process. patJetsReapplyJEC +process.patJetAK8CorrFactorsReapplyJEC + process. patJetsAK8ReapplyJEC )

#process.ggNtuplizer.jecAK8PayloadNames=cms.vstring(jecLevels)
process.ggNtuplizer.runHFElectrons=cms.bool(False)
process.ggNtuplizer.isAOD=cms.bool(useAOD)
process.ggNtuplizer.doGenParticles=cms.bool(False)
process.ggNtuplizer.dumpSubJets=cms.bool(False)
process.ggNtuplizer.dumpJets=cms.bool(False)
process.ggNtuplizer.dumpTaus=cms.bool(False)
#process.ggNtuplizer.pfMETLabel=cms.InputTag("slimmedMETsMuEGClean", "", "ggKit")
process.ggNtuplizer.electronSrc=cms.InputTag("selectedPatElectrons")

switchOnVIDElectronIdProducer(process, dataFormat)
#switchOnVIDPhotonIdProducer(process, dataFormat)

# define which IDs we want to produce
#my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff',
#                 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronHLTPreselecition_Summer16_V1_cff',
#                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff',
#                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff',
#                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_HZZ_V1_cff']

my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V1_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
                                                                
my_phoid_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring16_V2p2_cff',
                    'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring16_nonTrig_V1_cff']

process.load("RecoEgamma.ElectronIdentification.ElectronIDValueMapProducer_cfi")
process.electronIDValueMapProducer.srcAOD = cms.InputTag('selectedPatElectrons')
process.electronMVAValueMapProducer.srcAOD = cms.InputTag('selectedPatElectrons')
#process.electronIDValueMapProducer.srcMiniAOD = cms.InputTag('slimmedElectrons')
#process.electronMVAValueMapProducer.srcMiniAOD = cms.InputTag('slimmedElectrons')

#process.photonIDValueMapProducer.srcAOD = cms.InputTag('selectedPatPhotons')
#process.photonMVAValueMapProducer.srcAOD = cms.InputTag('selectedPatPhotons')

process.ttk = cms.EDProducer(
    'AddElectronTransientTrack',
    patEleSrc = cms.InputTag('selectedPatElectrons'),

)


#add them to the VID producer
#for idmod in my_phoid_modules:
  #setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

if useAOD == True:
  process.p = cms.Path(
      process.egmGsfElectronIDSequence*
      #process.egmPhotonIDSequence*
      #process.ggMETFiltersSequence*
      process.ttk*
      process.ggNtuplizer,
      patAlgosToolsTask
    )
else:
  process.p = cms.Path(
      process.HBHENoiseFilterResultProducer* # produces HBHE bools
#        process.ApplyBaselineHBHENoiseFilter*  # reject events 
      process.reapplyJEC*
      process.pfImpactParameterTagInfosAK8 *
      process.pfInclusiveSecondaryVertexFinderTagInfosAK8 *
      process.pfBoostedDoubleSecondaryVertexAK8BJetTags *
      process.calibratedPatElectrons*
      process.calibratedPatPhotons*
      process.egmGsfElectronIDSequence*
      process.egmPhotonIDSequence*
      process.ggNtuplizer
      )
  


