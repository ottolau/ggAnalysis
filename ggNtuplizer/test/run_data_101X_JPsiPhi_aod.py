import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask

process = cms.Process('ggKit')

patAlgosToolsTask = getPatAlgosToolsTask(process)

#process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(False) )

#process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load("Configuration.Geometry.GeometryIdeal_cff" )
#process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff" )
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

#process.load("Configuration.Geometry.GeometryRecoDB_cff")
#process.load("Configuration.StandardSequences.MagneticField_cff")

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
#process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_ReReco_EOY17_v2')
#process.GlobalTag = GlobalTag(process.GlobalTag, '101X_dataRun2_HLT_v7')

#process.Tracer = cms.Service("Tracer")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        #'/store/data/Run2018A/DoubleMuon/AOD/22May2018-v1/80000/C694B48F-595E-E811-9634-7CD30AD0A684.root'
        'file:3AAB8D5C-BF59-E811-9475-002590E7DFE4.root'
        )
                            )

#process.load("PhysicsTools.PatAlgos.patSequences_cff")

#process.load( "PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff" )
#process.load( "PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff" )
#process.load( "PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff" )

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


#from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *
runOnData( process,  names=['Photons', 'Electrons','Muons'], outputModules = [] )
#runOnData( process,  names=['Photons', 'Electrons','Muons','Taus','Jets'], outputModules = [] )
#runOnData( process, outputModules = [] )
removeMCMatching(process, names=['All'], outputModules=[])

process.TFileService = cms.Service("TFileService", fileName = cms.string('ggtree_data.root'))

process.load('JetMETCorrections.Configuration.JetCorrectors_cff')

#process.ak4PFchsCorrectedJets   = cms.EDProducer('CorrectedPFJetProducer',
#    src         = cms.InputTag('ak4PFchsJets'),
#    correctors  = cms.VInputTag('ak4PFCHSL1FastL2L3Corrector')
#    )

'''
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import *
runMETCorrectionsAndUncertainties = RunMETCorrectionsAndUncertainties()
runMETCorrectionsAndUncertainties(process, metType="PF",
                                      correctionLevel=["T1"],
                                      computeUncertainties=True,
                                      produceIntermediateCorrections=False,
                                      addToPatDefaultSequence=False,
                                      runOnData=True,
                                      onMiniAOD=False,
                                      reapplyJEC=True,
                                      reclusterJets=False,
                                      jetSelection="pt>15 && abs(eta)<9.9",
                                      recoMetFromPFCs=False,
                                      autoJetCleaning="LepClean",
                                      manualJetConfig=False,
                                      jetFlavor="AK4PFchs",
                                      jetCorLabelUpToL3="ak4PFCHSL1FastL2L3Corrector",
                                      jetCorLabelL3Res="ak4PFCHSL1FastL2L3ResidualCorrector",
                                      jecUncertaintyFile="",
                                      CHS=False,
                                      postfix="",
                                      )
'''
process.load("ggAnalysis.ggNtuplizer.ggNtuplizer_cfi")
#process.load("ggAnalysis.ggNtuplizer.ggPhotonIso_CITK_PUPPI_cff")
process.ggNtuplizer.dumpSoftDrop= cms.bool(False)
process.ggNtuplizer.runHFElectrons=cms.bool(False)
process.ggNtuplizer.isAOD=cms.bool(True)
process.ggNtuplizer.doGenParticles=cms.bool(False)
process.ggNtuplizer.dumpSubJets=cms.bool(False)
process.ggNtuplizer.dumpJets=cms.bool(False)
process.ggNtuplizer.dumpTaus=cms.bool(False)

process.ttk = cms.EDProducer(
    'AddElectronTransientTrack',
    patEleSrc = cms.InputTag('selectedPatElectrons'),
)

process.p = cms.Path(
    ###process.reapplyJEC*
    ###process.pfImpactParameterTagInfosAK8 *
    ###process.pfInclusiveSecondaryVertexFinderTagInfosAK8 *
    ###process.pfBoostedDoubleSecondaryVertexAK8BJetTags *        
    #process.fullPatMetSequence* 
    #process.egcorrMET*
    #process.calibratedPatElectrons*
    #process.calibratedPatPhotons*
    #process.content*
    #process.egmGsfElectronIDSequence*
    #process.egmPhotonIDSequence*
    process.ttk*
    #process.ak4PFchsCorrectedJets*
    process.ggNtuplizer,
    patAlgosToolsTask
)

