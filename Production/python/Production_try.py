# Produce EventTuple for all channels.
# This file is part of https://github.com/hh-italian-group/h-tautau.

import re
import importlib
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
options.register('sampleType', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                        "Indicates the sample type: Spring15MC, Run2015B, Run2015C, Run2015D")
options.register('applyTriggerMatch', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                        "Apply trigger matching for signal objects. Default: True")
options.register('applyTriggerMatchCut', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                        "Apply trigger matching Cut not for signal samples. Default: True")
options.register('fileList', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                        "List of root files to process.")
options.register('fileNamePrefix', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                        "Prefix to add to input file names.")
options.register('anaChannels', 'all', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                        "Analysis channels to run.")
options.register('tupleOutput', 'eventTuple.root', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                        "Event tuple file.")
options.register('runSVfit', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                        "Run SVfit algorithm on the selected tau pair.")
options.register('runKinFit', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                        "Run HHKinFit algorithm for on the selected tau pair and all possible jet combinations.")
options.register('applyTriggerCut', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                        "Apply trigger cut for signal objects. Default: True")
options.register('storeLHEinfo', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                        "Store LHE information. Default: False")
options.register('applyRecoilCorr', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                        "Apply Met Recoil Corrections")
options.register('nJetsRecoilCorr', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                        "Number of Additional Jets for Recoil Correction")
options.register('lumiFile', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                        "JSON file with lumi mask.")
options.register('eventList', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                        "List of events to process.")
options.register('saveGenTopInfo', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                        "Save generator-level information for top quarks.")
options.register('saveGenBosonInfo', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                        "Save generator-level information for bosons.")
options.register('saveGenJetInfo', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                        "Save generator-level information for jets.")
options.register('saveGenParticleInfo', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                        "Save generator-level information for particles.")
options.register('isEmbedded', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                        "Is DY embedded sample.")
options.register('dumpPython', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                        "Dump full config into stdout.")
options.register('numberOfThreads', 1, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                        "Number of threads.")

options.parseArguments()

sampleConfig = importlib.import_module('h-tautau.Production.sampleConfig')
isData = sampleConfig.IsData(options.sampleType)
period = sampleConfig.GetPeriod(options.sampleType)
triggerCfg = sampleConfig.GetTriggerCfg(options.sampleType)

processName = 'tupleProduction'
process = cms.Process(processName)
process.options = cms.untracked.PSet()
process.options.wantSummary = cms.untracked.bool(False)
process.options.allowUnscheduled = cms.untracked.bool(True)
process.options.numberOfThreads = cms.untracked.uint32(options.numberOfThreads)
process.options.numberOfStreams=cms.untracked.uint32(0)

process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.GlobalTag.globaltag = sampleConfig.GetGlobalTag(options.sampleType)
#from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, options.globalTag, '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')
process.source = cms.Source('PoolSource', fileNames = cms.untracked.vstring())
process.TFileService = cms.Service('TFileService', fileName = cms.string(options.tupleOutput) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(0) )

if len(options.fileList) > 0:
    from AnalysisTools.Run.readFileList import *
    readFileList(process.source.fileNames, options.fileList, options.fileNamePrefix)
    process.maxEvents.input = options.maxEvents

if len(options.lumiFile) > 0:
    import FWCore.PythonUtilities.LumiList as LumiList
    process.source.lumisToProcess = LumiList.LumiList(filename = options.lumiFile).getVLuminosityBlockRange()

if options.eventList != '':
    process.source.eventsToProcess = cms.untracked.VEventRange(re.split(',', options.eventList))

## and add btag discriminator to the event content
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

btagVector = [
  'pfDeepFlavourJetTags:probb',
  'pfDeepFlavourJetTags:probbb',
  'pfDeepFlavourJetTags:problepb',
  'pfDeepFlavourJetTags:probc',
  'pfDeepFlavourJetTags:probuds',
  'pfDeepFlavourJetTags:probg'
   ]

if period == 'Run2016':
    btagVector2016 = [
        'pfDeepCSVJetTags:probudsg',
        'pfDeepCSVJetTags:probb',
        'pfDeepCSVJetTags:probc',
        'pfDeepCSVJetTags:probbb'
    ]
    btagVector.append(btagVector2016)

jec_levels = ['L1FastJet', 'L2Relative', 'L3Absolute']

if(not options.isEmbedded):
    jec_levels.append('L2L3Residual')

updateJetCollection(
    process,
    jetSource = cms.InputTag('slimmedJets'),
    pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
    svSource = cms.InputTag('slimmedSecondaryVertices'),
    jetCorrections = ('AK4PFchs', cms.vstring(jec_levels), 'None'),
    btagDiscriminators = btagVector,
    postfix='NewDFTraining'
)
process.jecSequence = cms.Sequence(process.patJetCorrFactorsNewDFTraining *
                                   process.updatedPatJetsNewDFTraining *
                                   process.patJetCorrFactorsTransientCorrectedNewDFTraining *
                                   process.pfImpactParameterTagInfosNewDFTraining *
                                   process.pfInclusiveSecondaryVertexFinderTagInfosNewDFTraining *
                                   process.pfDeepCSVTagInfosNewDFTraining *
                                   process.pfDeepFlavourTagInfosNewDFTraining *
                                   process.pfDeepFlavourJetTagsNewDFTraining *
                                   process.updatedPatJetsTransientCorrectedNewDFTraining *
                                   process.selectedUpdatedPatJetsNewDFTraining)

process.load('RecoMET.METFilters.ecalBadCalibFilter_cfi')

baddetEcallist = cms.vuint32(
    [872439604,872422825,872420274,872423218,
     872423215,872416066,872435036,872439336,
     872420273,872436907,872420147,872439731,
     872436657,872420397,872439732,872439339,
     872439603,872422436,872439861,872437051,
     872437052,872420649,872422436,872421950,
     872437185,872422564,872421566,872421695,
     872421955,872421567,872437184,872421951,
     872421694,872437056,872437057,872437313])


process.ecalBadCalibReducedMINIAODFilter = cms.EDFilter(
    "EcalBadCalibFilter",
    EcalRecHitSource = cms.InputTag("reducedEgamma:reducedEERecHits"),
    ecalMinEt        = cms.double(50.),
    baddetEcal    = baddetEcallist,
    taggingMode = cms.bool(True),
    debug = cms.bool(False)
    )

from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

runMetCorAndUncFromMiniAOD (
        process,
        isData = isData, # false for MC
        fixEE2017 = True,
        fixEE2017Params = {'userawPt': True, 'ptThreshold':50.0, 'minEtaThreshold':2.65, 'maxEtaThreshold': 3.139} ,
        postfix = "ModifiedMET"
)

MetInputTag = cms.InputTag('slimmedMETsModifiedMET', '', processName)

# Update electron ID following recommendations from https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
ele_era = { 'Run2016': '2016-Legacy', 'Run2017': '2017-Nov17ReReco'} #add 2018 'Run2018': '2018-Prompt'
setupEgammaPostRecoSeq(process, runVID=True, runEnergyCorrections=False, era=ele_era[period])

if options.anaChannels == 'all':
    channels = [ 'eTau', 'muTau', 'tauTau', 'muMu' ]
else:
    channels = re.split(',', options.anaChannels)


for channel in channels:
    producerName = 'tupleProducer_{}'.format(channel)
    producerClassName = 'TupleProducer_{}'.format(channel)

    setattr(process, producerName, cms.EDAnalyzer(producerClassName,
        electronSrc             = cms.InputTag('slimmedElectrons'),
        tauSrc                  = cms.InputTag('slimmedTaus'),
        muonSrc                 = cms.InputTag('slimmedMuons'),
        vtxSrc                  = cms.InputTag('offlineSlimmedPrimaryVertices'),
        jetSrc                  = cms.InputTag('selectedUpdatedPatJetsNewDFTraining'),
        fatJetSrc               = cms.InputTag('slimmedJetsAK8'),
        PUInfo                  = cms.InputTag('slimmedAddPileupInfo'),
        pfMETSrc                = MetInputTag,
        genMetSrc               = cms.InputTag('genMetTrue'),
        prescales               = cms.InputTag('patTrigger'),
        objects                 = cms.InputTag('slimmedPatTrigger'),
        lheEventProducts        = cms.InputTag('externalLHEProducer'),
        genEventInfoProduct     = cms.InputTag('generator'),
        topGenEvent             = cms.InputTag('genEvt'),
        genParticles            = cms.InputTag('prunedGenParticles'),
        genJets                 = cms.InputTag('slimmedGenJets'),
        l1JetParticleProduct    = cms.InputTag('l1extraParticles', 'IsoTau'),
        isMC                    = cms.bool(not isData),
        applyTriggerMatch       = cms.bool(options.applyTriggerMatch),
        applyTriggerMatchCut    = cms.bool(options.applyTriggerMatchCut),
        runSVfit                = cms.bool(options.runSVfit),
        runKinFit               = cms.bool(options.runKinFit),
        applyTriggerCut         = cms.bool(options.applyTriggerCut),
        storeLHEinfo            = cms.bool(options.storeLHEinfo),
        applyRecoilCorr         = cms.bool(options.applyRecoilCorr),
        nJetsRecoilCorr         = cms.int32(options.nJetsRecoilCorr),
        period                  = cms.string(period),
        triggerCfg              = cms.string(triggerCfg),
        saveGenTopInfo          = cms.bool(options.saveGenTopInfo),
        saveGenBosonInfo        = cms.bool(options.saveGenBosonInfo),
        saveGenJetInfo          = cms.bool(options.saveGenJetInfo),
        saveGenParticleInfo     = cms.bool(options.saveGenParticleInfo),
        isEmbedded              = cms.bool(options.isEmbedded),
        rho                     = cms.InputTag('fixedGridRhoAll'),
    ))


if period == 'Run2017':
    process.p = cms.Path(
        process.egmGsfElectronIDSequence *
        process.egammaPostRecoSeq *
        process.jecSequence *
        #process.rerunMvaIsolationSequence *
        #getattr(process, updatedTauName) *
        #process.fullPatMetSequence +
        process.fullPatMetSequenceModifiedMET
    )

if options.dumpPython:
    print process.dumpPython()