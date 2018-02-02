# Configurations dependent on the sample type.
# This file is part of https://github.com/hh-italian-group/h-tautau.

import sys
from sets import Set
import FWCore.ParameterSet.Config as cms

mcSampleTypes = Set([ 'Summer16MC', 'Fall17MC' ])
dataSampleTypes = Set([ 'Run2016' , 'Run2017' ])

hltPaths_eTau_Run2016 = cms.VPSet(
    cms.PSet( pattern = cms.string("HLT_Ele23_WPLoose_Gsf_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltEle23WPLooseGsfTrackIsoFilter") ),
    cms.PSet( pattern = cms.string("HLT_Ele24_eta2p1_WPLoose_Gsf_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltSingleEle24WPLooseGsfTrackIsoFilter") ),
    cms.PSet( pattern = cms.string("HLT_Ele25_WPTight_Gsf_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltEle25WPTightGsfTrackIsoFilter") ),
    cms.PSet( pattern = cms.string("HLT_Ele25_eta2p1_WPLoose_Gsf_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltEle25erWPLooseGsfTrackIsoFilter") ),
    cms.PSet( pattern = cms.string("HLT_Ele25_eta2p1_WPTight_Gsf_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltEle25erWPTightGsfTrackIsoFilter") ),
    cms.PSet( pattern = cms.string("HLT_Ele27_WPLoose_Gsf_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltEle27noerWPLooseGsfTrackIsoFilter") ),
    cms.PSet( pattern = cms.string("HLT_Ele27_WPTight_Gsf_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltEle27WPTightGsfTrackIsoFilter") ),
    cms.PSet( pattern = cms.string("HLT_Ele27_eta2p1_WPLoose_Gsf_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltEle27erWPLooseGsfTrackIsoFilter") ),
    cms.PSet( pattern = cms.string("HLT_Ele27_eta2p1_WPTight_Gsf_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltEle27erWPTightGsfTrackIsoFilter") ),
    cms.PSet( pattern = cms.string("HLT_Ele32_eta2p1_WPTight_Gsf_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltEle32WPTightGsfTrackIsoFilter") ),
    cms.PSet( pattern = cms.string("HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltEle22WPLooseL1SingleIsoEG20erGsfTrackIsoFilter",
                                               "hltOverlapFilterSingleIsoEle22WPLooseGsfLooseIsoPFTau20"),
              filters2 = cms.untracked.vstring("hltPFTau20TrackLooseIso",
                                               "hltOverlapFilterSingleIsoEle22WPLooseGsfLooseIsoPFTau20") ),
    cms.PSet( pattern = cms.string("HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltEle24WPLooseL1SingleIsoEG22erGsfTrackIsoFilter",
                                               "hltOverlapFilterSingleIsoEle24WPLooseGsfLooseIsoPFTau20"),
              filters2 = cms.untracked.vstring("hltPFTau20TrackLooseIso",
                                               "hltOverlapFilterSingleIsoEle24WPLooseGsfLooseIsoPFTau20") ),
    cms.PSet( pattern = cms.string("HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltEle24WPLooseL1IsoEG22erTau20erGsfTrackIsoFilter",
                                               "hltOverlapFilterIsoEle24WPLooseGsfLooseIsoPFTau20"),
              filters2 = cms.untracked.vstring("hltPFTau20TrackLooseIso",
                                               "hltOverlapFilterIsoEle24WPLooseGsfLooseIsoPFTau20") ),
    cms.PSet( pattern = cms.string("HLT_Ele27_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltEle27erWPLooseGsfTrackIsoFilter",
                                               "hltOverlapFilterIsoEle27WPLooseGsfLooseIsoPFTau20"),
              filters2 = cms.untracked.vstring("hltPFTau20TrackLooseIso",
                                               "hltOverlapFilterIsoEle27WPLooseGsfLooseIsoPFTau20") ),
)

hltPaths_eTau_Run2017 = cms.VPSet(
    cms.PSet( pattern = cms.string("HLT_Ele32_WPTight_Gsf_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltEle32WPTightGsfTrackIsoFilter") ),
    cms.PSet( pattern = cms.string("HLT_Ele35_WPTight_Gsf_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltEle35noerWPTightGsfTrackIsoFilter") ),
    cms.PSet( pattern = cms.string("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltEle24erWPTightGsfTrackIsoFilterForTau",
                                               "hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30"),
              filters2 = cms.untracked.vstring("hltSelectedPFTau30LooseChargedIsolationL1HLTMatched",
                                              "hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30") ),
)

hltPaths_eTau = {
    'Run2016'   : hltPaths_eTau_Run2016,
    'Run2017'   : hltPaths_eTau_Run2017
}

hltPaths_muTau_Run2016 = cms.VPSet(
    cms.PSet( pattern = cms.string("HLT_IsoMu18_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09") ),
    cms.PSet( pattern = cms.string("HLT_IsoMu20_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sMu18L1f0L2f10QL3f20QL3trkIsoFiltered0p09") ),
    cms.PSet( pattern = cms.string("HLT_IsoMu22_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09") ),
    cms.PSet( pattern = cms.string("HLT_IsoMu22_eta2p1_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sSingleMu20erL1f0L2f10QL3f22QL3trkIsoFiltered0p09") ),
    cms.PSet( pattern = cms.string("HLT_IsoMu24_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09") ),
    cms.PSet( pattern = cms.string("HLT_IsoMu27_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p09") ),
    cms.PSet( pattern = cms.string("HLT_IsoTkMu18_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltL3fL1sMu16L1f0Tkf18QL3trkIsoFiltered0p09") ),
    cms.PSet( pattern = cms.string("HLT_IsoTkMu20_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltL3fL1sMu18L1f0Tkf20QL3trkIsoFiltered0p09") ),
    cms.PSet( pattern = cms.string("HLT_IsoTkMu22_eta2p1_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltL3fL1sMu20erL1f0Tkf22QL3trkIsoFiltered0p09") ),
    cms.PSet( pattern = cms.string("HLT_IsoTkMu22_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltL3fL1sMu20L1f0Tkf22QL3trkIsoFiltered0p09") ),
    cms.PSet( pattern = cms.string("HLT_IsoTkMu24_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09") ),
    cms.PSet( pattern = cms.string("HLT_IsoTkMu27_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltL3fL1sMu22Or25L1f0Tkf27QL3trkIsoFiltered0p09") ),
    cms.PSet( pattern = cms.string("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sSingleMu16erL1f0L2f10QL3f17QL3trkIsoFiltered0p09",
                                               "hltOverlapFilterSingleIsoMu17LooseIsoPFTau20"),
              filters2 = cms.untracked.vstring("hltPFTau20TrackLooseIsoAgainstMuon",
                                               "hltOverlapFilterSingleIsoMu17LooseIsoPFTau20") ),
    cms.PSet( pattern = cms.string("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sMu16erTauJet20erL1f0L2f10QL3f17QL3trkIsoFiltered0p09",
                                               "hltOverlapFilterIsoMu17LooseIsoPFTau20"),
              filters2 = cms.untracked.vstring("hltPFTau20TrackLooseIsoAgainstMuon",
                                               "hltOverlapFilterIsoMu17LooseIsoPFTau20") ),
    cms.PSet( pattern = cms.string("HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sSingleMu18erIorSingleMu20erL1f0L2f10QL3f19QL3trkIsoFiltered0p09",
                                               "hltOverlapFilterSingleIsoMu19LooseIsoPFTau20"),
              filters2 = cms.untracked.vstring("hltPFTau20TrackLooseIsoAgainstMuon",
                                               "hltOverlapFilterSingleIsoMu19LooseIsoPFTau20") ),
    cms.PSet( pattern = cms.string("HLT_IsoMu19_eta2p1_LooseIsoPFTau20_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sMu18erTauJet20erL1f0L2f10QL3f19QL3trkIsoFiltered0p09",
                                               "hltOverlapFilterIsoMu19LooseIsoPFTau20"),
              filters2 = cms.untracked.vstring("hltPFTau20TrackLooseIsoAgainstMuon",
                                               "hltOverlapFilterIsoMu19LooseIsoPFTau20") ),
    cms.PSet( pattern = cms.string("HLT_IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sSingleMu20erIorSingleMu22erL1f0L2f10QL3f21QL3trkIsoFiltered0p09",
                                               "hltOverlapFilterSingleIsoMu21LooseIsoPFTau20"),
              filters2 = cms.untracked.vstring("hltPFTau20TrackLooseIsoAgainstMuon",
                                               "hltOverlapFilterSingleIsoMu21LooseIsoPFTau20") ),
)

hltPaths_muTau_Run2017 = cms.VPSet(
    cms.PSet( pattern = cms.string("HLT_IsoMu24_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07") ),
    cms.PSet( pattern = cms.string("HLT_IsoMu27_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07") ),
    cms.PSet( pattern = cms.string("HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_eta2p1_SingleL1_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltPFTau20TrackLooseChargedIsoAgainstMuon",
                                               "hltOverlapFilterIsoMu24LooseChargedIsoPFTau20"),
              filters2 = cms.untracked.vstring("hltL3crIsoL1sSingleMu22erL1f0L2f10QL3f24QL3trkIsoFiltered0p07",
                                               "hltOverlapFilterIsoMu24LooseChargedIsoPFTau20") ),
    cms.PSet( pattern = cms.string("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07",
                                               "hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded"),
              filters2 = cms.untracked.vstring("hltSelectedPFTau27LooseChargedIsolationAgainstMuonL1HLTMatched",
                                               "hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded") ),
)


hltPaths_muTau = {
    'Run2016'   : hltPaths_muTau_Run2016,
    'Run2017'   : hltPaths_muTau_Run2017
}

hltPaths_tauTau_Run2016 = cms.VPSet(
    cms.PSet( pattern = cms.string("HLT_DoubleMediumIsoPFTau32_Trk1_eta2p1_Reg_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltDoublePFTau32TrackPt1MediumIsolationDz02Reg"),
              filters2 = cms.untracked.vstring("hltDoublePFTau32TrackPt1MediumIsolationDz02Reg") ),
    cms.PSet( pattern = cms.string("HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltDoublePFTau35TrackPt1MediumIsolationDz02Reg"),
              filters2 = cms.untracked.vstring("hltDoublePFTau35TrackPt1MediumIsolationDz02Reg") ),
    cms.PSet( pattern = cms.string("HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltDoublePFTau40TrackPt1MediumIsolationDz02Reg"),
              filters2 = cms.untracked.vstring("hltDoublePFTau40TrackPt1MediumIsolationDz02Reg") ),
    cms.PSet( pattern = cms.string("HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v"),
              nLegs = cms.untracked.uint32(2)),
    cms.PSet( pattern = cms.string("HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg_v"),
              nLegs = cms.untracked.uint32(2)),
    cms.PSet( pattern = cms.string("HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_v"),
              nLegs = cms.untracked.uint32(2)),
)

hltPaths_tauTau_Run2017 = cms.VPSet(
    cms.PSet( pattern = cms.string("HLT_DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltDoublePFTau35TrackPt1MediumChargedIsolationDz02Reg"),
              filters2 = cms.untracked.vstring("hltDoublePFTau35TrackPt1MediumChargedIsolationDz02Reg") ),
    cms.PSet( pattern = cms.string("HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg"),
              filters2 = cms.untracked.vstring("hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg") ),
    cms.PSet( pattern = cms.string("HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltDoublePFTau40TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg"),
              filters2 = cms.untracked.vstring("hltDoublePFTau40TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg") ),
    cms.PSet( pattern = cms.string("HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v"),
              nLegs = cms.untracked.uint32(2),
              filters1 = cms.untracked.vstring("hltDoublePFTau40TrackPt1TightChargedIsolationDz02Reg"),
              filters2 = cms.untracked.vstring("hltDoublePFTau40TrackPt1TightChargedIsolationDz02Reg") ),

)

hltPaths_tauTau = {
    'Run2016'   : hltPaths_tauTau_Run2016,
    'Run2017'   : hltPaths_tauTau_Run2017
}

hltPaths_muMu_Run2016 = cms.VPSet(
    cms.PSet( pattern = cms.string("HLT_IsoMu22_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09") )
)

hltPaths_muMu_Run2017 = cms.VPSet(
    cms.PSet( pattern = cms.string("HLT_IsoMu22_v"),
              nLegs = cms.untracked.uint32(1),
              filters1 = cms.untracked.vstring("hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09") )
)

hltPaths_muMu = {
    'Run2016'   : hltPaths_muMu_Run2016,
    'Run2017'   : hltPaths_muMu_Run2017
}

hltPaths = { 'eTau' : hltPaths_eTau, 'muTau' : hltPaths_muTau, 'tauTau' : hltPaths_tauTau, 'muMu' : hltPaths_muMu }
for channelName in hltPaths:
    hltPaths[channelName]['Summer16MC'] = hltPaths[channelName]['Run2016']
    hltPaths[channelName]['Fall17MC'] = hltPaths[channelName]['Run2017']

def IsData(sampleType):
    isData = sampleType in dataSampleTypes
    if not isData and not sampleType in mcSampleTypes:
        print "ERROR: unknown sample type = '{}'".format(sampleType)
        sys.exit(1)
    return isData

def GetHltPaths(channelName, sampleType):
    if not channelName in hltPaths or not sampleType in hltPaths[channelName]:
        print "ERROR: no HLT paths found for sample type '{}' for channel '{}'.".format(sampleType, channel)
        sys.exit(1)
    return hltPaths[channelName][sampleType]
