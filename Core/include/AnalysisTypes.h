/*! Common simple types for analysis purposes.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/PhysicalValue.h"
#include "AnalysisTools/Core/include/EnumNameMap.h"

namespace analysis {

enum class Channel { ETau = 0, MuTau = 1, TauTau = 2, MuMu = 3 };
ENUM_NAMES(Channel) = {
    { Channel::ETau, "eTau" }, { Channel::MuTau, "muTau" }, { Channel::TauTau, "tauTau" }, { Channel::MuMu, "muMu" }
};
const EnumNameMap<Channel> __Channel_names_latex("ChannelLatex", {
    { Channel::ETau, "e#tau_{h}" }, { Channel::MuTau, "#mu#tau_{h}" }, { Channel::TauTau, "#tau_{h}#tau_{h}" },
    { Channel::MuMu, "#mu#mu" }
});

enum class LegType { e = 0, mu = 1, tau = 2, jet = 3 };
ENUM_NAMES(LegType) = {
    { LegType::e, "e" }, { LegType::mu, "mu" }, { LegType::tau, "tau" }, { LegType::jet, "jet" }
};

using ChannelLegTypes = std::pair<LegType, LegType>;
inline const ChannelLegTypes GetChannelLegTypes(Channel channel)
{
    static const std::map<Channel, ChannelLegTypes> leg_types {
        { Channel::ETau, { LegType::e, LegType::tau } }, { Channel::MuTau, { LegType::mu, LegType::tau } },
        { Channel::TauTau, { LegType::tau, LegType::tau } }, { Channel::MuMu, { LegType::mu, LegType::mu } },
    };
    auto iter = leg_types.find(channel);
    if(iter == leg_types.end())
        throw exception("Leg types are not defined for channel '%1%'") % channel;
    return iter->second;
}

enum class UncertaintyScale { Central = 0, Up = 1, Down = -1 };
ENUM_NAMES(UncertaintyScale) = {
    { UncertaintyScale::Central, "Central" }, { UncertaintyScale::Up, "Up" }, { UncertaintyScale::Down, "Down" }
};

enum class UncertaintySource { None = 0, TauES = 1, JetFull_Total = 2, TopPt = 3, JetFull_AbsoluteStat = 4, JetFull_AbsoluteScale = 5,
     JetFull_AbsoluteMPFBias = 6, JetFull_AbsoluteFlavMap = 7, JetFull_Fragmentation = 8, JetFull_SinglePionECAL = 9, JetFull_SinglePionHCAL = 10,
     JetFull_FlavorQCD = 11, JetFull_FlavorZJet = 12, JetFull_FlavorPhotonJet = 13, JetFull_FlavorPureGluon = 14, JetFull_FlavorPureQuark = 15, JetFull_FlavorPureCharm = 16,
     JetFull_FlavorPureBottom = 17, JetFull_TimePtEta = 18, JetFull_RelativeJEREC1 = 19, JetFull_RelativeJEREC2 = 20, JetFull_RelativeJERHF = 21,
     JetFull_RelativePtBB = 22, JetFull_RelativePtEC1 = 23, JetFull_RelativePtEC2 = 24, JetFull_RelativePtHF = 25, JetFull_RelativeBal = 26,
     JetFull_RelativeFSR = 27, JetFull_PileUpDataMC = 28, JetFull_PileUpPtRef = 29, JetFull_PileUpPtBB = 30, JetFull_PileUpPtEC1 = 31, JetFull_PileUpPtEC2 = 32,
     JetFull_PileUpPtHF = 33, JetFull_SubTotalPileUp = 34, JetFull_SubTotalRelative = 35, JetFull_SubTotalPt = 36, JetFull_SubTotalScale = 37,
     JetFull_SubTotalAbsolute = 38, JetFull_SubTotalMC = 39, JetFull_TotalNoFlavor = 40, JetFull_TotalNoTime = 41, JetFull_TotalNoFlavorNoTime = 42,
     EleFakingTauES = 43, JetReduced_Absolute = 44, JetReduced_Absolute_year = 45,
     JetReduced_BBEC1 = 46, JetReduced_BBEC1_year = 47, JetReduced_EC2 = 48, JetReduced_EC2_year = 49,
     JetReduced_FlavorQCD = 50, JetReduced_HF = 51, JetReduced_HF_year = 52, JetReduced_RelativeBal = 53,
     JetReduced_RelativeSample_year = 54, JetReduced_Total = 55 };
ENUM_NAMES(UncertaintySource) = {
    { UncertaintySource::None, "None" }, { UncertaintySource::TauES, "TauES" },
    { UncertaintySource::JetFull_Total, "JetFull_Total" }, { UncertaintySource::TopPt, "TopPt" },
    { UncertaintySource::JetFull_AbsoluteStat, "JetFull_AbsoluteStat" },
    { UncertaintySource::JetFull_AbsoluteScale, "JetFull_AbsoluteScale" },
    { UncertaintySource::JetFull_AbsoluteMPFBias, "JetFull_AbsoluteMPFBias" },
    { UncertaintySource::JetFull_AbsoluteFlavMap, "JetFull_AbsoluteFlavMap" },
    { UncertaintySource::JetFull_Fragmentation, "JetFull_Fragmentation" }, { UncertaintySource::JetFull_SinglePionECAL, "JetFull_SinglePionECAL" },
    { UncertaintySource::JetFull_SinglePionHCAL, "JetFull_SinglePionHCAL" },
    { UncertaintySource::JetFull_FlavorQCD, "JetFull_FlavorQCD" },
    { UncertaintySource::JetFull_FlavorZJet, "JetFull_FlavorZJet" },
    { UncertaintySource::JetFull_FlavorPhotonJet, "JetFull_FlavorPhotonJet" },
    { UncertaintySource::JetFull_FlavorPureGluon, "JetFull_FlavorPureGluon" },
    { UncertaintySource::JetFull_FlavorPureQuark, "JetFull_FlavorPureQuark" },
    { UncertaintySource::JetFull_FlavorPureCharm, "JetFull_FlavorPureCharm" },
    { UncertaintySource::JetFull_FlavorPureBottom, "JetFull_FlavorPureBottom" },
    { UncertaintySource::JetFull_TimePtEta, "JetFull_TimePtEta" },
    { UncertaintySource::JetFull_RelativeJEREC1, "JetFull_RelativeJEREC1" },
    { UncertaintySource::JetFull_RelativeJEREC2, "JetFull_RelativeJEREC2" },
    { UncertaintySource::JetFull_RelativeJERHF, "JetFull_RelativeJERHF" },
    { UncertaintySource::JetFull_RelativePtBB, "JetFull_RelativePtBB" },
    { UncertaintySource::JetFull_RelativePtEC1, "JetFull_RelativePtEC1" },
    { UncertaintySource::JetFull_RelativePtEC2, "JetFull_RelativePtEC2" },
    { UncertaintySource::JetFull_RelativePtHF, "JetFull_RelativePtHF" },
    { UncertaintySource::JetFull_RelativeBal, "JetFull_RelativeBal" },
    { UncertaintySource::JetFull_RelativeFSR, "JetFull_RelativeFSR" },
    { UncertaintySource::JetFull_PileUpDataMC, "JetFull_PileUpDataMC" },
    { UncertaintySource::JetFull_PileUpPtRef, "JetFull_PileUpPtRef" },
    { UncertaintySource::JetFull_PileUpPtBB, "JetFull_PileUpPtBB" },
    { UncertaintySource::JetFull_PileUpPtEC1, "JetFull_PileUpPtEC1" },
    { UncertaintySource::JetFull_PileUpPtEC2, "JetFull_PileUpPtEC2" },
    { UncertaintySource::JetFull_PileUpPtHF, "JetFull_PileUpPtHF" },
    { UncertaintySource::JetFull_SubTotalPileUp, "JetFull_SubTotalPileUp" },
    { UncertaintySource::JetFull_SubTotalRelative, "JetFull_SubTotalRelative" },
    { UncertaintySource::JetFull_SubTotalPt, "JetFull_SubTotalPt" },
    { UncertaintySource::JetFull_SubTotalScale, "JetFull_SubTotalScale" },
    { UncertaintySource::JetFull_SubTotalAbsolute, "JetFull_SubTotalAbsolute" },
    { UncertaintySource::JetFull_SubTotalMC, "JetFull_SubTotalMC" },
    { UncertaintySource::JetFull_TotalNoFlavor, "JetFull_TotalNoFlavor"},
    { UncertaintySource::JetFull_TotalNoTime, "JetFull_TotalNoTime"},
    { UncertaintySource::JetFull_TotalNoFlavorNoTime, "JetFull_TotalNoFlavorNoTime"},
    { UncertaintySource::EleFakingTauES, "EleFakingTauES"},
    { UncertaintySource::JetReduced_Absolute, "JetReduced_Absolute"},
    { UncertaintySource::JetReduced_Absolute_year, "JetReduced_Absolute_year"},
    { UncertaintySource::JetReduced_BBEC1, "JetReduced_BBEC1"},
    { UncertaintySource::JetReduced_BBEC1_year, "JetReduced_BBEC1_year"},
    { UncertaintySource::JetReduced_EC2, "JetReduced_EC2"},
    { UncertaintySource::JetReduced_EC2_year, "JetReduced_EC2_year"},
    { UncertaintySource::JetReduced_FlavorQCD, "JetReduced_FlavorQCD"},
    { UncertaintySource::JetReduced_HF, "JetReduced_HF"},
    { UncertaintySource::JetReduced_HF_year, "JetReduced_HF_year"},
    { UncertaintySource::JetReduced_RelativeBal, "JetReduced_RelativeBal"},
    { UncertaintySource::JetReduced_RelativeSample_year, "JetReduced_RelativeSample_year"},
    { UncertaintySource::JetReduced_Total, "JetReduced_Total"}

};

const std::set<UncertaintyScale>& GetAllUncertaintyScales();
const std::set<UncertaintyScale>& GetActiveUncertaintyScales(UncertaintySource unc_source);

enum class DiscriminatorWP { VVVLoose = 0, VVLoose = 1, VLoose = 2, Loose = 3, Medium = 4, Tight = 5,
                             VTight = 6, VVTight = 7, VVVTight = 8 };
ENUM_NAMES(DiscriminatorWP) = {
    { DiscriminatorWP::VVVLoose, "VVVLoose" }, { DiscriminatorWP::VVLoose, "VVLoose" },
    { DiscriminatorWP::VLoose, "VLoose" }, { DiscriminatorWP::Loose, "Loose" }, { DiscriminatorWP::Medium, "Medium" },
    { DiscriminatorWP::Tight, "Tight" }, { DiscriminatorWP::VTight, "VTight" }, { DiscriminatorWP::VVTight, "VVTight" },
    { DiscriminatorWP::VVVTight, "VVVTight" }
};
const EnumNameMap<DiscriminatorWP> __DiscriminatorWP_short_names("ShortWPNames", {
    { DiscriminatorWP::VVVLoose, "VVVL" }, { DiscriminatorWP::VVLoose, "VVL" }, { DiscriminatorWP::VLoose, "VL" },
    { DiscriminatorWP::Loose, "L" }, { DiscriminatorWP::Medium, "M" }, { DiscriminatorWP::Tight, "T" },
    { DiscriminatorWP::VTight, "VT" }, { DiscriminatorWP::VVTight, "VVT" }, { DiscriminatorWP::VVVTight, "VVVT" }
});

enum class MetType { PF, MVA, PUPPI };
ENUM_NAMES(MetType) = {
    { MetType::PF, "PF" }, { MetType::MVA, "MVA" }, { MetType::PUPPI, "PUPPI" }
};

enum class Period { Run2015, Run2016, Run2017, Run2018 };
ENUM_NAMES(Period) = {
    { Period::Run2015, "Run2015" },
    { Period::Run2016, "Run2016" },
    { Period::Run2017, "Run2017" },
    { Period::Run2018, "Run2018" },
};

enum class GenLeptonMatch { Electron = 1, Muon = 2, TauElectron = 3,  TauMuon = 4, Tau = 5, NoMatch = 6 };
ENUM_NAMES(GenLeptonMatch) = {
    { GenLeptonMatch::Electron, "gen_electron" },
    { GenLeptonMatch::Muon, "gen_muon" },
    { GenLeptonMatch::TauElectron, "gen_electron_from_tau" },
    { GenLeptonMatch::TauMuon, "gen_muon_from_tau" },
    { GenLeptonMatch::Tau, "gen_tau" },
    { GenLeptonMatch::NoMatch, "no_gen_match" }
};

enum class GenEventType { Other = 0, TTbar_Hadronic = 1, TTbar_SemiLeptonic = 2, TTbar_Leptonic = 3 };
ENUM_NAMES(GenEventType) = {
    { GenEventType::Other, "other" },
    { GenEventType::TTbar_Hadronic, "TTbar_Hadronic" },
    { GenEventType::TTbar_SemiLeptonic, "TTbar_SemiLeptonic" },
    { GenEventType::TTbar_Leptonic, "TTbar_Leptonic" },
};

} // namespace analysis
