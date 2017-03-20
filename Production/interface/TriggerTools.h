/*! Tools for trigger selection and matching.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "h-tautau/Cuts/include/H_tautau_2016_baseline.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"
#include "AnalysisTools/Core/include/EnumNameMap.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "h-tautau/Analysis/include/TriggerResults.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1Trigger/interface/Tau.h"

namespace analysis {

enum class CMSSW_Process { SIM, HLT, RECO, PAT, TEST };
ENUM_NAMES(CMSSW_Process) = {
    { CMSSW_Process::SIM, "SIM" },
    { CMSSW_Process::HLT, "HLT" },
    { CMSSW_Process::RECO, "RECO" },
    { CMSSW_Process::PAT, "PAT" },
    { CMSSW_Process::TEST, "TEST" }
};

namespace detail {
template<typename PatObject>
const std::set<trigger::TriggerObjectType>& GetTriggerObjectTypes(const PatObject&);

template<>
inline const std::set<trigger::TriggerObjectType>& GetTriggerObjectTypes<pat::Electron>(const pat::Electron&)
{
    static const std::set<trigger::TriggerObjectType> types = { trigger::TriggerElectron, trigger::TriggerCluster };
    return types;
}

template<>
inline const std::set<trigger::TriggerObjectType>& GetTriggerObjectTypes<pat::Muon>(const pat::Muon&)
{
    static const std::set<trigger::TriggerObjectType> types = { trigger::TriggerMuon };
    return types;
}

template<>
inline const std::set<trigger::TriggerObjectType>& GetTriggerObjectTypes<pat::Tau>(const pat::Tau&)
{
    static const std::set<trigger::TriggerObjectType> types = { trigger::TriggerTau };
    return types;
}

} // namespace detail

class TriggerTools {
public:
    //using L1ParticlePtrSet = std::set<const l1extra::L1JetParticle*>;
    using L1ParticlePtrSet = std::set<const l1t::Tau*>;
    template<typename T> using EDGetTokenT = edm::EDGetTokenT<T>;
    template<typename T> using Handle = edm::Handle<T>;
    using TriggerObjectSet = std::set<const pat::TriggerObjectStandAlone*>;

    TriggerTools(EDGetTokenT<edm::TriggerResults>&& _triggerResultsSIM_token,
                 EDGetTokenT<edm::TriggerResults>&& _triggerResultsHLT_token,
                 EDGetTokenT<edm::TriggerResults>&& _triggerResultsRECO_token,
                 EDGetTokenT<edm::TriggerResults>&& _triggerResultsPAT_token,
                 EDGetTokenT<pat::PackedTriggerPrescales>&& _triggerPrescales_token,
                 EDGetTokenT<pat::TriggerObjectStandAloneCollection>&& _triggerObjects_token,
//                 EDGetTokenT<std::vector<l1extra::L1JetParticle>>&& _l1JetParticles_token);
                 EDGetTokenT<BXVector<l1t::Tau>>&& _l1JetParticles_token);

    TriggerTools(const edm::ParameterSet& iConfig);

    void Initialize(const edm::Event& iEvent);
    
    template<typename HiggsCandidate>
    std::set<const l1t::Tau*> PrintL1Particles(const HiggsCandidate& candidate)
    {
        std::set<const l1t::Tau*> matches;
        const BXVector<l1t::Tau>& l1taus = *l1JetParticles.product();
       
        const auto& candidateMomentums = candidate.GetDaughterMomentums();
        
        std::cout << " l1JetParticle size: " << l1taus.size(0) << std::endl;
        for (unsigned n = 0; n < l1taus.size(0); ++n){
            //std::cout << "n" << n << " - l1JetParticle elem: " << l1taus.at(0,n) << std::endl;
            const l1t::Tau& l1tau = l1taus.at(0,n);
            std::cout << "n" << n << " - l1tau pt: " << l1tau.pt() << std::endl;
            for (const auto& candidateMomentum : candidateMomentums){
                const double deltaR2 = std::pow(0.5, 2);
                if(ROOT::Math::VectorUtil::DeltaR2(l1tau.p4(), candidateMomentum) >= deltaR2) continue;
                matches.insert(&l1tau);
                break;
            }
        }
        return matches;
    }

    void SetTriggerAcceptBits(const analysis::TriggerDescriptors& descriptors, analysis::TriggerResults& results);

    TriggerObjectSet FindMatchingTriggerObjects(const analysis::TriggerDescriptors& descriptors, size_t path_index,
            const std::set<trigger::TriggerObjectType>& objectTypes, const LorentzVector& candidateMomentum,
            size_t leg_id, double deltaR_Limit);

    template<typename Candidate>
    TriggerObjectSet FindMatchingTriggerObjects(const analysis::TriggerDescriptors& descriptors, size_t path_index,
            const Candidate& candidate, size_t leg_id, double deltaR_Limit)
    {
        return FindMatchingTriggerObjects(descriptors, path_index, detail::GetTriggerObjectTypes(*candidate),
                                          candidate.GetMomentum(), leg_id, deltaR_Limit);
    }

    template<typename HiggsCandidate>
    void SetTriggerMatchBits(const analysis::TriggerDescriptors& descriptors, analysis::TriggerResults& results,
                             const HiggsCandidate& candidate, double deltaR_Limit, bool can_flip = false)
    {
        std::cout << "descriptors size: " <<  descriptors.size() << std::endl;
        for(size_t n = 0; n < descriptors.size(); ++n) {
            const size_t n_legs = descriptors.GetNumberOfLegs(n);
            if(n_legs > 2 || n_legs == 0)
                throw exception("Unsupported number of legs = %1%.") % n_legs;
            bool match_found = false;
            const size_t max_flip = can_flip ? 2 : 1;
            for(size_t flip = 0; !match_found && flip < max_flip; ++flip) {
                std::map<size_t, TriggerObjectSet> matches;
                const size_t first = (flip % 2) + 1, second = ((flip + 1) % 2) + 1;
                matches[first] = FindMatchingTriggerObjects(descriptors, n, candidate.GetFirstDaughter(), first,
                                                            deltaR_Limit);
                std::cout << "Found first trigger match " <<  std::endl;
                matches[second] = FindMatchingTriggerObjects(descriptors, n, candidate.GetSecondDaughter(), second,
                                                             deltaR_Limit);

                std::cout << "Found second trigger match " <<  std::endl;
                std::vector<const pat::TriggerObjectStandAlone*> comb_match;
                std::set_union(matches[1].begin(), matches[1].end(), matches[2].begin(), matches[2].end(),
                               std::back_inserter(comb_match));

                match_found = matches[1].size() >= 1 && matches[2].size() >= n_legs - 1 && comb_match.size() >= n_legs;
                std::cout << "Found match " <<  std::endl;
            }
            results.SetMatch(n, match_found);
            
            std::cout << "Set results " <<  std::endl;
        }
    }

    bool TryGetTriggerResult(CMSSW_Process process, const std::string& name, bool& result) const;
    bool GetTriggerResult(CMSSW_Process process, const std::string& name) const;
    bool TryGetAnyTriggerResult(const std::string& name, bool& result) const;
    bool GetAnyTriggerResult(const std::string& name) const;

private:
    std::map<CMSSW_Process, EDGetTokenT<edm::TriggerResults>> triggerResults_tokens;
    EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_token;
    EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_token;
//    EDGetTokenT<std::vector<l1extra::L1JetParticle>> l1JetParticles_token;
    EDGetTokenT<BXVector<l1t::Tau>> l1JetParticles_token;

    const edm::Event* iEvent;
    std::map<CMSSW_Process, Handle<edm::TriggerResults>> triggerResultsMap;
    edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
    edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
//    edm::Handle<std::vector<l1extra::L1JetParticle>> l1JetParticles;
    edm::Handle<BXVector<l1t::Tau>> l1JetParticles;
};

} // namespace analysis
