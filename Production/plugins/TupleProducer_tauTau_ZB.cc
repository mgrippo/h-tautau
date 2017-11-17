/*! Implementation of an event tuple producer for the tau-tau channel.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "../interface/TupleProducer_tauTau_ZB.h"
#include "../interface/GenTruthTools.h"
#include "../interface/TriggerTools.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

void TupleProducer_tauTau_ZB::ProcessEvent(Cutter& cut)
{
    using namespace cuts::H_tautau_2016::TauTau;
    using TriggerObjectSet = std::set<const pat::TriggerObjectStandAlone*>;

    SelectionResults selection(eventId, eventEnergyScale);
    if(applyTriggerMatch) {
        triggerTools.SetTriggerAcceptBits(triggerDescriptors, selection.triggerResults);
        cut(selection.triggerResults.AnyAccpet(), "trigger");
        std::cout << "Trigger match done " <<  std::endl;
    
        triggerTools.SetTriggerMatchBits_data(triggerDescriptors, selection.triggerResults, true);
    }
    
    std::cout << "Trigger Match done " <<  std::endl;


    FillEventTuple(selection);
    std::cout << "Filled eventNtuples done " <<  std::endl;
    
    if(eventEnergyScale == analysis::EventEnergyScale::Central){
        previous_selection = SelectionResultsPtr(new SelectionResults(selection));
        std::cout << "End " <<  std::endl;
    }
    std::cout << "End " <<  std::endl;
}

void TupleProducer_tauTau_ZB::FillEventTuple(const SelectionResults& selection)
{
    using Channel = analysis::Channel;
    using EventPart = ntuple::StorageMode::EventPart;

    BaseTupleProducer::FillEventTuple(selection, previous_selection.get());
    eventTuple().channelId = static_cast<int>(Channel::TauTau_ZB);

    if (applyTriggerMatch){
        std::cout << "Pippo" << std::endl;
        using TriggerObjectSet = std::set<const pat::TriggerObjectStandAlone*>;
        std::cout << "Set TriggerMatchObject - matches size: " << selection.triggerResults.GetTriggerMatchObjects().size() << std::endl;
        std::vector<double> hlt_matches_pt_1;
        std::vector<double> hlt_matches_pt_2;
        for (const auto match : selection.triggerResults.GetTriggerMatchObjects()){
            std::cout << "Match: " << match.first << ", " << match.second.size() << std::endl;
            std::map<size_t, TriggerObjectSet> matches = match.second;
            std::cout << "Match_obj 1: " << matches.at(1).size() << std::endl;
            std::cout << "Match_obj 2: " << matches.at(2).size() << std::endl;
            if (matches.at(1).size() != 0 && matches.at(2).size() != 0){
                
                for (const auto& hlt_tau1 : matches.at(1)){
                    hlt_matches_pt_1.push_back(hlt_tau1->pt());
                }
                
                for (const auto& hlt_tau2 : matches.at(2)){
                    hlt_matches_pt_2.push_back(hlt_tau2->pt());
                }
            }
        }
        
        std::sort(hlt_matches_pt_1.begin(), hlt_matches_pt_1.end(),TupleProducer_tauTau_ZB::PtComparator);
        std::sort(hlt_matches_pt_2.begin(), hlt_matches_pt_2.end(),TupleProducer_tauTau_ZB::PtComparator);
        
        for (unsigned n = 0; n < hlt_matches_pt_1.size(); ++n){
            std::cout << "n: " << n << ", pt: " << hlt_matches_pt_1.at(n) << std::endl;
        }
        for (unsigned n = 0; n < hlt_matches_pt_2.size(); ++n){
            std::cout << "n: " << n << ", pt: " << hlt_matches_pt_2.at(n) << std::endl;
        }
        
        LorentzVector hlt_tau1_momentum, hlt_tau2_momentum;
        for (const auto match : selection.triggerResults.GetTriggerMatchObjects()){
            std::map<size_t, TriggerObjectSet> matches = match.second;
            if (matches.at(1).size() != 0 && matches.at(2).size() != 0){
                
                for (const auto& hlt_tau1 : matches.at(1)){
                    if (hlt_matches_pt_1.at(0) == hlt_tau1->pt()){
                        std::cout << "hlt_matches_pt_1: " << hlt_matches_pt_1.at(0) << ", hlt_tau1_pt(): " << hlt_tau1->pt() << std::endl;
                        std::cout << "Found 1" << std::endl;
                        eventTuple().hlt_match_p4_1.push_back(ntuple::LorentzVectorE(hlt_tau1->p4()));
                        hlt_tau1_momentum = hlt_tau1->p4();
                    }
                }
                
                for (const auto& hlt_tau2 : matches.at(2)){
                    if (hlt_matches_pt_2.at(1) == hlt_tau2->pt()){
                        std::cout << "hlt_matches_pt_2: " << hlt_matches_pt_2.at(1) << ", hlt_tau2_pt(): " << hlt_tau2->pt() << std::endl;
                        std::cout << "Found 2" << std::endl;
                        eventTuple().hlt_match_p4_2.push_back(ntuple::LorentzVectorE(hlt_tau2->p4()));
                        hlt_tau2_momentum = hlt_tau2->p4();
                    }
                }
            }
        }
        
        const auto L1tauMatches = triggerTools.MatchL1Particles_data(hlt_tau1_momentum, hlt_tau2_momentum);
        std::cout << "L1tauMatches: " << L1tauMatches.size() << std::endl;
        const std::set<const l1t::Tau*> l1tauMatches_1 = L1tauMatches.at(0);
        const std::set<const l1t::Tau*> l1tauMatches_2 = L1tauMatches.at(1);
        std::cout << "l1tauMatches_1: " << l1tauMatches_1.size() << ", l1tauMatches_2: " << l1tauMatches_2.size() << std::endl;
        if (l1tauMatches_1.size() != 0 && l1tauMatches_2.size() != 0){
            for (const auto& tau1 : l1tauMatches_1){
                std::cout << "Found 1 l1" << std::endl;
                eventTuple().l1_match_p4_1.push_back(ntuple::LorentzVectorE(tau1->p4()));
                eventTuple().l1_hwIso_1.push_back(tau1->hwIso());
            }
            
            for (const auto& tau2 : l1tauMatches_2){
                std::cout << "Found 2 l1" << std::endl;
                eventTuple().l1_match_p4_2.push_back(ntuple::LorentzVectorE(tau2->p4()));
                eventTuple().l1_hwIso_2.push_back(tau2->hwIso());
            }
        }


    }

    eventTuple.Fill();
}

bool TupleProducer_tauTau_ZB::PtComparator(double pt_1, double pt_2)
{
    return pt_1 > pt_2;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TupleProducer_tauTau_ZB);
