/*! Implementation of an event tuple producer for the tau-tau channel.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "../interface/TupleProducer_tauTau.h"
#include "../interface/GenTruthTools.h"

void TupleProducer_tauTau::ProcessEvent(Cutter& cut)
{
    using namespace cuts::H_tautau_2016::TauTau;

    SelectionResults selection(eventId, eventEnergyScale);
    cut(primaryVertex.isNonnull(), "vertex");

    std::cout << "Vertex done " <<  std::endl;
    if(applyTriggerMatch) {
        triggerTools.SetTriggerAcceptBits(triggerDescriptors, selection.triggerResults);
        cut(selection.triggerResults.AnyAccpet(), "trigger");
        std::cout << "Trigger match done " <<  std::endl;
    }

    const auto selectedTaus = CollectSignalTaus();
    cut(selectedTaus.size(), "taus");
    std::cout << "Taus done " <<  std::endl;

    const double DeltaR_betweenSignalObjects = productionMode == ProductionMode::hh
            ? cuts::hh_bbtautau_2016::DeltaR_betweenSignalObjects
            : cuts::H_tautau_2016::DeltaR_betweenSignalObjects;
    auto higgses = FindCompatibleObjects(selectedTaus, selectedTaus, DeltaR_betweenSignalObjects, "H_tau_tau");
    cut(higgses.size(), "tau_tau_pair");
    std::cout << "Higgs pair done " <<  std::endl;

    std::sort(higgses.begin(), higgses.end(), &HiggsComparitor<HiggsCandidate>);
    auto selected_higgs = higgses.front();
    if (selected_higgs.GetFirstDaughter().GetMomentum().Pt() < selected_higgs.GetSecondDaughter().GetMomentum().Pt())
        selected_higgs = HiggsCandidate(selected_higgs.GetSecondDaughter(), selected_higgs.GetFirstDaughter());
    
    std::cout << "Higgs candidate selected " <<  std::endl;

    if(applyTriggerMatch){
        triggerTools.SetTriggerMatchBits(triggerDescriptors, selection.triggerResults, selected_higgs,
                                         cuts::H_tautau_2016::DeltaR_triggerMatch, true);
        
        const auto L1tauMatches = triggerTools.MatchL1Particles(selected_higgs,selection.triggerResults);
        
    }
    
    std::cout << "Trigger Match done " <<  std::endl;

    selection.SetHiggsCandidate(selected_higgs);
    std::cout << "Set Higgs candidate done " <<  std::endl;

    //Third-Lepton Veto
    //const auto electronVetoCollection = CollectVetoElectrons();
    const auto muonVetoCollection = CollectVetoMuons();
    //selection.electronVeto = electronVetoCollection.size();
    selection.electronVeto = false;
    selection.muonVeto = muonVetoCollection.size();
    std::cout << "Vetoes done " <<  std::endl;

    ApplyBaseSelection(selection, selection.higgs->GetDaughterMomentums());
    std::cout << "Applied baseline selection " <<  std::endl;
    if(runSVfit)
        selection.svfitResult = svfitProducer->Fit(*selection.higgs, *met);
    FillEventTuple(selection);
    std::cout << "Filled eventNtuples done " <<  std::endl;
    
    if(eventEnergyScale == analysis::EventEnergyScale::Central){
        previous_selection = SelectionResultsPtr(new SelectionResults(selection));
        std::cout << "End " <<  std::endl;
    }
    std::cout << "End " <<  std::endl;
}

std::vector<BaseTupleProducer::TauCandidate> TupleProducer_tauTau::CollectSignalTaus()
{
    using namespace std::placeholders;
    const auto base_selector = std::bind(&TupleProducer_tauTau::SelectSignalTau, this, _1, _2);
    return CollectObjects("SignalTaus", base_selector, taus);
}

void TupleProducer_tauTau::SelectSignalTau(const TauCandidate& tau, Cutter& cut) const
{
    using namespace cuts::H_tautau_2016::TauTau::tauID;

    cut(true, "gt0_tau_cand");
    const LorentzVector& p4 = tau.GetMomentum();
//    cut(p4.Pt() > pt, "pt", p4.Pt());
    cut(std::abs(p4.Eta()) < eta, "eta", p4.Eta());
    const auto dmFinding = tau->tauID("decayModeFinding");
    cut(dmFinding > decayModeFinding, "oldDecayMode", dmFinding);
    const auto packedLeadTauCand = dynamic_cast<const pat::PackedCandidate*>(tau->leadChargedHadrCand().get());
    cut(std::abs(packedLeadTauCand->dz()) < dz, "dz", packedLeadTauCand->dz());
    cut(std::abs(tau->charge()) == absCharge, "charge", tau->charge());
    if(productionMode == ProductionMode::hh) {
        cut(tau->tauID("againstElectronVLooseMVA6") > againstElectronVLooseMVA6, "againstElectron");
        cut(tau->tauID("againstMuonLoose3") > againstMuonLoose3, "againstMuon");
    }
}

void TupleProducer_tauTau::FillEventTuple(const SelectionResults& selection)
{
    using Channel = analysis::Channel;
    using EventPart = ntuple::StorageMode::EventPart;

    BaseTupleProducer::FillEventTuple(selection, previous_selection.get());
    eventTuple().channelId = static_cast<int>(Channel::TauTau);

    ntuple::StorageMode storageMode(eventTuple().storageMode);
    const bool store_tauIds_1 = !previous_selection || !selection.HaveSameFirstLegOrigin(*previous_selection);
    const bool store_tauIds_2 = !previous_selection || !selection.HaveSameSecondLegOrigin(*previous_selection);
    storageMode.SetPresence(EventPart::FirstTauIds, store_tauIds_1);
    storageMode.SetPresence(EventPart::SecondTauIds, store_tauIds_2);
    eventTuple().storageMode = storageMode.Mode();

    FillTauLeg(1, selection.higgs->GetFirstDaughter(), store_tauIds_1);
    FillTauLeg(2, selection.higgs->GetSecondDaughter(), store_tauIds_2);
    
    using TriggerObjectSet = std::set<const pat::TriggerObjectStandAlone*>;
    if (applyTriggerMatch){
        
        //for (unsigned n = 0; n < selection.triggerResults.GetL1TriggerMatchObjects().size(); ++n){}
            //std::cout << " l1tau matches size: " << l1tauMatch.size() << std::endl;
        const std::set<const l1t::Tau*> l1tauMatches_1 = selection.triggerResults.GetL1TriggerMatchObjects().at(0);
        const std::set<const l1t::Tau*> l1tauMatches_2 = selection.triggerResults.GetL1TriggerMatchObjects().at(1);
        if (l1tauMatches_1.size() != 0 && l1tauMatches_2.size() != 0){
            for (const auto& tau1 : l1tauMatches_1){
                eventTuple().l1_match_p4_1.push_back(ntuple::LorentzVectorE(tau1->p4()));
                eventTuple().l1_hwIso_1.push_back(tau1->hwIso());
            }
        
            for (const auto& tau2 : l1tauMatches_2){
                eventTuple().l1_match_p4_2.push_back(ntuple::LorentzVectorE(tau2->p4()));
                eventTuple().l1_hwIso_2.push_back(tau2->hwIso());
            }
        }

        std::cout << "Set TriggerMatchObject - matches size: " << selection.triggerResults.GetTriggerMatchObjects().size() << std::endl;
        for (const auto match : selection.triggerResults.GetTriggerMatchObjects()){
            std::cout << "Match: " << match.first << ", " << match.second.size() << std::endl;
            std::map<size_t, TriggerObjectSet> matches = match.second;
            std::cout << "Match_obj 1: " << matches.at(1).size() << std::endl;
            std::cout << "Match_obj 2: " << matches.at(2).size() << std::endl;
            if (matches.at(1).size() != 0 && matches.at(2).size() != 0){
                for (const auto& hlt_tau1 : matches.at(1)){
                    eventTuple().hlt_match_p4_1.push_back(ntuple::LorentzVectorE(hlt_tau1->p4()));
                }
                
                for (const auto& hlt_tau2 : matches.at(2)){
                    eventTuple().hlt_match_p4_2.push_back(ntuple::LorentzVectorE(hlt_tau2->p4()));
                }
            }
        }
    }

    eventTuple.Fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TupleProducer_tauTau);
