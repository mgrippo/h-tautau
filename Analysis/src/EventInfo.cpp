/*! Definiton of analysis::FlatEventInfo class.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/Analysis/include/EventInfo.h"

namespace analysis {

SummaryInfo::SummaryInfo(const ProdSummary& _summary, const Channel& _channel, const std::string& _uncertainties_source,
                         const std::string& _trigger_cfg) :
                         summary(_summary)
{
    if(!_uncertainties_source.empty())
        jecUncertainties = std::make_shared<jec::JECUncertaintiesWrapper>(_uncertainties_source);
    if(!_trigger_cfg.empty()){
        triggerDescriptors[_channel] = TriggerDescriptorCollection::Load(_trigger_cfg,_channel);
    }

}

std::shared_ptr<const TriggerDescriptorCollection> SummaryInfo::GetTriggerDescriptors(Channel channel) const
{
    if(!triggerDescriptors.count(channel))
        throw exception("Information for channel %1% not found.") % channel;
    return triggerDescriptors.at(channel);
}

const SummaryInfo::ProdSummary& SummaryInfo::operator*() const { return summary; }
const SummaryInfo::ProdSummary* SummaryInfo::operator->() const { return &summary; }

const jec::JECUncertaintiesWrapper& SummaryInfo::GetJecUncertainties() const
{
    if(!jecUncertainties)
        throw exception("Jec Uncertainties not stored.");
    return *jecUncertainties;
}

std::array<size_t,2> EventInfoBase::GetSelectedBjetIndices() const
{
    std::array<size_t,2> bjet_indexes;
    bjet_indexes[0] = selected_signal_jets.selectedBjetPair.first;
    bjet_indexes[1] = selected_signal_jets.selectedBjetPair.second;
    return bjet_indexes;
}

std::set<size_t> EventInfoBase::GetSelectedBjetIndicesSet() const
{
    std::set<size_t> bjet_indexes;
    bjet_indexes.insert(selected_signal_jets.selectedBjetPair.first);
    bjet_indexes.insert(selected_signal_jets.selectedBjetPair.second);
    return bjet_indexes;
}

const analysis::LepCandidate& EventInfoBase::GetFirstLeg()
{
    Lock lock(*mutex);
    if(!leg1) {
        tuple_leg1 = std::make_shared<ntuple::TupleLepton>(*event, GetLegIndex(1));
        leg1 = std::make_shared<analysis::LepCandidate>(*tuple_leg1, tuple_leg1->iso());
    }
    return *leg1;
}

const analysis::LepCandidate& EventInfoBase::GetSecondLeg()
{
    Lock lock(*mutex);
    if(!leg2) {
        tuple_leg2 = std::make_shared<ntuple::TupleLepton>(*event, GetLegIndex(2));
        leg2 = std::make_shared<analysis::LepCandidate>(*tuple_leg2, tuple_leg2->iso());
    }
    return *leg2;
}

std::shared_ptr<EventInfoBase> EventInfoBase::ApplyShift(UncertaintySource uncertainty_source,
    UncertaintyScale scale)
{
    std::shared_ptr<EventInfoBase> shifted_event_info = std::make_shared<EventInfoBase>(*this);
    const SummaryInfo& summaryInfo = shifted_event_info->GetSummaryInfo();
    const jec::JECUncertaintiesWrapper& jecUncertainties = summaryInfo.GetJecUncertainties();
    const JetCollection& jets = shifted_event_info->GetJets();
    const auto& other_jets_p4 = event->other_jets_p4;
    auto shifted_met_p4(shifted_event_info->GetMET().GetMomentum());
    const JetCollection& corrected_jets = jecUncertainties.ApplyShift(jets,uncertainty_source,scale,&other_jets_p4,&shifted_met_p4);
    shifted_event_info->SetJets(corrected_jets);
    shifted_event_info->SetMetMomentum(shifted_met_p4);
    return shifted_event_info;
}

EventInfoBase::EventInfoBase(const Event& _event, const SummaryInfo* _summaryInfo,
                             size_t _selected_htt_index, const SignalObjectSelector::SelectedSignalJets& _selected_signal_jets,
                             Period _period, JetOrdering _jet_ordering) :
event(&_event), summaryInfo(_summaryInfo), selected_htt_index(_selected_htt_index), eventIdentifier(_event.run, _event.lumi, _event.evt),
 selected_signal_jets(_selected_signal_jets), period(_period), jet_ordering(_jet_ordering)
{
    mutex = std::make_shared<Mutex>();
    triggerResults.SetAcceptBits(event->trigger_accepts);
    triggerResults.SetMatchBits(event->trigger_matches.at(selected_htt_index));
    if(summaryInfo)
        triggerResults.SetDescriptors(summaryInfo->GetTriggerDescriptors(EventInfoBase::GetChannel()));
}


const EventInfoBase::Event& EventInfoBase::operator*() const { return *event; }
const EventInfoBase::Event* EventInfoBase::operator->() const { return event; }

const EventIdentifier& EventInfoBase::GetEventId() const { return eventIdentifier; }
EventEnergyScale EventInfoBase::GetEnergyScale() const
{
    return static_cast<EventEnergyScale>(event->eventEnergyScale);
}

const TriggerResults& EventInfoBase::GetTriggerResults() const { return triggerResults; }
const SummaryInfo& EventInfoBase::GetSummaryInfo() const
{
    if(!summaryInfo)
        throw exception("SummaryInfo was not provided for this event.");
    return *summaryInfo;
}

const kin_fit::FitProducer& EventInfoBase::GetKinFitProducer()
{
    static kin_fit::FitProducer kinfitProducer;
    return kinfitProducer;
}

const sv_fit_ana::FitProducer& EventInfoBase::GetSVFitProducer()
{
    static sv_fit_ana::FitProducer svfitProducer;
    return svfitProducer;
}

// const AnalysisObject& EventInfoBase::GetLeg(size_t /*leg_id*/) { throw exception("Method not supported."); }
// LorentzVector EventInfoBase::GetHiggsTTMomentum(bool /*useSVfit*/) { throw exception("Method not supported."); }

size_t EventInfoBase::GetNJets() const { return event->jets_p4.size(); }
size_t EventInfoBase::GetNFatJets() const { return event->fatJets_p4.size(); }

const SignalObjectSelector::SelectedSignalJets& EventInfoBase::GetSelectedSignalJets() const { return selected_signal_jets; }
Period EventInfoBase::GetPeriod() const { return period; }
JetOrdering EventInfoBase::GetJetOrdering() const {return jet_ordering; }

const EventInfoBase::JetCollection& EventInfoBase::GetJets()
{
    Lock lock(*mutex);
    if(!jets) {
        jets = std::shared_ptr<JetCollection>(new JetCollection());
        tuple_jets = std::make_shared<std::list<ntuple::TupleJet>>();
        for(size_t n = 0; n < GetNJets(); ++n) {
            tuple_jets->push_back(ntuple::TupleJet(*event, n));
            jets->push_back(JetCandidate(tuple_jets->back()));
        }
    }
    return *jets;
}

void EventInfoBase::SetJets(const JetCollection& new_jets)
{
    Lock lock(*mutex);
    jets = std::make_shared<JetCollection>(new_jets);
}

EventInfoBase::JetCollection EventInfoBase::SelectJets(double pt_cut, double eta_cut, bool applyPu,
                                                       bool passBtag, JetOrdering jet_ordering,
                                                       const std::set<size_t>& jet_to_exclude_indexes,
                                                       double low_eta_cut)
{
    Lock lock(*mutex);
    BTagger bTagger(period,jet_ordering);
    const JetCollection& all_jets = GetJets();
    JetCollection selected_jets;
    std::vector<analysis::jet_ordering::JetInfo<LorentzVector>> jet_info_vector;
    for(size_t n = 0; n < all_jets.size(); ++n) {
        const JetCandidate& jet = all_jets.at(n);
        if(ROOT::Math::VectorUtil::DeltaR(GetLeg(1).GetMomentum(), jet.GetMomentum()) <= cuts::H_tautau_2016::DeltaR_betweenSignalObjects) continue;
        if(ROOT::Math::VectorUtil::DeltaR(GetLeg(2).GetMomentum(), jet.GetMomentum()) <= cuts::H_tautau_2016::DeltaR_betweenSignalObjects) continue;
        if(!SignalObjectSelector::PassEcalNoiceVetoJets(jet.GetMomentum(), period, event->jets_pu_id.at(n) )) continue;
        if(jet_to_exclude_indexes.count(n)) continue;
        if(applyPu && (event->jets_pu_id.at(n) & (1 << 2)) == 0) continue;
        if(std::abs(jet.GetMomentum().eta()) < low_eta_cut) continue;
        if(passBtag && !bTagger.Pass(*event,n,DiscriminatorWP::Medium)) continue;

        jet_info_vector.emplace_back(jet.GetMomentum(),n,bTagger.BTag(*event,n));
    }
    auto jets_ordered = jet_ordering::OrderJets(jet_info_vector,true,pt_cut,eta_cut);
    for(size_t h = 0; h < jets_ordered.size(); ++h){
        const JetCandidate& jet = all_jets.at(jets_ordered.at(h).index);
        selected_jets.push_back(jet);
    }
    return selected_jets;
}

double EventInfoBase::GetHT(bool includeHbbJets, bool apply_eta_cut)
{
    static constexpr double other_jets_min_pt = 20;
    static constexpr double other_jets_max_eta = 4.7;
    static const std::set<size_t> empty_set = {};

    const double eta_cut = apply_eta_cut ? other_jets_max_eta : 5;
    const std::set<size_t>& jets_to_exclude = includeHbbJets ? empty_set : GetSelectedBjetIndicesSet();

    double ht = 0;
    const auto& jets = SelectJets(other_jets_min_pt,eta_cut,true,false,JetOrdering::DeepCSV,jets_to_exclude);
    for(size_t n = 0; n < jets.size(); ++n) {
        const auto& jet = jets.at(n);
        ht += jet.GetMomentum().pt();
    }
    return ht;
}

const EventInfoBase::FatJetCollection& EventInfoBase::GetFatJets()
{
    Lock lock(*mutex);
    if(!fatJets) {
        fatJets = std::shared_ptr<FatJetCollection>(new FatJetCollection());
        tuple_fatJets = std::make_shared<std::list<ntuple::TupleFatJet>>();
        for(size_t n = 0; n < GetNFatJets(); ++n) {
            tuple_fatJets->push_back(ntuple::TupleFatJet(*event, n));
            fatJets->push_back(FatJetCandidate(tuple_fatJets->back()));
        }
    }
    return *fatJets;
}

bool EventInfoBase::HasBjetPair() const { return selected_signal_jets.HasBjetPair(GetNJets()); }
bool EventInfoBase::HasVBFjetPair() const { return selected_signal_jets.HasVBFPair(GetNJets()); }

const JetCandidate& EventInfoBase::GetVBFJet(const size_t index)
{
    if(!HasVBFjetPair() || (index != 1 && index != 2))
        throw exception("VBF jet not found.");
    if(index == 1)
        return GetJets().at(selected_signal_jets.selectedVBFjetPair.first);
    return GetJets().at(selected_signal_jets.selectedVBFjetPair.second);
}

const JetCandidate& EventInfoBase::GetBJet(const size_t index)
{
    if(!HasBjetPair() || (index != 1 && index != 2) )
        throw exception("B jet not found.");
    if(index == 1)
        return GetJets().at(selected_signal_jets.selectedBjetPair.first);
    return GetJets().at(selected_signal_jets.selectedBjetPair.second);
}

const EventInfoBase::HiggsBBCandidate& EventInfoBase::GetHiggsBB()
{
    Lock lock(*mutex);
    if(!HasBjetPair())
        throw exception("Can't create H->bb candidate.");
    if(!higgs_bb) {
        const auto& jets = GetJets();
        higgs_bb = std::make_shared<HiggsBBCandidate>(jets.at(selected_signal_jets.selectedBjetPair.first),
                                                      jets.at(selected_signal_jets.selectedBjetPair.second));
    }
    return *higgs_bb;
}

const MET& EventInfoBase::GetMET()
{
    Lock lock(*mutex);
    if(!met) {
        tuple_met = std::shared_ptr<ntuple::TupleMet>(new ntuple::TupleMet(*event, MetType::PF));
        met = std::shared_ptr<MET>(new MET(*tuple_met, tuple_met->cov()));
    }
    return *met;
}

size_t EventInfoBase::GetLegIndex(size_t leg_id)
{
    if(leg_id == 1) return event->first_daughter_indexes.at(selected_htt_index);
    if(leg_id == 2) return event->second_daughter_indexes.at(selected_htt_index);
    throw exception("Invalid leg id = %1%.") % leg_id;
}

const kin_fit::FitResults& EventInfoBase::GetKinFitResults()
{
    Lock lock(*mutex);
    if(!HasBjetPair())
        throw exception("Can't retrieve KinFit results.");
    if(!kinfit_results) {
        const size_t pairId = ntuple::CombinationPairToIndex(selected_signal_jets.selectedBjetPair, GetNJets());
        const auto iter = std::find(event->kinFit_jetPairId.begin(), event->kinFit_jetPairId.end(), pairId);
        kinfit_results = std::make_shared<kin_fit::FitResults>();
        if(iter == event->kinFit_jetPairId.end()){
            double energy_resolution_1 = GetBJet(1)->resolution() * GetBJet(1).GetMomentum().E();
            double energy_resolution_2 = GetBJet(2)->resolution() * GetBJet(2).GetMomentum().E();
            const auto& kinfitProducer = GetKinFitProducer();
            const auto& result = kinfitProducer.Fit(GetLeg(1).GetMomentum(), GetLeg(2).GetMomentum(),
                                                    GetBJet(1).GetMomentum(), GetBJet(2).GetMomentum(),
                                                    GetMET(), energy_resolution_1, energy_resolution_2);
            kinfit_results->convergence = result.convergence;
            kinfit_results->chi2 = result.chi2;
            kinfit_results->probability = TMath::Prob(result.chi2, 2);
            kinfit_results->mass = result.mass;
        }
        else {
            const size_t index = static_cast<size_t>(std::distance(event->kinFit_jetPairId.begin(), iter));
            kinfit_results->convergence = event->kinFit_convergence.at(index);
            kinfit_results->chi2 = event->kinFit_chi2.at(index);
            kinfit_results->probability = TMath::Prob(kinfit_results->chi2, 2);
            kinfit_results->mass = event->kinFit_m.at(index);
        }
    }
    return *kinfit_results;
}

const sv_fit_ana::FitResults& EventInfoBase::GetSVFitResults()
{
    Lock lock(*mutex);
    if(!svfit_results){
        const auto iter = std::find(event->SVfit_Higges_indexes.begin(), event->SVfit_Higges_indexes.end(), selected_htt_index);
        svfit_results = std::make_shared<sv_fit_ana::FitResults>();
        if(iter == event->SVfit_Higges_indexes.end()){
            const auto& svfitProducer = GetSVFitProducer();
            const auto& result = svfitProducer.Fit(GetLeg(1),GetLeg(2),GetMET());
            svfit_results->has_valid_momentum = result.has_valid_momentum;
            svfit_results->momentum = result.momentum;
            svfit_results->momentum_error = result.momentum_error;
            svfit_results->transverseMass = result.transverseMass;
            svfit_results->transverseMass_error = result.transverseMass_error;
        }
        else {
            const size_t index = static_cast<size_t>(std::distance(event->SVfit_Higges_indexes.begin(), iter));
            svfit_results->has_valid_momentum = event->SVfit_is_valid.at(index);
            svfit_results->momentum = event->SVfit_p4.at(index);
            svfit_results->momentum_error = event->SVfit_p4_error.at(index);
            svfit_results->transverseMass = event->SVfit_mt.at(index);
            svfit_results->transverseMass_error = event->SVfit_mt_error.at(index);
        }
    }
    return *svfit_results;
}

LorentzVector EventInfoBase::GetResonanceMomentum(bool useSVfit, bool addMET)
{
    Lock lock(*mutex);
    if(useSVfit && addMET)
        throw exception("Can't add MET and with SVfit applied.");
    LorentzVector p4 = GetHiggsTTMomentum(useSVfit) + GetHiggsBB().GetMomentum();
    if(addMET)
        p4 += GetMET().GetMomentum();
    return p4;
}

double EventInfoBase::GetMT2()
{
    Lock lock(*mutex);
    if(!mt2.is_initialized()) {
        mt2 = Calculate_MT2(GetLeg(1).GetMomentum(), GetLeg(2).GetMomentum(),
                            GetHiggsBB().GetFirstDaughter().GetMomentum(),
                            GetHiggsBB().GetSecondDaughter().GetMomentum(), event->pfMET_p4);
    }
    return *mt2;
}

const FatJetCandidate* EventInfoBase::SelectFatJet(double mass_cut, double deltaR_subjet_cut)
{
    Lock lock(*mutex);
    using FatJet = ntuple::TupleFatJet;
    using SubJet = ntuple::TupleSubJet;
    if(!HasBjetPair()) return nullptr;
    for(const FatJetCandidate& fatJet : GetFatJets()) {
        if(fatJet->m(FatJet::MassType::SoftDrop) < mass_cut) continue;
        if(fatJet->subJets().size() < 2) continue;
        std::vector<SubJet> subJets = fatJet->subJets();
        std::sort(subJets.begin(), subJets.end(), [](const SubJet& j1, const SubJet& j2) -> bool {
            return j1.p4().Pt() > j2.p4().Pt(); });
        std::vector<double> deltaR;
        for(size_t n = 0; n < 2; ++n) {
            for(size_t k = 0; k < 2; ++k) {
                const auto dR = ROOT::Math::VectorUtil::DeltaR(subJets.at(n).p4(),
                                                               GetHiggsBB().GetDaughterMomentums().at(k));
                deltaR.push_back(dR);
            }
        }
        if((deltaR.at(0) < deltaR_subjet_cut && deltaR.at(3) < deltaR_subjet_cut)
                || (deltaR.at(1) < deltaR_subjet_cut && deltaR.at(2) < deltaR_subjet_cut))
            return &fatJet;
    }
    return nullptr;
}

void EventInfoBase::SetMvaScore(double _mva_score)
{
    Lock lock(*mutex);
    mva_score = _mva_score;
}

double EventInfoBase::GetMvaScore() const { return mva_score; }

boost::optional<EventInfoBase> CreateEventInfo(const ntuple::Event& event,
                                               const SignalObjectSelector& signalObjectSelector,
                                               const SummaryInfo* summaryInfo,
                                               Period period,
                                               JetOrdering jet_ordering)
{
    boost::optional<size_t> selected_higgs_index = signalObjectSelector.GetHiggsCandidateIndex(event);
    if(!selected_higgs_index.is_initialized()) return boost::optional<EventInfoBase>();
    SignalObjectSelector::SelectedSignalJets selected_signal_jets  = signalObjectSelector.SelectSignalJets(event,period,jet_ordering,*selected_higgs_index);
    return EventInfoBase(event,summaryInfo,*selected_higgs_index,selected_signal_jets,period,jet_ordering);

}



} // namespace analysis
