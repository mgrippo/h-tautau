/*! Classes that represent analysis objects at the tuple level.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisMath.h"
#include "AnalysisTypes.h"
#include "EventTuple.h"

namespace ntuple {

class TupleObject {
public:
    using MetType = analysis::MetType;
    using DiscriminatorWP = analysis::DiscriminatorWP;
    using LorentzVectorE = analysis::LorentzVector;
    using LorentzVectorM = analysis::LorentzVectorM;
    using DiscriminatorResult = float;

    TupleObject(const ntuple::Event& _event) : event(&_event) {}

protected:
    const Event* event;
};

class TupleLepton : public TupleObject {
public:
    TupleLepton(const ntuple::Event& _event, size_t _leg_id)
        : TupleObject(_event), leg_id(_leg_id)
    {
        if(leg_id < 1 || leg_id > 2)
            throw analysis::exception("Invalid leg id = %1%.") % leg_id;
    }

    const LorentzVectorM& p4() const { return leg_id == 1 ? event->p4_1 : event->p4_2; }
    int charge() const { return leg_id == 1 ? event->q_1 : event->q_2; }
    double d0() const { return leg_id == 1 ? event->d0_1 : event->d0_2; }
    double dZ() const { return leg_id == 1 ? event->dZ_1 : event->dZ_2; }

    double mt(MetType metType = MetType::PF) const
    {
        if(metType == MetType::PF)
            return leg_id == 1 ? event->pfmt_1 : event->pfmt_2;
        if(metType == MetType::MVA)
            return leg_id == 1 ? event->mt_1 : event->mt_2;
        if(metType == MetType::PUPPI)
            return leg_id == 1 ? event->puppimt_1 : event->puppimt_2;
        throw analysis::exception("Unsupported MET type.");
    }

    double iso() const { return leg_id == 1 ? event->iso_1 : event->iso_2; }
    int gen_match() const { return leg_id == 1 ? event->gen_match_1 : event->gen_match_2; }

protected:
    size_t leg_id;
};

class TupleElectron : public TupleLepton {
public:
    explicit TupleElectron(const ntuple::Event& _event, size_t _leg_id = 1) : TupleLepton(_event, _leg_id) {}
};

class TupleMuon : public TupleLepton {
public:
    explicit TupleMuon(const ntuple::Event& _event, size_t _leg_id = 1) : TupleLepton(_event, _leg_id) {}
};

class TupleTau : public TupleLepton {
public:
    using TupleLepton::TupleLepton;

    DiscriminatorResult tauID(const std::string& discriminator) const
    {
        const auto& tauIDs = leg_id == 1 ? event->tauIDs_1 : event->tauIDs_2;
        if(!tauIDs.count(discriminator))
            throw analysis::exception("TauID discriminator '%1%' not found.") % discriminator;
        return tauIDs.at(discriminator);
    }

    DiscriminatorResult againstElectronMVA6(DiscriminatorWP wp) const
    {
        std::ostringstream ss_name;
        ss_name << "againstElectron" << wp << "MVA6";
        return tauID(ss_name.str());
    }

    DiscriminatorResult againstMuon3(DiscriminatorWP wp) const
    {
        std::ostringstream ss_name;
        ss_name << "againstMuon" << wp << "3";
        return tauID(ss_name.str());
    }

    DiscriminatorResult byCombinedIsolationDeltaBetaCorrRaw3Hits() const
    {
        return tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
    }

    DiscriminatorResult byIsolationMVA3raw(bool use_new_dm, bool use_lifetime) const
    {
        const std::string dm_str = use_new_dm ? "new" : "old";
        const std::string lt_str = use_lifetime ? "w" : "wo";
        std::ostringstream ss_name;
        ss_name << "byIsolationMVArun2v1DB" << dm_str << "DM" << lt_str << "LTraw";
        return tauID(ss_name.str());
    }
};

class TupleJet : public TupleObject {
public:
    TupleJet(const ntuple::Event& _event, size_t _jet_id)
        : TupleObject(_event), jet_id(_jet_id)
    {
        if(jet_id >= event->jets_p4.size())
            throw analysis::exception("Jet id = %1% is out of range.") % jet_id;
    }

    const LorentzVectorE& p4() const { return event->jets_p4.at(jet_id); }
    DiscriminatorResult mva() const { return event->jets_mva.at(jet_id); }
    DiscriminatorResult csv() const { return event->jets_csv.at(jet_id); }
    DiscriminatorResult partonFlavour() const { return event->jets_partonFlavour.at(jet_id); }

private:
    size_t jet_id;
};

class TupleSubJet : public TupleObject {
public:
    TupleSubJet(const ntuple::Event& _event, size_t _jet_id)
        : TupleObject(_event), jet_id(_jet_id)
    {
        if(jet_id >= event->subJets_p4.size())
            throw analysis::exception("Fat sub-jet id = %1% is out of range.") % jet_id;
    }

    const LorentzVectorE& p4() const { return event->subJets_p4.at(jet_id); }
    DiscriminatorResult csv() const { return event->subJets_csv.at(jet_id); }

private:
    size_t jet_id;
};

class TupleFatJet : public TupleObject {
public:
    enum class MassType { Pruned, Filtered, Trimmed, SoftDrop };

    TupleFatJet(const ntuple::Event& _event, size_t _jet_id)
        : TupleObject(_event), jet_id(_jet_id)
    {
        if(jet_id >= event->fatJets_p4.size())
            throw analysis::exception("Fat jet id = %1% is out of range.") % jet_id;

        for(size_t n = 0; n < event->subJets_p4.size(); ++n) {
            if(event->subJets_parentIndex.at(n) == jet_id)
                sub_jets.push_back(TupleSubJet(_event, n));
        }
    }

    const LorentzVectorE& p4() const { return event->fatJets_p4.at(jet_id); }
    DiscriminatorResult csv() const { return event->fatJets_csv.at(jet_id); }

    double m(MassType massType) const
    {
        if(massType == MassType::Pruned) return event->fatJets_m_pruned.at(jet_id);
        if(massType == MassType::Filtered) return event->fatJets_m_filtered.at(jet_id);
        if(massType == MassType::Trimmed) return event->fatJets_m_trimmed.at(jet_id);
        if(massType == MassType::SoftDrop) return event->fatJets_m_softDrop.at(jet_id);
        throw analysis::exception("Unsupported fat jet mass type");
    }

    DiscriminatorResult n_subjettiness(size_t tau_index) const
    {
        if(tau_index == 1) return event->fatJets_n_subjettiness_tau1.at(jet_id);
        if(tau_index == 2) return event->fatJets_n_subjettiness_tau2.at(jet_id);
        if(tau_index == 3) return event->fatJets_n_subjettiness_tau3.at(jet_id);
        throw analysis::exception("Unsupported tau index = %1% for fat jet subjettiness.") % tau_index;
    }

    const std::vector<TupleSubJet>& subJets() const { return sub_jets; }

private:
    size_t jet_id;
    std::vector<TupleSubJet> sub_jets;
};

class TupleMet : public TupleObject {
public:
    using CovMatrix = analysis::SquareMatrix<2>;
    TupleMet(const ntuple::Event& _event, MetType _met_type)
        : TupleObject(_event), met_type(_met_type)
    {
        static const std::set<MetType> supported_types = { MetType::PF, MetType::MVA, MetType::PUPPI };
        if(!supported_types.count(met_type))
            throw analysis::exception("Unsupported met type.");
    }

    MetType type() const { return met_type; }

    const LorentzVectorM& p4() const
    {
        if(met_type == MetType::PF) return event->pfMET_p4;
        if(met_type == MetType::MVA) return event->mvaMET_p4;
        return event->puppiMET_p4;
    }

    const CovMatrix& cov() const
    {
        if(met_type == MetType::PF) return event->pfMET_cov;
        if(met_type == MetType::MVA) return event->mvaMET_cov;
        return event->puppiMET_cov;
    }

    double pt() const { return p4().pt(); }
    double phi() const { return p4().phi(); }

private:
    MetType met_type;
};

} // namespace ntuple
