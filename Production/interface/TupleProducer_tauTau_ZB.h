/*! Definition of an event tuple producer for the tau-tau channel.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include <memory>
#include <vector>
#include "TTree.h"
#include "Math/VectorUtil.h"

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "h-tautau/Production/interface/BaseTupleProducer.h"

class TupleProducer_tauTau_ZB: public BaseTupleProducer {
public:
    using SelectionResults = analysis::SelectionResults<TauCandidate, TauCandidate>;
    using SelectionResultsPtr = std::shared_ptr<SelectionResults>;
    using HiggsCandidate = SelectionResults::HiggsCandidate;

public:
    TupleProducer_tauTau_ZB(const edm::ParameterSet& iConfig) : BaseTupleProducer(iConfig, "tauTau_ZB") {}

private:
    virtual void ProcessEvent(Cutter& cut) override;
    void FillEventTuple(const SelectionResults& selection);

    std::vector<TauCandidate> CollectSignalTaus();
    void SelectSignalTau(const TauCandidate& tau, Cutter& cut) const;
    static bool PtComparator(double pt_1, double pt_2);

private:
    SelectionResultsPtr previous_selection;
};
