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

class TupleProducer_tauTau: public BaseTupleProducer {
public:
    using SelectionResults = analysis::SelectionResults<TauCandidate>;
    using HiggsCandidate = SelectionResults::HiggsCandidate;

public:
    TupleProducer_tauTau(const edm::ParameterSet& iConfig) : BaseTupleProducer(iConfig, "tauTau") {}

private:
    virtual void ProcessEvent(Cutter& cut) override;

    std::vector<TauCandidate> CollectSignalTaus();
    void SelectSignalTau(const TauCandidate& tau, Cutter& cut) const;
    void FillEventTuple(const SelectionResults& selection);
    void FillSyncTuple(const SelectionResults& selection);

};
