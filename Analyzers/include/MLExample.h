#ifndef MLExample_h
#define MLExample_h

#include "AnalyzerCore.h"
#include "SystematicHelper.h"
#include "MLHelper.h"

class MLExample : public AnalyzerCore {
public:
    void initializeAnalyzer();
    void executeEvent();
    void TreeMaker();
    void Infer();


    bool RunTreeMake;
    bool RunInfer;

    TString IsoMuTriggerName;
    float TriggerSafePtCut;

    Event ev;
    RVec<TString> MuonIDSFKeys;
    RVec<Muon> AllMuons;
    RVec<Electron> AllElectrons;
    RVec<Jet> AllJets;

    RVec<Muon> Muons;
    RVec<Muon> Veto_Muons;
    RVec<Electron> Veto_Electrons;
    RVec<Jet> Jets;

    float weight_Prefire;

    unique_ptr<SystematicHelper> systHelper;
    unique_ptr<MLHelper> mlHelper_xgboost;
    unique_ptr<MLHelper> mlHelper_attention;

    FloatArray input_xgboost;
    FloatArray input_attention;


    MLExample();
    ~MLExample();
};

#endif

