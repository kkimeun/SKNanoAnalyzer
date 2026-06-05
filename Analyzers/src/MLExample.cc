#include "MLExample.h"

MLExample::MLExample() {}
MLExample::~MLExample() {}

void MLExample::initializeAnalyzer() {


    // Select the mode of the Analyzer, Using the user flags
    RunTreeMake = HasFlag("RunTreeMake");
    RunInfer = HasFlag("RunInfer");
    if(RunTreeMake && RunInfer){
        cerr << "[MLExample::initializeAnalyzer] RunTreeMake and RunInfer cannot be used together" << endl;
        exit(EXIT_FAILURE);
    }
    else if(!(RunTreeMake || RunInfer)){
        cerr << "[MLExample::initializeAnalyzer] RunTreeMake or RunInfer has to be used" << endl;
        exit(EXIT_FAILURE);
    }

    if(RunTreeMake){
        cout << "[MLExample::initializeAnalyzer] RunTreeMake" << endl;
    }
    else if(RunInfer){
        cout << "[MLExample::initializeAnalyzer] RunInfer" << endl;
    }

    if (DataEra == "2022EE") {
        IsoMuTriggerName = "HLT_IsoMu24";
        TriggerSafePtCut = 26.;
    }
    else if (DataEra == "2024") {
	IsoMuTriggerName = "HLT_IsoMu24";
	TriggerSafePtCut = 26.;
    }
    else {
        cerr << "[MLExample::initializeAnalyzer] Only 2022EE, 2024 are available" << endl;
        exit(EXIT_FAILURE);
    }

    cout << "[MLExample::initializeAnalyzer] IsoMuTriggerName = " << IsoMuTriggerName << endl;
    cout << "[MLExample::initializeAnalyzer] TriggerSafePtCut = " << TriggerSafePtCut << endl;

    // init B-Tagging (DeepJet Medium WP example)
    myCorr = new MyCorrection(DataEra, DataPeriod, IsDATA ? DataStream : MCSample, IsDATA);
    myCorr->SetTaggingParam(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerWP::Medium);

    // init SystematicHelper
    string SKNANO_HOME = getenv("SKNANO_HOME");
    if (IsDATA) {
        systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/docs/noSyst.yaml", DataStream, DataEra);
    } else {
        //systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/docs/ExampleSystematic.yaml", MCSample, DataEra);
        systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/docs/JesTotal.yaml", MCSample, DataEra);
    }
    
    
    if(RunTreeMake){
        // We can create a new Tree using NewTree function
        // Second argument is "keep", which is a list of branches in NanoAOD that you want to keep
        // Third argument is "drop", which is a list of branches in NanoAOD that you want to drop
        // Thus below line will create the empty tree
        NewTree("Training_Tree", {""}, {"*"});
    }
    else{
        string model_xgboost = "/data9/Users/eunsu/MLExample/Models/xgboost_model.onnx";
        string model_attention = "/data9/Users/eunsu/MLExample/Models/attention_model.onnx";
        if (model_xgboost == "" || model_attention == "") {
            throw runtime_error("Model file is not set");
        }
        mlHelper_xgboost = make_unique<MLHelper>(model_xgboost, MLHelper::ModelType::ONNX);
        mlHelper_attention = make_unique<MLHelper>(model_attention, MLHelper::ModelType::ONNX);
        input_xgboost.reserve(28);
        input_attention.reserve(28);
    }

    
}

void MLExample::executeEvent() {
    //==== Example 1
    //==== Dimuon Z-peak events with two muon IDs, with systematics
    
    // *IMPORTANT TO SAVE CPU TIME*
    // Every GetMuon() funtion first collect ALL NANOAOD muons with GetAllMuons(), and then check ID booleans.
    // GetAllMuons not only loops over all NANOAOD muons, but also actually CONSTRUCT muon objects for each muons.
    // We are now running systematics, and you don't want to do this for every systematic sources
    // So, I defined "RVec<Muon> AllMuons;" in Analyzers/include/MLExample.h,
    // and save muons objects at the very beginning of executeEvent().
    // Later, do "SelectMuons(AllMuons, ID, pt, eta)" to get muons with ID cuts
    AllMuons = GetAllMuons();
    AllElectrons = GetAllElectrons();
    AllJets = GetAllJets();
    ev = GetEvent();
    
    //Below are the baseline selections.
    Veto_Muons = SelectMuons(AllMuons, Muon::MuonID::POG_MEDIUM, 15., 2.4);
    Muons = SelectMuons(Veto_Muons, Muon::MuonID::POG_MEDIUM, 26., 2.4);
    Veto_Electrons = SelectElectrons(AllElectrons, Electron::ElectronID::POG_MVAISO_WP90, 15., 2.4);
    Jets = SelectJets(AllJets, Jet::JetID::TIGHT, 30., 2.4);
    //Require OS dimuon pair
    if(!(Muons.size() == 2 && Veto_Muons.size() == 2 && Veto_Electrons.size() == 0)) return;
    if(Muons.at(0).Charge() == Muons.at(1).Charge()) return;
    // Do lepton Cleaning
    Jets = JetsVetoLeptonInside(Jets, Veto_Electrons, Muons, 0.4);
    // Require at least 2 jets
    if(Jets.size() < 2) return;
    
    // Above selections allow to Dileptonic ttbar events and Z->mumu+Jets events
    // Of course these two categories have very easy even if we do not use ML, will clear differences in 
    // M_{ll}, b-jet multiplicity, MET.
    // But this is just an example.

    if(RunTreeMake){
        // RunTreeMaker() is a function that fills the tree
        TreeMaker();
    }
    else{
        // RunInfer() is a function that fills the Infer
        input_xgboost.clear();
        input_attention.clear();
        Infer();
    }
       
}

void MLExample::TreeMaker(){
    // In this function, Fill the tree with variables you want to use for ML training.
    // Let's use fixed size of Objects for simplicity.
    // If you want to use variable size of objects, It will requires more complicated techniques like Masking.
    // First sort Jets and Muons by Pt
    std::sort(Jets.begin(), Jets.end(), PtComparing);
    std::sort(Muons.begin(), Muons.end(), PtComparing);
    JetTagging::JetFlavTagger FlavTagger = JetTagging::JetFlavTagger::ParT;
    // Fill the tree
    // SetBranch function will get the value, and internaly save these values in the deques.
    // And link the address of each element of the deque to the branch of the tree.
    // For now, only float and int, bool are supported.
    SetBranch("Training_Tree", "Muon0_Pt", Muons.at(0).Pt());
    SetBranch("Training_Tree", "Muon0_Eta", Muons.at(0).Eta());
    SetBranch("Training_Tree", "Muon0_CosPhi", std::cos(Muons.at(0).Phi()));
    SetBranch("Training_Tree", "Muon0_SinPhi", std::sin(Muons.at(0).Phi()));
    SetBranch("Training_Tree", "Muon0_M", Muons.at(0).M());
    SetBranch("Training_Tree", "Muon0_BvsC", 0.f);
    SetBranch("Training_Tree", "Muon0_isJet", 0);
    SetBranch("Training_Tree", "Muon1_Pt", Muons.at(1).Pt());
    SetBranch("Training_Tree", "Muon1_Eta", Muons.at(1).Eta());
    SetBranch("Training_Tree", "Muon1_CosPhi", std::cos(Muons.at(1).Phi()));
    SetBranch("Training_Tree", "Muon1_SinPhi", std::sin(Muons.at(1).Phi()));
    SetBranch("Training_Tree", "Muon1_M", Muons.at(1).M());
    SetBranch("Training_Tree", "Muon1_BvsC", 0.f);
    SetBranch("Training_Tree", "Muon1_isJet", 0);
    SetBranch("Training_Tree", "Jet0_Pt", Jets.at(0).Pt());
    SetBranch("Training_Tree", "Jet0_Eta", Jets.at(0).Eta());
    SetBranch("Training_Tree", "Jet0_CosPhi", std::cos(Jets.at(0).Phi()));
    SetBranch("Training_Tree", "Jet0_SinPhi", std::sin(Jets.at(0).Phi()));
    SetBranch("Training_Tree", "Jet0_M", Jets.at(0).M());
    SetBranch("Training_Tree", "Jet0_BvsC", Jets.at(0).GetTaggerResult(FlavTagger, JetTagging::JetFlavTaggerScoreType::B));
    SetBranch("Training_Tree", "Jet0_isJet", 1);
    SetBranch("Training_Tree", "Jet1_Pt", Jets.at(1).Pt());
    SetBranch("Training_Tree", "Jet1_Eta", Jets.at(1).Eta());
    SetBranch("Training_Tree", "Jet1_CosPhi", std::cos(Jets.at(1).Phi()));
    SetBranch("Training_Tree", "Jet1_SinPhi", std::sin(Jets.at(1).Phi()));
    SetBranch("Training_Tree", "Jet1_M", Jets.at(1).M());
    SetBranch("Training_Tree", "Jet1_BvsC", Jets.at(1).GetTaggerResult(FlavTagger, JetTagging::JetFlavTaggerScoreType::B));
    SetBranch("Training_Tree", "Jet1_isJet", 1);

    // After setting all the branches, Fill the tree
    // if you pass no argument to FillTrees, it will fill all trees that you have created.
    FillTrees("Training_Tree");
    //FillTree Function will Fill and clear the deques.
    //So be careful to this function not making an unintended behavior.
    //You might need placeholders if You only fill some trees in the event.
}


void MLExample::Infer(){
    std::unordered_map<std::string, std::vector<int>> input_shape_xgboost;
    // input_shape_xgboost["input"] = {1, 24}; // batch size, input shape
    // Order must be same as the order of the input nodes in the model
    /* Jet0_BvsC', 'Jet0_Eta', 'Jet0_M', 'Jet0_Phi', 'Jet0_Pt', 'Jet0_isJet',
       'Jet1_BvsC', 'Jet1_Eta', 'Jet1_M', 'Jet1_Phi', 'Jet1_Pt', 'Jet1_isJet',
       'Muon0_BvsC', 'Muon0_Eta', 'Muon0_M', 'Muon0_Phi', 'Muon0_Pt',
       'Muon0_isJet', 'Muon1_BvsC', 'Muon1_Eta', 'Muon1_M', 'Muon1_Phi',
       'Muon1_Pt', 'Muon1_isJet'*/
    // input_xgboost = {Jets.at(0).GetTaggerResult(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerScoreType::B), Jets.at(0).Eta(), Jets.at(0).M(), Jets.at(0).Phi(), Jets.at(0).Pt(), 1,
    //                 Jets.at(1).GetTaggerResult(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerScoreType::B), Jets.at(1).Eta(), Jets.at(1).M(), Jets.at(1).Phi(), Jets.at(1).Pt(), 1,
    //                 0.f, Muons.at(0).Eta(), Muons.at(0).M(), Muons.at(0).Phi(), Muons.at(0).Pt(), 0,
    //                 0.f, Muons.at(1).Eta(), Muons.at(1).M(), Muons.at(1).Phi(), Muons.at(1).Pt(), 0};

    input_shape_xgboost["input"] = {1, 28};

    input_xgboost = {
        // Jet0
        Jets.at(0).GetTaggerResult(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerScoreType::B),
        float(std::cos(Jets.at(0).Phi())),
        float(Jets.at(0).Eta()),
        float(Jets.at(0).M()),
        float(Jets.at(0).Pt()),
        float(std::sin(Jets.at(0).Phi())),
        1.f,

        // Jet1
        Jets.at(1).GetTaggerResult(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerScoreType::B),
        float(std::cos(Jets.at(1).Phi())),
        float(Jets.at(1).Eta()),
        float(Jets.at(1).M()),
        float(Jets.at(1).Pt()),
        float(std::sin(Jets.at(1).Phi())),
        1.f,

        // Muon0
        0.f,
        float(std::cos(Muons.at(0).Phi())),
        float(Muons.at(0).Eta()),
        float(Muons.at(0).M()),
        float(Muons.at(0).Pt()),
        float(std::sin(Muons.at(0).Phi())),
        0.f,

        // Muon1
        0.f,
        float(std::cos(Muons.at(1).Phi())),
        float(Muons.at(1).Eta()),
        float(Muons.at(1).M()),
        float(Muons.at(1).Pt()),
        float(std::sin(Muons.at(1).Phi())),
        0.f
    };


    // Run the model
    std::unordered_map<std::string, VariousArray> input_data_xgboost;
    // VariousArray is a union of FloatArray, IntArray, and BoolArray(uint8_t)
    input_data_xgboost["input"] = input_xgboost;
    std::unordered_map<std::string, FloatArray> output_xgboost = mlHelper_xgboost->Run_ONNX_Model(input_data_xgboost, input_shape_xgboost);

    std::unordered_map<std::string, std::vector<int>> input_shape_attention;
    input_shape_attention["input"] = {1, 4, 7}; // batch size, input shape
    
    // This model is symmetric under the permutation of the objects, but the order of features should be maintained
    // Though the input dimension is 2D, We must flatten the input as 1D array with row-major order
    // Order of features: '_M', '_SinPhi', '_CosPhi', '_Pt', '_BvsC', '_Eta', '_isJet'
    //input_attention = {Jets.at(0).GetTaggerResult(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerScoreType::B), Jets.at(0).Phi(), Jets.at(0).Eta(), Jets.at(0).M(), Jets.at(0).Pt(), 1,
    //                   Jets.at(1).GetTaggerResult(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerScoreType::B), Jets.at(1).Phi(), Jets.at(1).Eta(), Jets.at(1).M(), Jets.at(1).Pt(), 1,
    //                   0.f, Muons.at(0).Phi(), Muons.at(0).Eta(), Muons.at(0).M(), Muons.at(0).Pt(), 0,
    //                   0.f, Muons.at(1).Phi(), Muons.at(1).Eta(), Muons.at(1).M(), Muons.at(1).Pt(), 0};

    input_attention = {
        // Jet0: BvsC, CosPhi, Eta, M, Pt, SinPhi, isJet
        Jets.at(0).GetTaggerResult(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerScoreType::B),
        float(std::cos(Jets.at(0).Phi())),
        float(Jets.at(0).Eta()),
        float(Jets.at(0).M()),
        float(Jets.at(0).Pt()),
        float(std::sin(Jets.at(0).Phi())),
        1.f,

        // Jet1
        Jets.at(1).GetTaggerResult(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerScoreType::B),
        float(std::cos(Jets.at(1).Phi())),
        float(Jets.at(1).Eta()),
        float(Jets.at(1).M()),
        float(Jets.at(1).Pt()),
        float(std::sin(Jets.at(1).Phi())),
        1.f,

        // Muon0
        0.f,
        float(std::cos(Muons.at(0).Phi())),
        float(Muons.at(0).Eta()),
        float(Muons.at(0).M()),
        float(Muons.at(0).Pt()),
        float(std::sin(Muons.at(0).Phi())),
        0.f,

        // Muon1
        0.f,
        float(std::cos(Muons.at(1).Phi())),
        float(Muons.at(1).Eta()),
        float(Muons.at(1).M()),
        float(Muons.at(1).Pt()),
        float(std::sin(Muons.at(1).Phi())),
        0.f
    };

    cout << "[MLExample::Infer] input_xgboost.size() = " << input_xgboost.size() << endl;
    cout << "[MLExample::Infer] input_attention.size() = " << input_attention.size() << endl;

    // Run the model
    std::unordered_map<std::string, VariousArray> input_data_attention;
    input_data_attention["input"] = input_attention;
    std::unordered_map<std::string, FloatArray> output_attention = mlHelper_attention->Run_ONNX_Model(input_data_attention, input_shape_attention);

    // Do something with the output. name of the output is can be checked by when you save the model.
    float xgboost_output = output_xgboost.at("probabilities")[1];
    FloatArray attention_output = output_attention.at("output");
    // for attention model, we have to do softmax by ourselves
    if (!attention_output.empty()) {
        std::vector<float> softmax_output;
        for(auto val : attention_output) {
            softmax_output.push_back(std::exp(val));
        }
        
        float sum = std::accumulate(softmax_output.begin(), softmax_output.end(), 0.0f);
        if (sum > 0) {
            std::transform(softmax_output.begin(), softmax_output.end(), softmax_output.begin(), 
                            [sum](float val) { return val / sum; });
        }

        // Histogram filling with boundary checks
        if (softmax_output.size() > 1) {
            FillHist("xgboost_output", xgboost_output, 1., 50, 0., 1.);
            FillHist("attention_output", softmax_output[1], 1., 50, 0., 1.);
        }
    }
    // Fill the histogram
    FillHist("Jet0_BvsC", Jets.at(0).GetTaggerResult(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerScoreType::B), 1., 50, 0., 1.);
    FillHist("Jet1_BvsC", Jets.at(1).GetTaggerResult(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerScoreType::B), 1., 50, 0., 1.);
    FillHist("Muon0_Pt", Muons.at(0).Pt(), 1., 50, 0., 200.);
    FillHist("Muon1_Pt", Muons.at(1).Pt(), 1., 50, 0., 200.);
    FillHist("DiMuon_M", (Muons.at(0) + Muons.at(1)).M(), 1., 50, 0., 200.);
}
