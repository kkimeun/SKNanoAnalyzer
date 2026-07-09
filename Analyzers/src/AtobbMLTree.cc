#include "AtobbMLTree.h"

AtobbMLTree::AtobbMLTree() :
  const_top_mass(172.5),
  const_top_width(1.5),
  const_w_mass(80.4),
  const_w_width(2.085)
{
}

AtobbMLTree::~AtobbMLTree() {}

void AtobbMLTree::initializeAnalyzer() {

  MuonIDs = {
    Muon::MuonID::POG_TIGHT,
    Muon::MuonID::POG_MEDIUM_PROMPT,
    Muon::MuonID::POG_MVA_MU_TIGHT
  };

  MuonIDISOSFKeys = {
    "NUM_TightID_DEN_TrackerMuons",
    "NUM_TightPFIso_DEN_TightID",
    "NUM_MediumPromptID_DEN_TrackerMuons",
    "NUM_TightPFIso_DEN_MediumPromptID",
    "NUM_TightMvaMuID_DEN_TrackerMuons"
  };

  if (DataEra == "2016preVFP" || DataEra == "2016postVFP" || DataEra == "2018") {
    IsoMuTriggerName = "HLT_IsoMu24";
    TriggerSafePtCut = 26.;
  }
  else if (DataEra == "2017") {
    IsoMuTriggerName = "HLT_IsoMu27";
    TriggerSafePtCut = 29.;
  }
  else if (DataEra == "2022" || DataEra == "2022EE") {
    IsoMuTriggerName = "HLT_IsoMu24";
    TriggerSafePtCut = 26.;
  }
  else if (DataEra == "2023") {
    IsoMuTriggerName = "";
    TriggerSafePtCut = 26.;
  }
  else if (DataEra == "2024") {
    IsoMuTriggerName = "HLT_IsoMu24";
    TriggerSafePtCut = 26.;
  }
  else {
    cerr << "[AtobbMLTree::initializeAnalyzer] DataEra is not set properly" << endl;
    exit(EXIT_FAILURE);
  }

  cout << "[AtobbMLTree::initializeAnalyzer] IsoMuTriggerName = "
       << IsoMuTriggerName << endl;
  cout << "[AtobbMLTree::initializeAnalyzer] TriggerSafePtCut = "
       << TriggerSafePtCut << endl;

  myCorr = new MyCorrection(DataEra, DataPeriod, IsDATA ? DataStream : MCSample, IsDATA);
  myCorr->SetTaggingParam(
      JetTagging::JetFlavTagger::ParT,
      JetTagging::JetFlavTaggerWP::Medium
  );

  string SKNANO_HOME = getenv("SKNANO_HOME");
  if (IsDATA) {
    systHelper = std::make_unique<SystematicHelper>(
        SKNANO_HOME + "/docs/noSyst.yaml", DataStream, DataEra
    );
  }
  else {
    systHelper = std::make_unique<SystematicHelper>(
        SKNANO_HOME + "/docs/noSyst.yaml", MCSample, DataEra
    );
  }

  NewTree("Training_Tree", {""}, {"*"});
}

void AtobbMLTree::executeEvent() {

  AllMuonViews = GetAllMuonViews();
  AllElectronViews = GetAllElectronViews();
  AllJetViews = GetAllJetViews();
  AllGenViews = GetAllGenViews();
  ev = GetEvent();

  if (!PassMetFilter(AllJetViews, ev)) return;
  if (IsoMuTriggerName != "" && !(ev.PassTrigger(IsoMuTriggerName))) return;

  Particle METv = ev.GetMETVector(Event::MET_Type::PUPPI);

  bool use_pog_tight_muon_id = true;
  bool use_pog_mva_tight_muon_id = false;
  bool correct_b_jet_pt = true;
  bool use_UParT_JEC = false;
  bool apply_pu_id = false;

  Muon::MuonID this_muon_id = MuonIDs[0];
  TString this_muon_id_sf_key = MuonIDISOSFKeys[0];
  TString this_muon_iso_sf_key = MuonIDISOSFKeys[1];
  TString this_muon_trig_sf_key = "";

  if (use_pog_mva_tight_muon_id) {
    this_muon_id = MuonIDs[2];
    this_muon_id_sf_key = MuonIDISOSFKeys[4];
  }
  else if (use_pog_tight_muon_id) {
    this_muon_id = MuonIDs[0];
    this_muon_id_sf_key = MuonIDISOSFKeys[0];
    this_muon_iso_sf_key = MuonIDISOSFKeys[1];
    this_muon_trig_sf_key = "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight";
  }
  else {
    this_muon_id = MuonIDs[1];
    this_muon_id_sf_key = MuonIDISOSFKeys[2];
    this_muon_iso_sf_key = MuonIDISOSFKeys[3];
  }

  std::vector<size_t> SelectedMuonIndices_id_only =
      SelectMuonIndices(AllMuonViews, this_muon_id, 15., 2.4);

  std::vector<size_t> SelectedMuonIndices = {};
  if (use_pog_mva_tight_muon_id) {
    SelectedMuonIndices = SelectedMuonIndices_id_only;
  }
  else {
    SelectedMuonIndices = SelectMuonIndices(
        AllMuonViews,
        SelectedMuonIndices_id_only,
        Muon::MuonID::POG_PFISO_TIGHT,
        15.,
        2.4
    );
  }

  std::vector<size_t> SelectedElectronIndices =
      SelectElectronIndices(AllElectronViews, Electron::ElectronID::POG_LOOSE, 15., 2.5);

  if (SelectedMuonIndices.size() + SelectedElectronIndices.size() != 1) return;

  RVec<Muon> muons = MaterializeMuons(AllMuonViews, SelectedMuonIndices);
  RVec<Electron> electrons = MaterializeElectrons(AllElectronViews, SelectedElectronIndices);

  Lepton lepton;
  if (muons.size() == 1) {
    lepton = Lepton(muons.at(0));
    if (muons.at(0).Pt() <= TriggerSafePtCut) return;
  } else if (electrons.size() == 1) {
    lepton = Lepton(electrons.at(0));
    if (electrons.at(0).Pt() <= 32) return; //TODO: set electron trigger safe cut properly
  } else {
    return;
  }

//  if (muons.size() != 1) return;
//  if (electrons.size() != 0) return;
//  if (muons.at(0).Pt() <= TriggerSafePtCut) return;

  MyCorrection::variation jes_variation = MyCorrection::variation::nom;
  MyCorrection::variation btag_jes_variation = MyCorrection::variation::nom;
  TString btag_source = "total";

  auto jet_id = apply_pu_id ? Jet::JetID::PUID_LOOSE : Jet::JetID::TIGHT;

  std::vector<size_t> SelectedJetIndices =
      SelectJetIndices(AllJetViews, jet_id, 0., 5.191, jes_variation, MyCorrection::variation::nom);

  RVec<Jet> jets =
      MaterializeJets(AllJetViews, SelectedJetIndices, jes_variation, MyCorrection::variation::nom);

  jets = JetsVetoLeptonInside(jets, electrons, muons, 0.3);

  // std::sort(muons.begin(), muons.end(), PtComparing);
  std::sort(jets.begin(), jets.end(), PtComparing);

  if (jets.size() < 6) return;
  if (!PassJetVetoMap(AllJetViews, AllMuonViews, "jetvetomap_fpix")) return;

  float btag_wp_cut = myCorr->GetBTaggingWP();

  for (auto& jet : jets) {
    double bscore = jet.GetTaggerResult(
        JetTagging::JetFlavTagger::ParT,
        JetTagging::JetFlavTaggerScoreType::B
    );

    if (bscore > btag_wp_cut) {
      double current_pt  = jet.Pt();
      double raw_pt      = jet.GetRawPt();
      double current_eta = jet.Eta();
      double current_phi = jet.Phi();
      double current_m   = jet.M();

      double L1L2L3Res = myCorr->GetJESSF(
          jet.GetArea(),
          current_eta,
          raw_pt,
          current_phi,
          Rho_fixedGridRhoFastjetAll,
          ev.run()
      );

      double L2L3Residual = myCorr->GetJESSF(
          jet.GetArea(),
          current_eta,
          raw_pt,
          current_phi,
          Rho_fixedGridRhoFastjetAll,
          ev.run(),
          "L2L3Residual"
      );

      double UParTAK4RegPtRawCorr = jet.UParTAK4RegPtRawCorr();
      double UParTAK4RegPtRawCorrNeutrino = jet.UParTAK4RegPtRawCorrNeutrino();
      double UParT_ratio = UParTAK4RegPtRawCorrNeutrino / UParTAK4RegPtRawCorr;

      double modified_pt = current_pt;
      double modified_m  = current_m;

      if (use_UParT_JEC) {
        modified_pt = current_pt * (1. / L1L2L3Res)
                    * UParTAK4RegPtRawCorrNeutrino
                    * L2L3Residual;
        modified_m = current_m * modified_pt / current_pt;
      }
      else {
        modified_pt = current_pt * UParT_ratio;
        modified_m  = current_m  * UParT_ratio;
      }

      if (correct_b_jet_pt) {
        jet.SetPtEtaPhiM(modified_pt, current_eta, current_phi, modified_m);
      }
    }
  }

  std::vector<bool> btag_vector;
  int n_bjets = 0;
  int n_jets_for_count = 0;

  for (unsigned int ij = 0; ij < jets.size(); ij++) {
    double bscore = jets.at(ij).GetTaggerResult(
        JetTagging::JetFlavTagger::ParT,
        JetTagging::JetFlavTaggerScoreType::B
    );

    bool is_btag = bscore > btag_wp_cut;
    btag_vector.push_back(is_btag);

    if (is_btag) {
      n_bjets++;
      if (jets.at(ij).Pt() > 25.) n_jets_for_count++;
    }
    else {
      if (jets.at(ij).Pt() > 30.) n_jets_for_count++;
    }
  }

  if (n_bjets < 3) return;
  if (n_jets_for_count < 6) return;

  float weight_mc = 1.f;
  float weight_trigger_lumi = 1.f;
  float weight_pu = 1.f;
  float weight_btag = 1.f;
  float weight_muon_id = 1.f;
  float weight_muon_iso = 1.f;
  float weight_muon_trig = 1.f;

  if (!IsDATA) {
    weight_mc = MCweight();
    weight_trigger_lumi = ev.GetTriggerLumi("Full");

    weight_muon_id = myCorr->GetMuonIDSF(
        this_muon_id_sf_key,
        muons,
        MyCorrection::variation::nom
    );

    if (!use_pog_mva_tight_muon_id) {
      weight_muon_iso = myCorr->GetMuonIDSF(
          this_muon_iso_sf_key,
          muons,
          MyCorrection::variation::nom
      );
    }

    if (use_pog_tight_muon_id) {
      weight_muon_trig = myCorr->GetMuonTriggerSF(
          this_muon_trig_sf_key,
          muons,
          MyCorrection::variation::nom
      );
    }

    weight_pu = myCorr->GetPUWeight(ev.nTrueInt(), MyCorrection::variation::nom);

    RVec<Jet> jets_2p5;
    for (auto jet : jets) {
      if (std::abs(jet.Eta()) < 2.499 && jet.Pt() > 20.) {
        jets_2p5.push_back(jet);
      }
    }

    weight_btag = myCorr->GetBTaggingSF(
        jets_2p5,
        JetTagging::JetFlavTagger::ParT,
        JetTagging::JetFlavTaggerWP::Medium,
        JetTagging::JetTaggingSFMethod::comb,
        btag_jes_variation,
        btag_source
    );
  }

  float weight_train = 1.f;
  if (!IsDATA) {
    weight_train = weight_mc * weight_trigger_lumi * weight_pu * weight_btag;
  }

  float weight_nominal = 1.f;
  if (!IsDATA) {
    weight_nominal = weight_mc
                   * weight_trigger_lumi
                   * weight_pu
                   * weight_btag
                   * weight_muon_id
                   * weight_muon_iso
                   * weight_muon_trig;
  }

  float best_chi2 = -999.f;
  float had_W_mass = -999.f;
  float had_top_mass = -999.f;
  float lep_W_mass = -999.f;
  float lep_top_mass = -999.f;

  std::vector<unsigned int> top_b_jet_candidates;
  std::vector<unsigned int> had_W_candidates;

  unsigned int n_check = std::min<unsigned int>(6, jets.size());

  for (unsigned int ij = 0; ij < n_check; ij++) {
    if (btag_vector.at(ij)) {
      if (top_b_jet_candidates.size() < 2) {
        top_b_jet_candidates.push_back(ij);
      }
    }
    else {
      if (had_W_candidates.size() < 3) {
        had_W_candidates.push_back(ij);
      }
    }
  }

  if (top_b_jet_candidates.size() >= 2 && had_W_candidates.size() >= 2) {
    std::vector<ttCombinatoric> combinatorics;

    for (unsigned int w1 = 0; w1 < had_W_candidates.size() - 1; w1++) {
      for (unsigned int w2 = w1 + 1; w2 < had_W_candidates.size(); w2++) {
        for (unsigned int b_idx = 0; b_idx < 2; b_idx++) {
          ttCombinatoric comb;
          comb.lepton = &lepton;
          comb.jets = &jets;
          comb.met = &METv;
          comb.had_W_jet_idx_1 = had_W_candidates.at(w1);
          comb.had_W_jet_idx_2 = had_W_candidates.at(w2);
          comb.had_top_b_jet_idx = top_b_jet_candidates.at(b_idx);
          comb.lep_top_b_jet_idx = top_b_jet_candidates.at(1 - b_idx);

          this->EvalChi2(comb);
          combinatorics.push_back(comb);
        }
      }
    }

    if (!combinatorics.empty()) {
      ttCombinatoric* best = &combinatorics.at(0);
      for (unsigned int i = 1; i < combinatorics.size(); i++) {
        if (combinatorics.at(i).best_chi2 < best->best_chi2) {
          best = &combinatorics.at(i);
        }
      }

      best_chi2 = best->best_chi2;
      had_W_mass = best->had_W_mass;
      had_top_mass = best->had_top_mass;
      lep_W_mass = best->best_lep_W_mass;
      lep_top_mass = best->best_lep_top_mass;
    }
  }

  

  FillTreeBranches(
      // muons,
      lepton,
      jets,
      METv,
      btag_vector,
      n_bjets,
      n_jets_for_count,
      btag_wp_cut,
      weight_train,
      weight_nominal,
      weight_mc * weight_trigger_lumi,
      weight_pu,
      weight_btag,
      weight_muon_id,
      weight_muon_iso,
      weight_muon_trig,
      best_chi2,
      had_W_mass,
      had_top_mass,
      lep_W_mass,
      lep_top_mass
  );
}


void AtobbMLTree::FillJetBranches(const TString& prefix, const Jet& jet, float btag_wp_cut) {
  float bscore = jet.GetTaggerResult(
      JetTagging::JetFlavTagger::ParT,
      JetTagging::JetFlavTaggerScoreType::B
  );

  SetBranch("Training_Tree", prefix + "_Pt",     float(jet.Pt()));
  SetBranch("Training_Tree", prefix + "_Eta",    float(jet.Eta()));
  SetBranch("Training_Tree", prefix + "_CosPhi", float(std::cos(jet.Phi())));
  SetBranch("Training_Tree", prefix + "_SinPhi", float(std::sin(jet.Phi())));
  SetBranch("Training_Tree", prefix + "_M",      float(jet.M()));
  SetBranch("Training_Tree", prefix + "_Bscore", bscore);
  SetBranch("Training_Tree", prefix + "_isBtag", int(bscore > btag_wp_cut));
  SetBranch("Training_Tree", prefix + "_isJet", 1);
}

void AtobbMLTree::FillDummyJetBranches(const TString& prefix) {
  SetBranch("Training_Tree", prefix + "_Pt",     -999.f);
  SetBranch("Training_Tree", prefix + "_Eta",    -999.f);
  SetBranch("Training_Tree", prefix + "_CosPhi", -999.f);
  SetBranch("Training_Tree", prefix + "_SinPhi", -999.f);
  SetBranch("Training_Tree", prefix + "_M",      -999.f);
  SetBranch("Training_Tree", prefix + "_Bscore", -999.f);
  SetBranch("Training_Tree", prefix + "_isBtag", 0);
  SetBranch("Training_Tree", prefix + "_isJet", 0);
}
void AtobbMLTree::FillAtoBBHighLevelBranches(
    RVec<Jet>& jets,
    const std::vector<bool>& btag_vector,
    float btag_wp_cut
) {
  const float target_MA = 60.f;

  float A_mass_best = -999.f;
  float A_dr_best = -999.f;
  float A_pt_best = -999.f;
  float A_eta_best = -999.f;
  float A_bscore_sum_best = -999.f;
  float A_bscore_min_best = -999.f;
  int A_pair_idx1 = -1;
  int A_pair_idx2 = -1;

  float bb_mass_min_all = 999999.f;
  float bb_mass_max_all = -999.f;
  float bb_mass_at_min_dr = -999.f;
  float bb_dr_min_all = 999999.f;
  float bb_dr_max_all = -999.f;
  float bb_mass_closest_MA60_all = 999999.f;
  float bb_mass_spread_all = -999.f;

  float bscore_max = -999.f;
  float bscore_2nd = -999.f;
  float bscore_3rd = -999.f;
  float bscore_sum = 0.f;
  float bscore_sum_top4 = 0.f;

  // A candidate pT / mass ratio, bb pair centrality, and mass*dr, which are expected to be somewhat correlated with the correct pairing and the A kinematics
  float A_pt_over_mass = -999.f;
  float A_eta_abs = -999.f;
  float A_mass_dr = -999.f;
  float A_pt_balance_best = -999.f;

  // if (A_mass_best > 0 && A_pt_best > 0 && A_dr_best > 0) {
  //  A_pt_over_mass = A_pt_best / A_mass_best;
  //  A_eta_abs = std::abs(A_eta_best);
  //  A_mass_dr = A_mass_best * A_dr_best;
  // }
  //

  std::vector<std::pair<float, unsigned int>> bscore_pairs;
  std::vector<unsigned int> bjet_indices;

  for (unsigned int i = 0; i < jets.size(); i++) {
    float bscore = jets.at(i).GetTaggerResult(
        JetTagging::JetFlavTagger::ParT,
        JetTagging::JetFlavTaggerScoreType::B
    );

    bscore_sum += bscore;
    bscore_pairs.push_back({bscore, i});

    if (bscore > btag_wp_cut) {
      bjet_indices.push_back(i);
    }
  }

  std::sort(
      bscore_pairs.begin(),
      bscore_pairs.end(),
      [](const auto& a, const auto& b) { return a.first > b.first; }
  );

  if (bscore_pairs.size() > 0) bscore_max = bscore_pairs.at(0).first;
  if (bscore_pairs.size() > 1) bscore_2nd = bscore_pairs.at(1).first;
  if (bscore_pairs.size() > 2) bscore_3rd = bscore_pairs.at(2).first;

  for (unsigned int i = 0; i < std::min<unsigned int>(4, bscore_pairs.size()); i++) {
    bscore_sum_top4 += bscore_pairs.at(i).first;
  }

  float best_mass_diff = 999999.f;

  for (unsigned int a = 0; a < bjet_indices.size(); a++) {
    for (unsigned int b = a + 1; b < bjet_indices.size(); b++) {
      unsigned int i = bjet_indices.at(a);
      unsigned int j = bjet_indices.at(b);

      TLorentzVector b1 = static_cast<TLorentzVector>(jets.at(i));
      TLorentzVector b2 = static_cast<TLorentzVector>(jets.at(j));
      TLorentzVector A = b1 + b2;

      float mass = A.M();
      float dr = b1.DeltaR(b2);

      float bscore1 = jets.at(i).GetTaggerResult(
          JetTagging::JetFlavTagger::ParT,
          JetTagging::JetFlavTaggerScoreType::B
      );
      float bscore2 = jets.at(j).GetTaggerResult(
          JetTagging::JetFlavTagger::ParT,
          JetTagging::JetFlavTaggerScoreType::B
      );

      if (mass < bb_mass_min_all) bb_mass_min_all = mass;
      if (mass > bb_mass_max_all) bb_mass_max_all = mass;

      if (dr < bb_dr_min_all) {
        bb_dr_min_all = dr;
        bb_mass_at_min_dr = mass;
      }

      if (dr > bb_dr_max_all) bb_dr_max_all = dr;

      float mass_diff = std::abs(mass - target_MA);
      if (mass_diff < bb_mass_closest_MA60_all) {
        bb_mass_closest_MA60_all = mass_diff;
      }

      if (mass_diff < best_mass_diff) {
        best_mass_diff = mass_diff;

        float pt1 = b1.Pt();
        float pt2 = b2.Pt();

        if ((pt1 + pt2) > 0) {
          A_pt_balance_best = std::abs(pt1 - pt2) / (pt1 + pt2);
        }

        A_mass_best = mass;
        A_dr_best = dr;
        A_pt_best = A.Pt();
        A_eta_best = A.Eta();
        A_bscore_sum_best = bscore1 + bscore2;
        A_bscore_min_best = std::min(bscore1, bscore2);
        A_pair_idx1 = i;
        A_pair_idx2 = j;
      }
    }
  }

  if (bb_mass_min_all < 999998.f && bb_mass_max_all > -998.f) {
    bb_mass_spread_all = bb_mass_max_all - bb_mass_min_all;
  } else {
    bb_mass_min_all = -999.f;
    bb_mass_max_all = -999.f;
    bb_mass_closest_MA60_all = -999.f;
    bb_dr_min_all = -999.f;
    bb_dr_max_all = -999.f;
  }

  if (A_mass_best > 0 && A_pt_best > 0 && A_dr_best > 0) {
  A_pt_over_mass = A_pt_best / A_mass_best;
  A_eta_abs = std::abs(A_eta_best);
  A_mass_dr = A_mass_best * A_dr_best;
}


  SetBranch("Training_Tree", "A_mass_best", A_mass_best);
  SetBranch("Training_Tree", "A_dr_best", A_dr_best);
  SetBranch("Training_Tree", "A_pt_best", A_pt_best);
  SetBranch("Training_Tree", "A_eta_best", A_eta_best);
  SetBranch("Training_Tree", "A_bscore_sum_best", A_bscore_sum_best);
  SetBranch("Training_Tree", "A_bscore_min_best", A_bscore_min_best);
  SetBranch("Training_Tree", "A_pair_idx1", A_pair_idx1);
  SetBranch("Training_Tree", "A_pair_idx2", A_pair_idx2);

  SetBranch("Training_Tree", "bb_mass_min_all", bb_mass_min_all);
  SetBranch("Training_Tree", "bb_mass_max_all", bb_mass_max_all);
  SetBranch("Training_Tree", "bb_mass_at_min_dr", bb_mass_at_min_dr);
  SetBranch("Training_Tree", "bb_dr_min_all", bb_dr_min_all);
  SetBranch("Training_Tree", "bb_dr_max_all", bb_dr_max_all);
  SetBranch("Training_Tree", "bb_mass_closest_MA60_all", bb_mass_closest_MA60_all);
  SetBranch("Training_Tree", "bb_mass_spread_all", bb_mass_spread_all);

  SetBranch("Training_Tree", "bscore_max", bscore_max);
  SetBranch("Training_Tree", "bscore_2nd", bscore_2nd);
  SetBranch("Training_Tree", "bscore_3rd", bscore_3rd);
  SetBranch("Training_Tree", "bscore_sum", bscore_sum);
  SetBranch("Training_Tree", "bscore_sum_top4", bscore_sum_top4);

  SetBranch("Training_Tree", "A_pt_over_mass", A_pt_over_mass);
  SetBranch("Training_Tree", "A_eta_abs", A_eta_abs);
  SetBranch("Training_Tree", "A_mass_dr", A_mass_dr);

  SetBranch("Training_Tree", "A_pt_balance_best", A_pt_balance_best);
}

void AtobbMLTree::FillTreeBranches(
    Lepton& lepton,
    // RVec<Muon>& muons,
    RVec<Jet>& jets,
    Particle& METv,
    const std::vector<bool>& btag_vector,
    int n_bjets,
    int n_jets_for_count,
    float btag_wp_cut,
    float weight_train,
    float weight_nominal,
    float weight_mc,
    float weight_pu,
    float weight_btag,
    float weight_muon_id,
    float weight_muon_iso,
    float weight_muon_trig,
    float best_chi2,
    float had_W_mass,
    float had_top_mass,
    float lep_W_mass,
    float lep_top_mass
) {
  int label = 0;
  if (MCSample.Contains("AToBB") || MCSample.Contains("HcToWA")) {
    label = 1;
  }

  SetBranch("Training_Tree", "label", label);
  SetBranch("Training_Tree", "isData", int(IsDATA));

  SetBranch("Training_Tree", "weight_train", weight_train);
  SetBranch("Training_Tree", "weight_nominal", weight_nominal);
  SetBranch("Training_Tree", "weight_mc", weight_mc);
  SetBranch("Training_Tree", "weight_pu", weight_pu);
  SetBranch("Training_Tree", "weight_btag", weight_btag);
  SetBranch("Training_Tree", "weight_muon_id", weight_muon_id);
  SetBranch("Training_Tree", "weight_muon_iso", weight_muon_iso);
  SetBranch("Training_Tree", "weight_muon_trig", weight_muon_trig);

  SetBranch("Training_Tree", "nJets", int(jets.size()));
  SetBranch("Training_Tree", "nJets_for_count", n_jets_for_count);
  SetBranch("Training_Tree", "nBJets", n_bjets);

  SetBranch("Training_Tree", "Lepton0_Pt",     float(lepton.Pt()));
  SetBranch("Training_Tree", "Lepton0_Eta",    float(lepton.Eta()));
  SetBranch("Training_Tree", "Lepton0_CosPhi", float(std::cos(lepton.Phi())));
  SetBranch("Training_Tree", "Lepton0_SinPhi", float(std::sin(lepton.Phi())));
  SetBranch("Training_Tree", "Lepton0_M",      float(lepton.M()));
  SetBranch("Training_Tree", "Lepton0_isJet",  0);

  SetBranch("Training_Tree", "MET_Pt",  float(METv.Pt()));
  SetBranch("Training_Tree", "MET_Phi", float(METv.Phi()));

  float HT = 0.f;
  for (const auto& jet : jets) {
    HT += float(jet.Pt());
  }
  SetBranch("Training_Tree", "HT", HT);
  SetBranch("Training_Tree", "ST", HT + float(lepton.Pt()) + float(METv.Pt()));

  for (unsigned int i = 0; i < 6; i++) {
    TString prefix = "Jet" + TString(std::to_string(i));
    if (i < jets.size()) FillJetBranches(prefix, jets.at(i), btag_wp_cut);
    else FillDummyJetBranches(prefix);
  }

  std::vector<unsigned int> bjet_indices;
  for (unsigned int i = 0; i < jets.size(); i++) {
    float bscore = jets.at(i).GetTaggerResult(
        JetTagging::JetFlavTagger::ParT,
        JetTagging::JetFlavTaggerScoreType::B
    );
    if (bscore > btag_wp_cut) bjet_indices.push_back(i);
  }

  float bb_mass_01 = -999.f;
  float bb_mass_02 = -999.f;
  float bb_mass_12 = -999.f;
  float bb_dr_01 = -999.f;
  float bb_dr_02 = -999.f;
  float bb_dr_12 = -999.f;

  if (bjet_indices.size() >= 2) {
    TLorentzVector b0 = static_cast<TLorentzVector>(jets.at(bjet_indices.at(0)));
    TLorentzVector b1 = static_cast<TLorentzVector>(jets.at(bjet_indices.at(1)));
    bb_mass_01 = float((b0 + b1).M());
    bb_dr_01 = float(b0.DeltaR(b1));
  }

  if (bjet_indices.size() >= 3) {
    TLorentzVector b0 = static_cast<TLorentzVector>(jets.at(bjet_indices.at(0)));
    TLorentzVector b1 = static_cast<TLorentzVector>(jets.at(bjet_indices.at(1)));
    TLorentzVector b2 = static_cast<TLorentzVector>(jets.at(bjet_indices.at(2)));

    bb_mass_02 = float((b0 + b2).M());
    bb_mass_12 = float((b1 + b2).M());

    bb_dr_02 = float(b0.DeltaR(b2));
    bb_dr_12 = float(b1.DeltaR(b2));
  }

  SetBranch("Training_Tree", "bb_mass_01", bb_mass_01);
  SetBranch("Training_Tree", "bb_mass_02", bb_mass_02);
  SetBranch("Training_Tree", "bb_mass_12", bb_mass_12);
  SetBranch("Training_Tree", "bb_dr_01", bb_dr_01);
  SetBranch("Training_Tree", "bb_dr_02", bb_dr_02);
  SetBranch("Training_Tree", "bb_dr_12", bb_dr_12);
  
  FillAtoBBHighLevelBranches(jets, btag_vector, btag_wp_cut);

  SetBranch("Training_Tree", "best_chi2", best_chi2);
  SetBranch("Training_Tree", "had_W_mass", had_W_mass);
  SetBranch("Training_Tree", "had_top_mass", had_top_mass);
  SetBranch("Training_Tree", "lep_W_mass", lep_W_mass);
  SetBranch("Training_Tree", "lep_top_mass", lep_top_mass);

  FillTrees("Training_Tree");
}

void AtobbMLTree::EvalChi2(ttCombinatoric& tt_combinatoric) {
  tt_combinatoric.EvalHadronicPart();
  tt_combinatoric.EvalLeptonicPart();

  tt_combinatoric.best_chi2 = 1e9;

  for (unsigned int i = 0; i < tt_combinatoric.neu_pz.size(); i++) {
    double chi2 = this->Chi2Function(
        tt_combinatoric.had_top_mass,
        tt_combinatoric.had_W_mass,
        tt_combinatoric.lep_top_mass.at(i),
        tt_combinatoric.lep_W_mass.at(i)
    );

    tt_combinatoric.chi2.push_back(chi2);

    if (chi2 < tt_combinatoric.best_chi2) {
      tt_combinatoric.best_lep_top_mass = tt_combinatoric.lep_top_mass.at(i);
      tt_combinatoric.best_lep_W_mass   = tt_combinatoric.lep_W_mass.at(i);
      tt_combinatoric.best_neu_pz       = tt_combinatoric.neu_pz.at(i);
      tt_combinatoric.best_chi2         = chi2;
    }
  }
}

double AtobbMLTree::Chi2Function(
    double had_top_mass,
    double had_W_mass,
    double lep_top_mass,
    double lep_W_mass
) {
  double chi2 = 0.;
  chi2 += TMath::Power((had_top_mass - const_top_mass) / const_top_width, 2);
  chi2 += TMath::Power((had_W_mass   - const_w_mass)   / const_w_width,   2);
  chi2 += TMath::Power((lep_top_mass - const_top_mass) / const_top_width, 2);
  chi2 += TMath::Power((lep_W_mass   - const_w_mass)   / const_w_width,   2);
  return chi2;
}

void AtobbMLTree::ttCombinatoric::EvalHadronicPart() {
  TLorentzVector had_W_vector =
      static_cast<TLorentzVector>(jets->at(had_W_jet_idx_1))
    + static_cast<TLorentzVector>(jets->at(had_W_jet_idx_2));

  TLorentzVector had_top_vector =
      had_W_vector
    + static_cast<TLorentzVector>(jets->at(had_top_b_jet_idx));

  had_W_mass = had_W_vector.M();
  had_top_mass = had_top_vector.M();
}

void AtobbMLTree::ttCombinatoric::EvalLeptonicPart() {
  lep_top_mass.clear();
  lep_W_mass.clear();
  neu_pz.clear();
  chi2.clear();

  double step = 5.;
  double pz = -700.;

  while (pz <= 700.) {
    neu_pz.push_back(pz);
    pz += step;
  }

  for (auto& a_neu_pz : neu_pz) {
    double E_neu = TMath::Sqrt(met->E() * met->E() + a_neu_pz * a_neu_pz);
    TLorentzVector neutrino_vector(met->Px(), met->Py(), a_neu_pz, E_neu);

    TLorentzVector lep_W_vector =
        neutrino_vector
      + static_cast<TLorentzVector>(*lepton);

    TLorentzVector lep_top_vector =
        lep_W_vector
      + static_cast<TLorentzVector>(jets->at(lep_top_b_jet_idx));

    lep_W_mass.push_back(lep_W_vector.M());
    lep_top_mass.push_back(lep_top_vector.M());
  }
}
