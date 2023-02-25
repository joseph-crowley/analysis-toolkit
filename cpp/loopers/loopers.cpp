/**
 * @file loopers.cpp
 * @brief Implements data analysis looper for the CMS Compact Muon Solenoid experiment.
 *
 * @author Joe Crowley
 * @date February 20, 2023 
 *
 * @see analysis.py
 * @see plotting.py
 *
 */

#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"
#include "TChain.h"
#include "TTreeCache.h"
#include "TTreeCacheUnzip.h"
#include "TTreePerfStats.h"
#include "TCanvas.h"
#include "TPad.h"
#include "THStack.h"
#include "TStyle.h"
#include "TText.h"
#include "TLine.h"
#include "TLegend.h"
#include "TRatioPlot.h"
#include "TLatex.h"
#include "TLorentzVector.h"

#include <iostream>
#include <iomanip>

using namespace std;

/**
 * @brief A struct to store the event data.
 * 
 * This struct stores the event data, including the event weight, the number of jets, the transverse momentum of the missing energy, the transverse momentum of the jets, the b-tagging status of the jets, the transverse momentum, pseudorapidity, and azimuthal angle of the leptons, and the minimum b-jet and lepton b-jet mass for each b-tagging working point.
 *  - event_wgt is the product of the luminosity, the cross section, and the number of events.
 *  - event_wgt_triggers_dilepton_matched is the product of the single lepton trigger efficiency and the dilepton trigger efficiency.
 *  - event_wgt_SFs_btagging is the product of the b-tagging scale factor for each b-tagged jet.
 *  - event_wgt_SFs_btagging is calculated using the b-tagging efficiency and the b-tagging mistag rate.
 *  - njet is the number of tight ak4jets.
 *  - PFMET_pt_final is the transverse momentum of the missing energy.
 *  - jet_pt is a vector of the transverse momentum of the jets.
 *  - jet_is_btagged is a vector of the b-tagging status of the jets.
 *  - lep_pt is a vector of the transverse momentum of the leptons.
 *  - lep_eta is a vector of the pseudorapidity of the leptons.
 *  - lep_phi is a vector of the azimuthal angle of the leptons.
 *  - lep_mass is a vector of the mass of the leptons.
 */
struct EventData {
    float event_wgt;
    float event_wgt_triggers_dilepton_matched;
    float event_wgt_TTxsec;
    float event_wgt_SFs_btagging;
    unsigned int njet;
    float PFMET_pt_final;
    std::vector<float> *jet_pt;
    std::vector<unsigned char> *jet_is_btagged;
    std::vector<float> *lep_pt;
    std::vector<float> *lep_eta;
    std::vector<float> *lep_phi;
    std::vector<float> *lep_mass;
};

/**
 * @brief Get the histograms used in the analysis.
 *
 *  Histograms are stored in a map, where the key is the name and the value is the histogram.
 *  - nbjet is the number of b-jets.
 *  - njet is the number of non-btagged jets.
 *  - lep1_pt is the transverse momentum of the leading lepton.
 *  - lep1_eta is the pseudorapidity of the leading lepton.
 *  - lep1_phi is the azimuthal angle of the leading lepton.
 *  - lep2_pt is the transverse momentum of the sub-leading lepton.
 *  - lep2_eta is the pseudorapidity of the sub-leading lepton.
 *  - lep2_phi is the azimuthal angle of the sub-leading lepton.
 *  
 *  Histograms are categorized by the number of b-jets and the b-tagging working point.
 * 
 *  @param btag_categories A vector of strings containing the b-tagging working points.
 *  @param nb_categories A vector of strings containing the number of b-jets.
 * 
 * @return A map containing the histograms used in the analysis.
 *
 */
std::unordered_map<std::string, TH1D *> get_histograms(std::vector<std::string> btag_categories, std::vector<std::string> nb_categories) {
    std::unordered_map<std::string, TH1D *> hists;


    // Initialize histograms for the jets and b-jets multiplicity and kinematics
    hists["jetpt"] = new TH1D("jetpt", "p_{T} of all non-btagged Jets", 20, 0, 20);
    for (const auto& category : btag_categories) {
        hists["bjetpt_" + category] = new TH1D(("bjetpt_" + category).c_str(), ("p_{T} of all " + category + " b-Jets").c_str(), 20, 0, 20);
        hists["nbjet_" + category] = new TH1D(("nbjet_" + category).c_str(), ("Number of " + category + " b-Jets").c_str(), 20, 0, 20);
        hists["njet_" + category] = new TH1D(("njet_" + category).c_str(), "Number of Jets", 20, 0, 20);
    }

    // Initialize histograms with fixed btag WP 
    for (const auto& category : nb_categories) {
        // Initialize histograms for the lepton kinematics
        hists["lep1_pt_" + category] = new TH1D(("lep1_pt_" + category).c_str(), "Leading Lepton p_{T}", 100, 0, 500);
        hists["lep1_eta_" + category] = new TH1D(("lep1_eta_" + category).c_str(), "Leading Lepton #eta", 100, -2.5, 2.5);
        hists["lep1_phi_" + category] = new TH1D(("lep1_phi_" + category).c_str(), "Leading Lepton #phi", 100, -3.14, 3.14);
        hists["lep2_pt_" + category] = new TH1D(("lep2_pt_" + category).c_str(), "Sub-leading Lepton p_{T}", 100, 0, 500);
        hists["lep2_eta_" + category] = new TH1D(("lep2_eta_" + category).c_str(), "Sub-leading Lepton #eta", 100, -2.5, 2.5);
        hists["lep2_phi_" + category] = new TH1D(("lep2_phi_" + category).c_str(), "Sub-leading Lepton #phi", 100, -3.14, 3.14);

        // Initialize histograms for the dilepton kinematics
        hists["m_ll_" + category] = new TH1D(("m_ll_" + category).c_str(), "Dilepton Mass", 100, 0, 1000);
        hists["pt_ll_" + category] = new TH1D(("pt_ll_" + category).c_str(), "Dilepton p_{T}", 100, 0, 500);

        // Initialize histograms for other kinematic quantities
        hists["PFMET_pt_final_" + category] = new TH1D(("PFMET_pt_final_" + category).c_str(), "PFMET p_{T}", 100, 0, 500);
        hists["Ht_" + category] = new TH1D("Ht", "H_{T}", 100, 0, 1000);
    }

    return hists;
}

/**
 * @brief A function to process a single event.
 *
 * This function takes a single event and processes it, filling histograms.
 * The histograms are filled calculated from the event weight, the trigger weight, and the b-tagging scale factor. 
 * The raw event weight is the product of the luminosity, the cross section, and the number of events.
 *     - trigger weight is the product of the single lepton trigger efficiency and the dilepton trigger efficiency.
 *     - b-tagging scale factor is the product of the b-tagging scale factor for each b-tagged jet.
 *     - b-tagging scale factor is calculated using the b-tagging efficiency and the mistag rate.
 * 
 * @param data An EventData object containing the event data. 
 * @param hists A vector of pairs of the histogram name and the histogram vector to be filled. Passed by reference.
 * @param btag_categories A vector of strings containing the b-tagging working points.
 * @param nb_categories A vector of strings containing the number of b-jets.
 * @param btag_WP The b-tagging working point.
 * @param PFMET_pt_final_threshold The threshold for the PFMET_pt_final variable.
 * @param pt_threshold_btagged The threshold for the b-tagged jet pT.
 * @param pt_threshold_unbtagged The threshold for the untagged jet pT.j
 * 
 */
void process_event(EventData data, std::unordered_map<std::string, TH1D *> &hists, std::vector<std::string> btag_categories, std::vector<std::string> nb_categories, unsigned int btag_WP, float PFMET_pt_final_threshold, float pt_threshold_btagged, float pt_threshold_unbtagged) {
    std::cout << "eventing..." << endl;
    double event_total_wgt = data.event_wgt * data.event_wgt_triggers_dilepton_matched * data.event_wgt_SFs_btagging * data.event_wgt_TTxsec;
    bool passMETCut = data.PFMET_pt_final > PFMET_pt_final_threshold;

    // Calculate variables needed for filling histograms
    std::vector<unsigned int> nbjet_ct;
    //std::vector<float> min_mlb;
    //std::vector<float> min_mbb;
    
    for (unsigned int i = 0; i < btag_categories.size(); i++) {
        nbjet_ct.push_back(0);
        //min_mlb.push_back(9999);
        //min_mbb.push_back(9999);
    }

    unsigned int njet_ct = 0;
    float Ht = 0;
    for (unsigned int i = 0; i < data.njet; i++) {
        // is_btagged is an unsigned char 
        // defined by int(is_btagged_loose) + int(is_btagged_medium) + int(is_btagged_tight)
        unsigned char is_btagged = data.jet_is_btagged->at(i);
        float const& pt = data.jet_pt->at(i);

        // check if the jet passes the pt threshold
        // note: pt_threshold has to be redefined after the b-tagging check since btagged jets require a different pt threshold
        float pt_threshold = (is_btagged ? pt_threshold_btagged : pt_threshold_unbtagged);
        if (is_btagged && pt < pt_threshold) is_btagged = 0;
        pt_threshold = (is_btagged ? pt_threshold_btagged : pt_threshold_unbtagged);

        if (pt > pt_threshold) {
            Ht += pt;
            
            // Fill histograms for the jet and b-jet kinematics
            if (!passMETCut) continue;

            // count the number of jets and b-jets
            if (is_btagged == 0) {
                njet_ct++;
                hists["jetpt"]->Fill(pt, event_total_wgt);
                // skip the rest of the loop if the jet is not b-tagged
                continue; 
            }

            // at this point, the jet is b-tagged
            // count the number of b-jets
            if (is_btagged >= 1) nbjet_ct.at(0)++;
            if (is_btagged >= 2) nbjet_ct.at(1)++;
            if (is_btagged >= 3) nbjet_ct.at(2)++;

            // fill the bjet pt for all btag categories less than or equal to is_btagged
            for (unsigned int is_btagged_category = 0; is_btagged_category < is_btagged; is_btagged_category++) {
                hists["bjetpt_" + btag_categories.at(is_btagged_category)]->Fill(pt, event_total_wgt);
            }
        }
    }

    // Calculate the dilepton kinematics
    TLorentzVector lep1;
    TLorentzVector lep2;
    lep1.SetPtEtaPhiM(data.lep_pt->at(0), data.lep_eta->at(0), data.lep_phi->at(0), data.lep_mass->at(0));
    lep2.SetPtEtaPhiM(data.lep_pt->at(1), data.lep_eta->at(1), data.lep_phi->at(1), data.lep_mass->at(1));
    TLorentzVector dilep = lep1 + lep2;

    // Fill histograms
    std::string category;

    for (unsigned int i_btag_category = 0; i_btag_category < btag_categories.size(); i_btag_category++){
        // loose medium tight categories
        category = btag_categories.at(i_btag_category);

        if (passMETCut) hists["nbjet_" + category]->Fill(nbjet_ct.at(i_btag_category), event_total_wgt);
    }

    for (unsigned int i_nb_category = 0; i_nb_category < nb_categories.size(); i_nb_category++){
        // Skip events that don't match the category
        category = nb_categories.at(i_nb_category);
        switch (i_nb_category) {
            case 0: // lt2
                if (nbjet_ct.at(btag_WP) >= 2) continue;
                break;
            case 1: // eq2
                if (nbjet_ct.at(btag_WP) != 2) continue;
                break;
            case 2: // gt2
                if (nbjet_ct.at(btag_WP) <= 2) continue;
                break;
            case 3: // eq1
                if (nbjet_ct.at(btag_WP) != 1) continue;
                break;
            case 4: // eq0
                if (nbjet_ct.at(btag_WP) != 0) continue;
                break;
            default:
                break;
        }

        if (passMETCut) {
            // jet multiplicity hists
            hists["nbjet_" + category ]->Fill(nbjet_ct.at(btag_WP), event_total_wgt);
            hists["njet_" + category]->Fill(njet_ct, event_total_wgt);

            // lepton hists
            hists["lep1_pt_" + category]->Fill(data.lep_pt->at(0), event_total_wgt);
            hists["lep1_eta_" + category]->Fill(data.lep_eta->at(0), event_total_wgt);
            hists["lep1_phi_" + category]->Fill(data.lep_phi->at(0), event_total_wgt);
            hists["lep2_pt_" + category]->Fill(data.lep_pt->at(1), event_total_wgt);
            hists["lep2_eta_" + category]->Fill(data.lep_eta->at(1), event_total_wgt);
            hists["lep2_phi_" + category]->Fill(data.lep_phi->at(1), event_total_wgt);
    
            // dilepton hists
            hists["m_ll_" + category]->Fill(dilep.M(), event_total_wgt);
            hists["pt_ll_" + category]->Fill(dilep.Pt(), event_total_wgt);
            //hists["m_lb_" + category]->Fill(min_mlb, event_total_wgt);
            //hists["m_bb_" + category]->Fill(min_mbb, event_total_wgt);
        }

        // No MET cut
        hists["PFMET_pt_final_" + category]->Fill(data.PFMET_pt_final, event_total_wgt);
        hists["Ht_" + category]->Fill(Ht, event_total_wgt);
    }
}

/**
 * @brief A function to process a TChain of files.
 * 
 * This function takes a TChain of files and processes them, producing histograms. 
 * The histograms are filled in the process_event function.
 * 
 * @param chain A pointer to the TChain of files to be processed.
 * @param sample_str A string containing the sample name.
 * 
 */
void process_chain(TChain *chain, std::string sample_str) {
    std::cout << "hello." << endl;

    // Set up the cuts
    constexpr unsigned int btag_WP = 1; // 0 = loose, 1 = medium, 2 = tight
    constexpr float PFMET_pt_final_threshold = 50.;
    constexpr float pt_threshold_btagged = 40.;
    constexpr float pt_threshold_unbtagged = 25.;

    // initialize the histograms
    std::vector<std::string> btag_categories = {"btagLooseWP", "btagMediumWP", "btagTightWP"};
    std::vector<std::string> nb_categories = {"nb_lt_2", "nb_eq_2", "nb_gt_2", "nb_eq_1", "nb_eq_0"};
    std::unordered_map<std::string, TH1D *> hists = get_histograms(btag_categories, nb_categories);  

    std::cout << "got some hists." << endl;

    // Set up the event data structure
    EventData eventData;
    chain->SetBranchAddress("event_wgt", &eventData.event_wgt);
    chain->SetBranchAddress("event_wgt_triggers_dilepton_matched", &eventData.event_wgt_triggers_dilepton_matched);
    chain->SetBranchAddress("event_wgt_SFs_btagging", &eventData.event_wgt_SFs_btagging);
    chain->SetBranchAddress("nak4jets_tight_pt25", &eventData.njet);
    chain->SetBranchAddress("pTmiss", &eventData.PFMET_pt_final);
    chain->SetBranchAddress("ak4jets_pt", &eventData.jet_pt);
    chain->SetBranchAddress("ak4jets_pass_btagging", &eventData.jet_is_btagged);
    chain->SetBranchAddress("leptons_pt", &eventData.lep_pt);
    chain->SetBranchAddress("leptons_eta", &eventData.lep_eta);
    chain->SetBranchAddress("leptons_phi", &eventData.lep_phi);
    chain->SetBranchAddress("leptons_mass", &eventData.lep_mass);

    std::cout << "set up event struct." << endl;

    // TODO: remove manual setting of weights for the TTbar samples and data
    if (sample_str.find("TT_") != std::string::npos) eventData.event_wgt_TTxsec = 0.826;
    if (sample_str.find("Data") != std::string::npos) eventData.event_wgt_SFs_btagging = 1.;

    // Event loop
    int n_entries = chain->GetEntries();
    for (int i = 0; i < n_entries; i++) {
        std::cout << "loopin..." << endl;
        chain->GetEntry(i);
        std::cout << "got entry..." << endl;
        process_event(eventData, hists, btag_categories, nb_categories, btag_WP, PFMET_pt_final_threshold, pt_threshold_btagged, pt_threshold_unbtagged); 
    }

    // Rebin and write histograms to file 
    std::string outfile_name = "hists_"+ std::string(sample_str) + "_" + std::string(btag_categories.at(btag_WP)) + "bpt" + std::to_string(pt_threshold_btagged) + "_jpt" + std::to_string(pt_threshold_unbtagged) + ".root";
    TFile *f = new TFile(outfile_name.c_str(), "RECREATE");
    TH1D *hist;
    int nbin;
    for (auto &pair : hists) {
        // hists is an unordered map, so take the second element of the pairs, which is the hist for the category
        hist = pair.second;
        nbin = hist->GetNbinsX();
        hist->SetBinContent(nbin, hist->GetBinContent(nbin + 1) + hist->GetBinContent(nbin));
        hist->SetBinError(nbin, std::sqrt(std::pow(hist->GetBinError(nbin + 1), 2) + std::pow(hist->GetBinError(nbin), 2)));
        hist->Write();
    }
    f->Close();
}

void loopers() { return; } 
