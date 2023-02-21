/**
 * @file loopers.cpp
 * @brief Implements a fast and optimized data analysis loop for the CMS Compact Muon Solenoid experiment.
 *
 * This file contains the implementation of a C++ data analysis loop for the CMS experiment. The loop uses advanced techniques for performance optimization, such as vectorization and multithreading, to quickly process large amounts of data and produce histograms for further analysis.
 *
 * @author Joe Crowley
 * @date February 20, 2023 
 *
 * @see analysis.py
 * @see plotting.py
 *
 * @copyright Copyright (c) University of California, Santa Barbara, 2023
 */

#include <vector>
#include <thread>
#include <algorithm>
#include <iostream>
#include <unordered_map>

#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"

#include "dependencies/tqdm.h"

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
 *  The histograms are stored in a map of vectors, where the key is the histogram name and the vector contains the histograms for each b-tagging working point.
 *  - nbjet is the number of b-jets.
 *  - njet is the number of non-btagged jets.
 *  - lep1_pt is the transverse momentum of the leading lepton.
 *  - lep1_eta is the pseudorapidity of the leading lepton.
 *  - lep1_phi is the azimuthal angle of the leading lepton.
 *  - lep2_pt is the transverse momentum of the sub-leading lepton.
 *  - lep2_eta is the pseudorapidity of the sub-leading lepton.
 *  - lep2_phi is the azimuthal angle of the sub-leading lepton.
 * 
 * @return A map of vectors containing the histograms used in the analysis.
 *
 */
std::unordered_map<std::string, std::vector<TH1D *>> get_histograms() {
    std::unordered_map<std::string, std::vector<TH1D *>> hists;
    std::vector<std::string> btag_categories = {"btagLooseWP", "btagMediumWP", "btagTightWP"};
    std::vector<std::string> nb_categories = {"nb_lt_2", "nb_eq_2", "nb_gt_2"};

    // Initialize histograms for the number of jets and b-jets
    for (const auto& category : btag_categories) {
        hists["h_nbjet_" + category].push_back(new TH1D(("h_nbjet_" + category).c_str(), "Number of b-Jets", 20, 0, 20));
        hists["h_njet_" + category].push_back(new TH1D(("h_njet_" + category).c_str(), "Number of Jets", 20, 0, 20));
    }

    // Initialize histograms for the lepton kinematics
    for (const auto& category : nb_categories) {
        hists["h_lep1_pt_" + category].push_back(new TH1D(("h_lep1_pt_" + category).c_str(), "Leading Lepton p_{T}", 100, 0, 500));
        hists["h_lep1_eta_" + category].push_back(new TH1D(("h_lep1_eta_" + category).c_str(), "Leading Lepton #eta", 100, -2.5, 2.5));
        hists["h_lep1_phi_" + category].push_back(new TH1D(("h_lep1_phi_" + category).c_str(), "Leading Lepton #phi", 100, -3.14, 3.14));
        hists["h_lep2_pt_" + category].push_back(new TH1D(("h_lep2_pt_" + category).c_str(), "Sub-leading Lepton p_{T}", 100, 0, 500));
        hists["h_lep2_eta_" + category].push_back(new TH1D(("h_lep2_eta_" + category).c_str(), "Sub-leading Lepton #eta", 100, -2.5, 2.5));
        hists["h_lep2_phi_" + category].push_back(new TH1D(("h_lep2_phi_" + category).c_str(), "Sub-leading Lepton #phi", 100, -3.14, 3.14));
    }

    // Initialize histograms for the dilepton kinematics
    for (const auto& category : nb_categories) {
        hists["h_m_ll_" + category].push_back(new TH1D(("h_m_ll_" + category).c_str(), "Dilepton Mass", 100, 0, 1000));
        hists["h_pt_ll_" + category].push_back(new TH1D(("h_pt_ll_" + category).c_str(), "Dilepton p_{T}", 100, 0, 500));
    }

    // Initialize histograms for other kinematic quantities
    for (const auto& category : nb_categories) {
        hists["h_PFMET_pt_final_" + category].push_back(new TH1D(("h_PFMET_pt_final_" + category).c_str(), "PFMET p_{T}", 100, 0, 500));
        hists["h_Ht"].push_back(new TH1D("h_Ht", "H_{T}", 100, 0, 1000));
    }

    return hists;
}

/**
 * @brief A function to process a single event.
 *
 * This function takes a single event and processes it, producing histograms.
 * The histograms are filled with the weights calculated from the event weight, the trigger weight, and the b-tagging scale factor. 
 * The raw event weight is the product of the luminosity, the cross section, and the number of events.
 *     - trigger weight is the product of the single lepton trigger efficiency and the dilepton trigger efficiency.
 *     - b-tagging scale factor is the product of the b-tagging scale factor for each b-tagged jet.
 *     - b-tagging scale factor is calculated using the b-tagging efficiency and the mistag rate.
 * 
 * @param data An EventData object containing the event data. 
 * @param hists A vector of pairs of the histogram name and the histogram vector to be filled. Passed by reference.
 * 
 */
void process_event(EventData data, std::unordered_map<std::string, std::vector<TH1D *>> &hists, std::vector<std::string> btag_categories, std::vector<std::string> nb_categories, unsigned int btag_WP) {
    double event_total_wgt = data.event_wgt * data.event_wgt_triggers_dilepton_matched * data.event_wgt_SFs_btagging * data.event_wgt_TTxsec;

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
        auto is_btagged = data.jet_is_btagged->at(i);
        auto const& pt = data.jet_pt->at(i);

        // TODO: add thresholds and WPs as arguments instead
        constexpr float pt_threshold_btagged = 40.;
        constexpr float pt_threshold_unbtagged = 25.;
        float pt_threshold = (is_btagged ? pt_threshold_btagged : pt_threshold_unbtagged);
        if (is_btagged && pt < pt_threshold) is_btagged = 0;
        pt_threshold = (is_btagged ? pt_threshold_btagged : pt_threshold_unbtagged);

        if (pt > pt_threshold) {
            Ht += pt;
            
            if (is_btagged == 0 && data.PFMET_pt_final > 50.) h_jetpt.front()->Fill(pt, event_total_wgt);
            if (is_btagged > 0 && data.PFMET_pt_final > 50.) h_bjetpt.at(is_btagged - 1)->Fill(pt, event_total_wgt);
            if (is_btagged == 0) njet_ct++;
            if (is_btagged >= 1) nbjet_ct.at(0)++;
            if (is_btagged >= 2) nbjet_ct.at(1)++;
            if (is_btagged >= 3) nbjet_ct.at(2)++;
        }
    }

    TLorentzVector lep1;
    TLorentzVector lep2;
    lep1.SetPtEtaPhiM(data.lep_pt->at(0), data.lep_eta->at(0), data.lep_phi->at(0), data.lep_mass->at(0));
    lep2.SetPtEtaPhiM(data.lep_pt->at(1), data.lep_eta->at(1), data.lep_phi->at(1), data.lep_mass->at(1));
    TLorentzVector dilep = lep1 + lep2;

    // Fill histograms

    for (unsigned int i_btag_category = 0; i_btag_category < 3; i_btag_category++){
        // loose medium tight categories
        if (data.PFMET_pt_final > 50.) h_nbjet.at(i_btag_category)->Fill(nbjet_ct.at(i_btag_category), event_total_wgt);
    }

    for (unsigned int i_nb_category = 0; i_nb_category < 3; i_nb_category++){
        // lt2 eq2 gt2 categories

        if ((i_nb_category == 0) && (nbjet_ct.at(btag_WP) >= 2)) continue;
        if ((i_nb_category == 1) && (nbjet_ct.at(btag_WP) != 2)) continue;
        if ((i_nb_category == 2) && (nbjet_ct.at(btag_WP) <= 2)) continue;

        if (data.PFMET_pt_final > 50.) {
            h_nbjet.at(i_nb_category)->Fill(nbjet_ct.at(btag_WP), event_total_wgt);
            h_njet.at(i_nb_category)->Fill(njet_ct, event_total_wgt);
            h_lep1_pt.at(i_nb_category)->Fill(data.lep_pt->at(0), event_total_wgt);
            h_lep1_eta.at(i_nb_category)->Fill(data.lep_eta->at(0), event_total_wgt);
            h_lep1_phi.at(i_nb_category)->Fill(data.lep_phi->at(0), event_total_wgt);
            h_lep2_pt.at(i_nb_category)->Fill(data.lep_pt->at(1), event_total_wgt);
            h_lep2_eta.at(i_nb_category)->Fill(data.lep_eta->at(1), event_total_wgt);
            h_lep2_phi.at(i_nb_category)->Fill(data.lep_phi->at(1), event_total_wgt);
    
            // dilepton hists
            h_m_ll.at(i_nb_category)->Fill(dilep.M(), event_total_wgt);
            h_pt_ll.at(i_nb_category)->Fill(dilep.Pt(), event_total_wgt);
            //h_m_lb.at(i_nb_category)->Fill(min_mlb, event_total_wgt);
            //h_m_bb.at(i_nb_category)->Fill(min_mbb, event_total_wgt);
        }
        h_met.at(i_nb_category)->Fill(data.PFMET_pt_final, event_total_wgt);
        h_Ht.at(i_nb_category)->Fill(Ht, event_total_wgt);
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
    std::unordered_map<std::string, std::vector<TH1D *>> hists = get_histograms();  

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

    // TODO: remove manual setting of weights for the TTbar samples and data
    if (sample_str.find("TT_") != std::string::npos) eventData.event_wgt_TTxsec = 0.826;
    if (sample_str.find("Data") != std::string::npos) eventData.event_wgt_SFs_btagging = 1.;

    int n_entries = chain->GetEntries();
    for (int i = 0; i < n_entries; i++) {
        chain->GetEntry(i);
        process_event(eventData, hists, btag_categories, nb_categories, btag_WP); 
    }

    // Rebin and write histograms to file 
    TFile *f = new TFile("hists_"+ sample_str + ".root", "RECREATE");
    // hists is an unordered map, so take the second element of the pairs, which is the vector hists for each category
    for (auto &pair : hists) {
        for (auto &hist : pair.second) {
            int nbin = hist->GetNbinsX();
            hist->SetBinContent(nbin, hist->GetBinContent(nbin + 1) + hist->GetBinContent(nbin));
            hist->SetBinError(nbin, std::sqrt(std::pow(hist->GetBinError(nbin + 1), 2) + std::pow(hist->GetBinError(nbin), 2)));
            hist->Write();
        }
    }
    f->Close();

    delete event;
}
