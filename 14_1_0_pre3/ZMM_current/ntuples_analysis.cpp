#include <filesystem>
#include <iostream>
#include <map>
#include <optional>
#include <string>
#include <string_view>

#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TPaveStats.h"
#include "TMath.h"
#include "TROOT.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"

struct muon_type
{
    muon_type(TTreeReader &reader, const char *name) : name_{name},
                                                       n_{reader, (string_n + name_).c_str()},
                                                       eta_{reader, (name_ + string_eta).c_str()},
                                                       phi_{reader, (name_ + string_phi).c_str()},
                                                       pt_{reader, (name_ + string_pt).c_str()}

    {
        // Initialize optional features for L3 muons
        if (name_.substr(0, 2) == "l3")
        {
            nhits_pixel_.emplace(reader, (name_ + string_n_pixel_hits).c_str());
            nhits_tracker_.emplace(reader, (name_ + string_n_tracker_hits).c_str());
        }

        if (debug)
        {
            std::cout << name_ << '\n'
                      << (string_n + name_).c_str() << '\n'
                      << (name_ + string_eta).c_str() << '\n'
                      << (name_ + string_phi).c_str() << '\n'
                      << (name_ + string_pt).c_str() << '\n'
                      << (name_ + string_n_pixel_hits).c_str() << '\n'
                      << (name_ + string_n_tracker_hits).c_str() << '\n';
        }
    }

    // Helper strings
    static inline const std::string string_n = "n";
    static inline const std::string string_eta = "_eta";
    static inline const std::string string_phi = "_phi";
    static inline const std::string string_pt = "_pt";
    static inline const std::string string_n_pixel_hits = "_nPixelHits";
    static inline const std::string string_n_tracker_hits = "_nTrkLays";

    // Features
    std::string name_;
    // Common to all types
    TTreeReaderValue<Int_t> n_;
    TTreeReaderArray<Float_t> eta_;
    TTreeReaderArray<Float_t> phi_;
    TTreeReaderArray<Float_t> pt_;
    // Only L3 and global muons have these
    std::optional<TTreeReaderArray<Short_t>> nhits_pixel_;
    std::optional<TTreeReaderArray<Short_t>> nhits_tracker_;

    std::vector<int> good_indexes_;

    bool debug = false;
};

int main()
{

    gROOT->ProcessLine("gErrorIgnoreLevel = 1001;");

    std::string FILE_PATH{"current_ZMM_PU_1k.root"};
    std::string RESULTS_FOLDER{"ntuples_plots/"};

    // Plot booking
    std::map<std::string, TH1 *> histos;
    // Gen
    histos["gen_mu_eta"] = new TH1F("gen_mu_eta", ";Gen #mu #eta;Entries", 100, -2.5, 2.5);
    histos["gen_mu_phi"] = new TH1F("gen_mu_phi", ";Gen #mu #phi;Entries", 100, -TMath::Pi(), TMath::Pi());
    std::vector<std::string> names = {"L1TkMu",
                                      "l3_tk_IO",
                                      "l2_mu_vtx",
                                      "l3_tk_OI",
                                      "l3_tk_merged",
                                      "l3_mu_global"};
    const int nbins = 100;
    const float max_deltaR_plot = 0.1;
    // Reco and Fake
    for (std::size_t i = 0; i != names.size(); i++)
    {
        std::string eta = names[i] + "_eta";
        histos[eta.c_str()] = new TH1F(eta.c_str(), "; #eta; Entries", nbins, -2.5, 2.5);
        if (names[i].substr(0, 2) == "L1")
        {
            std::string eta_reco = names[i] + "_eta_reco";
            histos[eta_reco.c_str()] = new TH1F(eta_reco.c_str(), "; Reco tracks #eta; Entries", nbins, -2.5, 2.5);
            std::string eta_fake = names[i] + "_eta_fake";
            histos[eta_fake.c_str()] = new TH1F(eta_fake.c_str(), "; Fake tracks #eta; Entries", nbins, -2.5, 2.5);
        }
        std::string phi = names[i] + "_phi";
        histos[phi.c_str()] = new TH1F(phi.c_str(), "; #phi; Entries", nbins, -TMath::Pi(), TMath::Pi());
        std::string delta_R = names[i] + "_delta_R";
        histos[delta_R.c_str()] = new TH1F(delta_R.c_str(), "; #Delta R; Entries", nbins, 0, max_deltaR_plot);
        std::string pt_fake = names[i] + "_pt_fake";
        histos[pt_fake.c_str()] = new TH1F(pt_fake.c_str(), "; Fake track pT; Entries", nbins, 0, 50);
        std::string pt_reco = names[i] + "_pt_reco";
        histos[pt_reco.c_str()] = new TH1F(pt_reco.c_str(), "; Reco track pT; Entries", nbins, 0, 100);
        std::string reco_per_event = names[i] + "_reco_per_event";
        histos[reco_per_event.c_str()] = new TH1I(reco_per_event.c_str(), "; Number of Reco objects per event; Occurrences", 6, 0, 6);
        std::string fake_per_event = names[i] + "_fake_per_event";
        histos[fake_per_event.c_str()] = new TH1I(fake_per_event.c_str(), "; Number of Fake objects per event; Occurrences", 18, 0, 18);

        if (names[i].substr(0, 2) == "l3")
        {
            std::string nhits_pixel_fake = names[i] + "_nPixelHits_fake";
            histos[nhits_pixel_fake.c_str()] = new TH1I(nhits_pixel_fake.c_str(), "; Fake Tracks number of hits in the Pixel; Entries", 15, 0, 15);
            std::string nhits_pixel_reco = names[i] + "_nPixelHits_reco";
            histos[nhits_pixel_reco.c_str()] = new TH1I(nhits_pixel_reco.c_str(), "; Reco Tracks number of hits in the Pixel; Entries", 15, 0, 15);
            std::string nhits_tracker_fake = names[i] + "_nTrkLays_fake";
            histos[nhits_tracker_fake.c_str()] = new TH1I(nhits_tracker_fake.c_str(), "; Fake Tracks number of hits in the Tracker; Entries", 15, 0, 15);
            std::string nhits_tracker_reco = names[i] + "_nTrkLays_reco";
            histos[nhits_tracker_reco.c_str()] = new TH1I(nhits_tracker_reco.c_str(), "; Reco Tracks number of hits in the Tracker; Entries", 15, 0, 15);
        }
    }
    // 2d plots
    std::map<std::string, TH2 *> histos_2d;
    histos_2d["n_l2_offline_vs_n_L1TkMu"] = new TH2F("l2_offline_vs_L1TkMu", "Number of L2 Offline Seeds vs Number of L1TkMu; # L1TkMu; # L2 offline seeds", 10, 0, 10, 30, 0, 30);
    histos_2d["n_l2_from_L1TkMu_vs_n_L1TkMu"] = new TH2F("l2_from_L1TkMu_vs_L1TkMu", "Number of L2 Seeds from L1TkMu vs Number of L1TkMu; # L1TkMu; # L2 seeds", 10, 0, 10, 10, 0, 10);
    histos_2d["n_l2_from_L1TkMu_vs_n_l2_offline"] = new TH2F("l2_from_L1TkMu_vs_l2_offline", "Number of L2 Seeds from L1TkMu vs Number of L2 Offline Seeds; # L2 offline seeds; # L2 seeds from L1TkMu", 30, 0, 30, 10, 0, 10);
    histos_2d["n_l2_muons_vs_n_L1TkMu"] = new TH2F("l2_muons_vs_L1TkMu", "Number of L2 Standalone Muons vs Number of L1TkMu; # L1TkMu; # L2 standalone muons", 10, 0, 10, 10, 0, 10);

    // Open the file containing the tree
    TFile *myFile = TFile::Open(FILE_PATH.c_str());
    // Create a TTreeReader for the tree
    TTreeReader reader{"Events", myFile};

    // GEN
    TTreeReaderValue<Int_t> n_gen{reader, "nGenPart"};
    TTreeReaderArray<Float_t> gen_eta{reader, "GenPart_eta"};
    TTreeReaderArray<Float_t> gen_phi{reader, "GenPart_phi"};
    TTreeReaderArray<Int_t> gen_pdgId{reader, "GenPart_pdgId"};
    TTreeReaderArray<Int_t> gen_status{reader, "GenPart_status"};

    // Reco Muons
    std::vector<muon_type> muon_types;
    for (std::size_t i = 0; i != names.size(); ++i)
    {
        muon_types.emplace_back(reader, names[i].c_str());
    }

    // L2 Seeds
    TTreeReaderValue<Int_t> n_l2_seed{reader, "nl2_seed"};
    TTreeReaderValue<Int_t> n_l2_seed_froml1{reader, "nl2_seed_froml1"};
    TTreeReaderValue<Int_t> n_l2_mu{reader, "nl2_mu"};

    int n_events = 1;
    while (reader.Next())
    {
        std::vector<int> gen_index;
        // Fill GEN plots eta, phi
        for (int i_gen = 0; i_gen < *n_gen; ++i_gen)
        {
            if (!(std::abs(gen_pdgId[i_gen]) == 13 && gen_status[i_gen] == 1))
            {
                continue;
            }
            histos["gen_mu_eta"]->Fill(gen_eta[i_gen]);
            histos["gen_mu_phi"]->Fill(gen_phi[i_gen]);
            gen_index.push_back(i_gen);
        }

        // Fill Reco plots eta, phi
        for (std::size_t i = 0; i != names.size(); ++i)
        {
            std::string eta_string = names[i] + "_eta";
            std::string phi_string = names[i] + "_phi";

            for (int i_mu = 0; i_mu < *(muon_types[i].n_); ++i_mu)
            {
                histos[eta_string.c_str()]->Fill(muon_types[i].eta_[i_mu]);
                histos[phi_string.c_str()]->Fill(muon_types[i].phi_[i_mu]);
            }
        }

        // Delta R for all muon types
        // Loop over gen_index to get delta_R
        for (auto i_gen : gen_index)
        {
            for (std::size_t i = 0; i != names.size(); ++i)
            {
                std::string delta_R_string = names[i] + "_delta_R";
                for (int i_mu = 0; i_mu < *(muon_types[i].n_); ++i_mu)
                {
                    float delta_eta = gen_eta[i_gen] - muon_types[i].eta_[i_mu];
                    float delta_phi = std::acos(std::cos(gen_phi[i_gen] - muon_types[i].phi_[i_mu]));
                    float delta_R = std::sqrt(delta_eta * delta_eta + delta_phi * delta_phi);
                    histos[delta_R_string.c_str()]->Fill(delta_R);
                    // Reco Tracks
                    float max_delta_R = names[i].substr(0, 2) == "l2" ? 0.05 : 0.01;
                    if (delta_R < max_delta_R)
                    {
                        // Save good index to filter out later
                        muon_types[i].good_indexes_.push_back(i_mu);
                        // Fill Reco Tracks pt, nPixelHits, and nTrackerHits plots
                        std::string pt_string = names[i] + "_pt_reco";
                        histos[pt_string.c_str()]->Fill(muon_types[i].pt_[i_mu]);
                        if (names[i].substr(0, 2) == "l3")
                        {
                            std::string nhits_pixel_string = names[i] + "_nPixelHits_reco";
                            histos[nhits_pixel_string.c_str()]->Fill(muon_types[i].nhits_pixel_->At(i_mu));
                            std::string nhits_tracker_string = names[i] + "_nTrkLays_reco";
                            histos[nhits_tracker_string.c_str()]->Fill(muon_types[i].nhits_tracker_->At(i_mu));
                        }
                        else if (names[i].substr(0, 2) == "L1")
                        {
                            std::string eta_reco_string = names[i] + "_eta_reco";
                            histos[eta_reco_string.c_str()]->Fill(muon_types[i].eta_[i_mu]);
                        }
                    }
                    // Wait to end the loop before considering fake tracks
                    // One reco could be fake for one of the two GEN Mu from Z
                }
            }
        } // End loop over gen_index

        // Fill Fake Tracks pt, nPixelHits, and nTrackerHits plots
        for (std::size_t i = 0; i != names.size(); ++i)
        {
            for (int i_mu = 0; i_mu != *(muon_types[i].n_); ++i_mu)
            {
                // Avoid good muons
                if (std::find(muon_types[i].good_indexes_.begin(), muon_types[i].good_indexes_.end(), i_mu) != muon_types[i].good_indexes_.end())
                {
                    continue;
                }
                // Fill Fake Tracks eta
                if (names[i].substr(0, 2) == "L1")
                {
                    std::string eta_fake_string = names[i] + "_eta_fake";
                    histos[eta_fake_string.c_str()]->Fill(muon_types[i].eta_[i_mu]);
                }
                // Fill Fake Tracks pt, nPixelHits, and nTrackerHits
                std::string pt_string = names[i] + "_pt_fake";
                histos[pt_string.c_str()]->Fill(muon_types[i].pt_[i_mu]);
                if (names[i].substr(0, 2) == "l3")
                {
                    std::string nhits_pixel_string = names[i] + "_nPixelHits_fake";
                    std::string nhits_tracker_string = names[i] + "_nTrkLays_fake";
                    histos[nhits_pixel_string.c_str()]->Fill(muon_types[i].nhits_pixel_->At(i_mu));
                    histos[nhits_tracker_string.c_str()]->Fill(muon_types[i].nhits_tracker_->At(i_mu));
                }
            }

            // Fill number of Reco and Fake objects per event
            std::string reco_per_event_string = names[i] + "_reco_per_event";
            histos[reco_per_event_string.c_str()]->Fill(muon_types[i].good_indexes_.size());
            std::string fake_per_event_string = names[i] + "_fake_per_event";
            histos[fake_per_event_string.c_str()]->Fill(*(muon_types[i].n_) - muon_types[i].good_indexes_.size());
            // Clear indexes event per event
            muon_types[i].good_indexes_.clear();
        } // End loop over muon types

        // Fill 2D L2 muons/seeds plots
        // #L2 Offline Seeds vs #L1TkMu
        histos_2d["n_l2_offline_vs_n_L1TkMu"]->Fill(*(muon_types[0].n_), *n_l2_seed);
        // #L2 Seeds From L1TkMu vs #L1TkMu
        histos_2d["n_l2_from_L1TkMu_vs_n_L1TkMu"]->Fill(*(muon_types[0].n_), *n_l2_seed_froml1);
        // #L2 Seeds From L1TkMu vs #L2 Offline Seeds
        histos_2d["n_l2_from_L1TkMu_vs_n_l2_offline"]->Fill(*n_l2_seed, *n_l2_seed_froml1);
        // #L2 Muons vs #L1TkMu
        histos_2d["n_l2_muons_vs_n_L1TkMu"]->Fill(*(muon_types[0].n_), *n_l2_mu);

        ++n_events;
    }

    // Save all histograms
    std::filesystem::create_directory(RESULTS_FOLDER);
    auto dt_color = TColor::GetColorTransparent(kOrange - 2, 0.5);
    bool print = true;
    // 1D
    for (const auto &[name, histo] : histos)
    {
        TCanvas c{name.c_str(), name.c_str(), 3000, 1500};
        c.SetGrid();

        histo->SetFillColor(dt_color);
        histo->Draw();
        histo->SetMinimum(0.0);

        if (print)
        {
            if (name.find("_eta") != std::string::npos)
            {
                std::cout << "Total number of " << name << ": " << histo->GetEntries() << '\n';
            }

            if (name.find("_pt_reco") != std::string::npos)
            {
                std::cout << "Number of Reco " << name << ": " << histo->GetEntries() << '\n';
            }

            if (name.find("_pt_fake") != std::string::npos)
            {
                std::cout << "Number of Fake " << name << ": " << histo->GetEntries() << '\n';
            }
        }
        c.SaveAs(Form("%s/%s.png", RESULTS_FOLDER.c_str(), name.c_str()));
        c.SaveAs(Form("%s/%s.pdf", RESULTS_FOLDER.c_str(), name.c_str()));
    }
    // 2D
    for (const auto &[name, histo] : histos_2d)
    {
        TCanvas c{name.c_str(), name.c_str(), 3000, 1500};
        c.SetGrid();

        histo->Draw("COLZ");
        c.Update();
        TPaveStats *stats = (TPaveStats *)c.GetPrimitive("stats");
        stats->SetX1NDC(0.75);
        stats->SetX2NDC(0.9);
        stats->SetY1NDC(0.65);
        stats->SetY2NDC(0.9);

        c.SaveAs(Form("%s/%s.png", RESULTS_FOLDER.c_str(), name.c_str()));
        c.SaveAs(Form("%s/%s.pdf", RESULTS_FOLDER.c_str(), name.c_str()));
    }

    std::cout << "Total number of events " << n_events << std::endl;
}