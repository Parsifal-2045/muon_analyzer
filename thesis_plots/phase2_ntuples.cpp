#include <filesystem>
#include <iostream>
#include <map>
#include <optional>
#include <string>
#include <string_view>
#include <numeric>

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

    std::vector<std::string> FILE_PATHS{"IO_first/data/IO_first_thesis_ZMM_15k.root", "OI_first/data/OI_first_thesis_ZMM_15k.root"};

    for (auto const &FILE_PATH : FILE_PATHS)
    {
        bool verbose = true;
        bool IO_first = FILE_PATH.substr(0, 2) == "IO";

        // Plot booking
        std::map<std::string, TH1 *> histos;
        // Gen
        histos["gen_mu_eta"] = new TH1F("gen_mu_eta", ";Gen #mu #eta;Occurrences", 100, -2.5, 2.5);
        histos["gen_mu_phi"] = new TH1F("gen_mu_phi", ";Gen #mu #phi;Occurrences", 100, -TMath::Pi(), TMath::Pi());
        histos["gen_mu_pt"] = new TH1F("gen_mu_pt", ";Gen #mu pt;Occurrences", 100, 0, 100);

        // L2 seeds from L1TkMu
        histos["l2_seed_pt"] = new TH1F("l2_seed_pt", ";pT;Occurrences", 100, 0., 100.);
        histos["l2_seed_delta_pt"] = new TH1F("l2_seed_delta_pt", "; #Delta pT;Occurrences", 100, 0, 5);
        histos["l2_seed_nhits"] = new TH1I("l2_seed_nhits", ";# of DT/CSC segments per L2 seed; Occurrences", 10, 0, 10);
        histos["phase2_l2_seeds"] = new TH1I("phase2_l2_seeds", ";# L2 seeds from L1 Tracker Muon (new); Occurrences", 30, 0, 30);

        // Reco and Fake
        std::vector<std::string> names = {"L1TkMu",
                                          "l2_mu_vtx",
                                          "l3_tk_IO",
                                          "l3_tk_OI",
                                          "l3_tk_merged",
                                          "l3_mu_global",
                                          "l3_mu_ID"};

        // Add plots for objects ro reuse in IO/OI first cases
        if (IO_first)
        {
            names.push_back("l2_mu_to_reuse");
            names.push_back("l3_tk_IO_filtered");
        }
        else
        {
            names.push_back("L1TkMu_to_reuse");
            names.push_back("l3_tk_OI_filtered");
        }

        if (verbose)
        {
            std::cout << "-------------------------------------------------------------" << std::endl;
            std::cout << (IO_first ? "Inside-Out reconstruction done first" : "Outside-In reconstruction done first") << std::endl;
            std::cout << "-------------------------------------------------------------" << std::endl;
            std::cout << "Objects in event: ";
            for (auto name : names)
            {
                std::cout << name << ", ";
            }
            std::cout << std::endl;
        }

        const int nbins = 100;
        const float max_deltaR_plot = 0.1;
        for (std::size_t i = 0; i != names.size(); i++)
        {
            std::string eta = names[i] + "_eta";
            histos[eta.c_str()] = new TH1F(eta.c_str(), "; #eta; Occurrences", nbins, -2.5, 2.5);
            if (names[i].substr(0, 2) == "L1")
            {
                std::string eta_reco = names[i] + "_eta_reco";
                histos[eta_reco.c_str()] = new TH1F(eta_reco.c_str(), "; Reco tracks #eta; Occurrences", nbins, -2.5, 2.5);
                std::string eta_fake = names[i] + "_eta_fake";
                histos[eta_fake.c_str()] = new TH1F(eta_fake.c_str(), "; Fake tracks #eta; Occurrences", nbins, -2.5, 2.5);
            }
            std::string phi = names[i] + "_phi";
            histos[phi.c_str()] = new TH1F(phi.c_str(), "; #phi; Occurrences", nbins, -TMath::Pi(), TMath::Pi());
            std::string delta_R = names[i] + "_delta_R";
            histos[delta_R.c_str()] = new TH1F(delta_R.c_str(), "; #Delta R; Occurrences", nbins, 0, max_deltaR_plot);
            std::string delta_pt = names[i] + "_delta_pt";
            histos[delta_pt.c_str()] = new TH1F(delta_pt.c_str(), "; #Delta pT; Occurrences", nbins, 0, 5);
            std::string pt_fake = names[i] + "_pt_fake";
            histos[pt_fake.c_str()] = new TH1F(pt_fake.c_str(), "; Fake track pT; Occurrences", nbins, 0, 50);
            std::string pt_reco = names[i] + "_pt_reco";
            histos[pt_reco.c_str()] = new TH1F(pt_reco.c_str(), "; Reco track pT; Occurrences", nbins, 0, 100);
            std::string reco_per_event = names[i] + "_reco_per_event";
            histos[reco_per_event.c_str()] = new TH1I(reco_per_event.c_str(), "; Number of Reco objects per event; Occurrences", 6, 0, 6);
            std::string fake_per_event = names[i] + "_fake_per_event";
            histos[fake_per_event.c_str()] = new TH1I(fake_per_event.c_str(), "; Number of Fake objects per event; Occurrences", 18, 0, 18);

            if (names[i].substr(0, 2) == "l3")
            {
                std::string nhits_pixel_fake = names[i] + "_nPixelHits_fake";
                histos[nhits_pixel_fake.c_str()] = new TH1I(nhits_pixel_fake.c_str(), "; Fake Tracks number of hits in the Pixel; Occurrences", 15, 0, 15);
                std::string nhits_pixel_reco = names[i] + "_nPixelHits_reco";
                histos[nhits_pixel_reco.c_str()] = new TH1I(nhits_pixel_reco.c_str(), "; Reco Tracks number of hits in the Pixel; Occurrences", 15, 0, 15);
                std::string nhits_tracker_fake = names[i] + "_nTrkLays_fake";
                histos[nhits_tracker_fake.c_str()] = new TH1I(nhits_tracker_fake.c_str(), "; Fake Tracks number of hits in the Tracker; Occurrences", 15, 0, 15);
                std::string nhits_tracker_reco = names[i] + "_nTrkLays_reco";
                histos[nhits_tracker_reco.c_str()] = new TH1I(nhits_tracker_reco.c_str(), "; Reco Tracks number of hits in the Tracker; Occurrences", 15, 0, 15);
            }
        }
        // 2d plots
        std::map<std::string, TH2 *> histos_2d;
        // phase2 seeds and l2 muons
        histos_2d["n_phase2_l2_from_L1TkMu_vs_n_L1TkMu"] = new TH2F("phase2_l2_from_L1TkMu_vs_L1TkMu", "Number of L2 Seeds from L1TkMu vs Number of L1TkMu; # L1TkMu; # Phase 2 L2 seeds", 10, 0, 10, 10, 0, 10);
        histos_2d["n_l2_muons_vs_n_L1TkMu"] = new TH2F("l2_muons_vs_L1TkMu", "Number of L2 Standalone Muons vs Number of L1TkMu; # L1TkMu; # L2 Standalone muons", 10, 0, 10, 10, 0, 10);
        histos_2d["phase2_l2_seed_pt_vs_eta"] = new TH2F("phase2_l2_seed_pt_vs_eta", "Phase 2 L2 seeds pT vs #eta; #eta; Phase 2 L2 seeds pT", 100, -2.5, 2.5, 100, 0, 100);
        histos_2d["n_phase2_l2_per_n_L1TkMu_vs_eta"] = new TH2F("n_phase2_l2_per_n_L1TkMu_vs_eta", "Number of L2 Seeds per L1TkMu vs L1TkMu #eta; L1TkMu #eta; # Phase 2 L2 seeds per L1TkMu", 100, -2.5, 2.5, 10, 0, 10);
        // cosmics
        // histos_2d["n_l2_cosmic_seed_vs_n_L1TkMu"] = new TH2F("n_l2_cosmic_seed_vs_n_L1TkMu", "Number of cosmics seeds (displaced) vs Number of L1TkMu; # L1TkMu; # Displaced seeds", 20, 0, 20, 20, 0, 20);
        // histos_2d["n_l2_cosmic_seed_from_l1_vs_n_L1TkMu"] = new TH2F("n_l2_cosmic_seed_from_l1_vs_n_L1TkMu", "Number of cosmics seeds from L1TMu (displaced) vs Number of L1TkMu; # L1TkMu; # Displaced seeds from L1TkMu", 10, 0, 10, 10, 0, 10);
        // histos_2d["n_l2_cosmic_mu_vs_n_L1TkMu"] = new TH2F("n_l2_cosmic_mu_vs_n_L1TkMu", "Number of L2 cosmic Muons (displaced) vs Number of L1TkMu; # L1TkMu; # L2 displaced muons", 10, 0, 10, 10, 0, 10);

        // Open the file containing the tree
        TFile *myFile = TFile::Open(FILE_PATH.c_str());
        // Create a TTreeReader for the tree
        TTreeReader reader{"Events", myFile};

        // GEN
        TTreeReaderValue<Int_t> n_gen{reader, "nGenPart"};
        TTreeReaderArray<Float_t> gen_eta{reader, "GenPart_eta"};
        TTreeReaderArray<Float_t> gen_phi{reader, "GenPart_phi"};
        TTreeReaderArray<Float_t> gen_pt{reader, "GenPart_pt"};
        TTreeReaderArray<Int_t> gen_pdgId{reader, "GenPart_pdgId"};
        TTreeReaderArray<Int_t> gen_status{reader, "GenPart_status"};

        // Reco Muons
        std::vector<muon_type> muon_types;
        for (std::size_t i = 0; i != names.size(); ++i)
        {
            muon_types.emplace_back(reader, names[i].c_str());
        }

        // L2 seeds from L1TkMu
        TTreeReaderValue<Int_t> n_phase2_l2_seed{reader, "nphase2_l2_seed"};
        TTreeReaderValue<Int_t> n_l2_mu{reader, "nl2_mu"};
        TTreeReaderArray<Float_t> l2_seed_pt{reader, "phase2_l2_seed_pt"};
        TTreeReaderArray<Short_t> l2_seed_nhits{reader, "phase2_l2_seed_nHits"}; // # hits for L2 seeds is # of DT/CSC segments matched
        TTreeReaderArray<Float_t> l2_seed_eta{reader, "phase2_l2_seed_eta"};

        // L2 cosmic seeds and muons
        // TTreeReaderValue<Int_t> n_l2_cosmic_seed{reader, "nl2_cosmic_seed"};
        // TTreeReaderValue<Int_t> n_l2_cosmic_seed_froml1{reader, "nl2_cosmic_seed_froml1"};
        // TTreeReaderValue<Int_t> n_l2_cosmic_mu{reader, "nl2_cosmic_mu"};

        int n_events = 0;
        int events_with_no_reuse = 0;
        std::vector<int> n_objects_to_reuse_vector;
        n_objects_to_reuse_vector.reserve(15000);

        while (reader.Next())
        {
            std::vector<int> gen_index;
            // Fill GEN plots eta, phi, pT
            for (int i_gen = 0; i_gen < *n_gen; ++i_gen)
            {
                if (!(std::abs(gen_pdgId[i_gen]) == 13 && gen_status[i_gen] == 1))
                {
                    continue;
                }
                histos["gen_mu_eta"]->Fill(gen_eta[i_gen]);
                histos["gen_mu_phi"]->Fill(gen_phi[i_gen]);
                histos["gen_mu_pt"]->Fill(gen_pt[i_gen]);
                gen_index.push_back(i_gen);
            }

            // Fill L2 seeds plots
            for (int i_seed = 0; i_seed < *n_phase2_l2_seed; ++i_seed)
            {
                histos["l2_seed_pt"]->Fill(l2_seed_pt[i_seed]);
                histos["l2_seed_nhits"]->Fill(l2_seed_nhits[i_seed]);
                histos_2d["phase2_l2_seed_pt_vs_eta"]->Fill(l2_seed_eta[i_seed], l2_seed_pt[i_seed]);
                histos_2d["n_phase2_l2_per_n_L1TkMu_vs_eta"]->Fill(l2_seed_eta[i_seed], *n_phase2_l2_seed / *(muon_types[0].n_));
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
            // Loop over gen_index to get delta_R and delta_pt
            for (auto i_gen : gen_index)
            {
                // delta_R and delta_pt for muons and muon tracks
                for (std::size_t i = 0; i != names.size(); ++i)
                {
                    std::string delta_R_string = names[i] + "_delta_R";
                    std::string delta_pt_string = names[i] + "_delta_pt";
                    for (int i_mu = 0; i_mu < *(muon_types[i].n_); ++i_mu)
                    {
                        float delta_eta = gen_eta[i_gen] - muon_types[i].eta_[i_mu];
                        float delta_phi = std::acos(std::cos(gen_phi[i_gen] - muon_types[i].phi_[i_mu]));
                        float delta_R = std::sqrt(delta_eta * delta_eta + delta_phi * delta_phi);
                        histos[delta_R_string.c_str()]->Fill(delta_R);
                        float delta_pt = std::abs(gen_pt[i_gen] - muon_types[i].pt_[i_mu]);
                        histos[delta_pt_string.c_str()]->Fill(delta_pt);
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
                // delta_pt for phase2 seeds
                for (int i_seed = 0; i_seed != *n_phase2_l2_seed; ++i_seed)
                {
                    float delta_pt = std::abs(gen_pt[i_gen] - l2_seed_pt[i_seed]);
                    histos["l2_seed_delta_pt"]->Fill(delta_pt);
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

            // Fill # of seeds plots
            histos["phase2_l2_seeds"]->Fill(*n_phase2_l2_seed);

            // Fill 2D L2 muons/seeds plots
            // #Phase2 L2 Seeds vs #L1TkMu
            histos_2d["n_phase2_l2_from_L1TkMu_vs_n_L1TkMu"]->Fill(*(muon_types[0].n_), *n_phase2_l2_seed);
            // #L2 Muons vs #L1TkMu
            histos_2d["n_l2_muons_vs_n_L1TkMu"]->Fill(*(muon_types[0].n_), *n_l2_mu);
            // Fill 2D cosmics plots
            // histos_2d["n_l2_cosmic_seed_vs_n_L1TkMu"]->Fill(*(muon_types[0].n_), *n_l2_cosmic_seed);
            // histos_2d["n_l2_cosmic_seed_from_l1_vs_n_L1TkMu"]->Fill(*(muon_types[0].n_), *n_l2_cosmic_seed_froml1);
            // histos_2d["n_l2_cosmic_mu_vs_n_L1TkMu"]->Fill(*(muon_types[0].n_), *n_l2_cosmic_mu);

            if (*(muon_types[7].n_) == 0)
            {
                ++events_with_no_reuse;
            }
            else if (*(muon_types[7].n_) != 0)
            {
                n_objects_to_reuse_vector.push_back(*(muon_types[7].n_));
            }
            ++n_events;
        }

        // Save all histograms
        std::string RESULTS_FOLDER;
        if (IO_first)
        {
            RESULTS_FOLDER = "IO_first/ntuples_analysis/";
        }
        else
        {
            RESULTS_FOLDER = "OI_first/ntuples_analysis/";
        }
        std::filesystem::create_directory(RESULTS_FOLDER);
        auto dt_color = TColor::GetColorTransparent(kOrange - 2, 0.5);

        std::string oFile = RESULTS_FOLDER + "histograms.root";
        TFile *outputFile = new TFile(oFile.c_str(), "RECREATE");

        // 1D
        for (const auto &[name, histo] : histos)
        {
            TCanvas c{name.c_str(), name.c_str(), 3000, 1500};
            c.SetGrid();
            histo->Write();
            histo->SetFillColor(dt_color);
            histo->Draw();
            histo->SetMinimum(0.0);

            if (verbose)
            {
                if (auto found_pos = name.find("_phi"); found_pos != std::string::npos)
                {
                    std::cout << "Total number of " << name.substr(0, found_pos) << ": " << histo->GetEntries() << '\n';
                }

                if (auto found_pos = name.find("_pt_reco"); found_pos != std::string::npos)
                {
                    std::cout << "Number of Reco " << name.substr(0, found_pos) << ": " << histo->GetEntries() << '\n';
                }

                if (auto found_pos = name.find("_pt_fake"); found_pos != std::string::npos)
                {
                    std::cout << "Number of Fake " << name.substr(0, found_pos) << ": " << histo->GetEntries() << '\n';
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
            histo->Write();
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

        std::cout << "Total number of events: " << n_events << std::endl;
        std::cout << "Events with no object to be reused: " << events_with_no_reuse << std::endl;
        std::cout << "Percentage of events with no reuse: " << static_cast<double>(events_with_no_reuse) * 100 / n_events << "%" << std::endl;

        double total_to_reuse = static_cast<double>(std::reduce(n_objects_to_reuse_vector.begin(), n_objects_to_reuse_vector.end()));
        double mean_to_reuse = total_to_reuse / n_objects_to_reuse_vector.size();
        std::cout << "Mean number of objects to reuse in events where necessary: " << mean_to_reuse << " +/- " << mean_to_reuse / std::sqrt(total_to_reuse) << std::endl;

        outputFile->Close();
        delete outputFile;
    }
}