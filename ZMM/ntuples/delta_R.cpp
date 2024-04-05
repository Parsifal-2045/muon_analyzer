#include <filesystem>
#include <iostream>
#include <map>

#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TROOT.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"

int main()
{
    gROOT->ProcessLine("gErrorIgnoreLevel = 1001;");

    std::string FILE_PATH{"step2_ZMM_full.root"};
    std::string RESULTS_FOLDER{"old_results/"};

    // Plot booking
    std::map<std::string, TH1 *> histos;

    // Gen
    histos["gen_mu_eta"] = new TH1F("gen_mu_eta", ";Gen #mu #eta;Entries", 100, -2.5, 2.5);
    histos["gen_mu_phi"] = new TH1F("gen_mu_phi", ";Gen #mu #phi;Entries", 100, -TMath::Pi(), TMath::Pi());

    // Eta and Phi
    histos["L1Tk_mu_eta"] = new TH1F("L1Tk_mu_eta", ";#mu #eta;Entries", 100, -2.5, 2.5);
    histos["L1Tk_mu_phi"] = new TH1F("L1Tk_mu_phi", ";#mu #phi;Entries", 100, -TMath::Pi(), TMath::Pi());
    histos["L3_mu_IO_eta"] = new TH1F("L3_mu_IO_eta", ";#mu #eta;Entries", 100, -2.5, 2.5);
    histos["L3_mu_IO_phi"] = new TH1F("L3_mu_IO_phi", ";#mu #phi;Entries", 100, -TMath::Pi(), TMath::Pi());
    histos["L2_mu_eta"] = new TH1F("L2_mu_eta", ";#mu #eta;Entries", 100, -2.5, 2.5);
    histos["L2_mu_phi"] = new TH1F("L2_mu_phi", ";#mu #phi;Entries", 100, -TMath::Pi(), TMath::Pi());
    histos["L3_mu_OI_eta"] = new TH1F("L3_mu_OI_eta", ";#mu #eta;Entries", 100, -2.5, 2.5);
    histos["L3_mu_OI_phi"] = new TH1F("L3_mu_OI_phi", ";#mu #phi;Entries", 100, -TMath::Pi(), TMath::Pi());
    histos["L3_mu_no_ID_eta"] = new TH1F("L3_mu_no_ID_eta", ";#mu #eta;Entries", 100, -2.5, 2.5);
    histos["L3_mu_no_ID_phi"] = new TH1F("L3_mu_no_ID_phi", ";#mu #phi;Entries", 100, -TMath::Pi(), TMath::Pi());
    histos["L3_mu_ID_eta"] = new TH1F("L3_mu_ID_eta", ";#mu #eta;Entries", 100, -2.5, 2.5);
    histos["L3_mu_ID_phi"] = new TH1F("L3_mu_ID_phi", ";#mu #phi;Entries", 100, -TMath::Pi(), TMath::Pi());

    // Delta R
    const int nbins = 100;
    const float deltaR_max = 0.5;
    histos["delta_R_L1Tk"] = new TH1F("delta_R_L1Tk", "; #Delta R; Entries", nbins, 0, deltaR_max);
    histos["delta_R_L3_mu_IO"] = new TH1F("delta_R_L3_mu_IO", "; #Delta R; Entries", nbins, 0, deltaR_max);
    histos["delta_R_L2_mu"] = new TH1F("delta_R_L2_mu", "; #Delta R; Entries", nbins, 0, deltaR_max);
    histos["delta_R_L3_mu_OI"] = new TH1F("delta_R_L3_mu_OI", "; #Delta R; Entries", nbins, 0, deltaR_max);
    histos["delta_R_L3_mu_no_ID"] = new TH1F("delta_R_L3_mu_no_ID", "; #Delta R; Entries", nbins, 0, deltaR_max);
    histos["delta_R_L3_mu_ID"] = new TH1F("delta_R_L3_mu_ID", "; #Delta R; Entries", nbins, 0, deltaR_max);

    std::vector<std::string> muon_types = {"L1Tk", "L3_mu_IO", "L2_mu", "L3_mu_OI", "L3_mu_no_ID", "L3_mu_ID"};
    for (int i = 0; i != muon_types.size(); i++)
    {
        std::string name = muon_types[i] + "_pt";
        histos[name.c_str()] = new TH1F(name.c_str(), "; p_T; Entries", 100, 0, 100);
    }

    // Open the file containing the tree
    TFile *myFile = TFile::Open(FILE_PATH.c_str());

    // Create a TTreeReader for the tree, for instance by passing the
    // TTree's name and the TDirectory / TFile it is in.
    TTreeReader reader{"Events", myFile};

    // Using TTreeReader, you just just set to read the branches you need

    // GEN
    TTreeReaderValue<Int_t> n_gen{reader, "nGenPart"};
    TTreeReaderArray<Float_t> gen_eta{reader, "GenPart_eta"};
    TTreeReaderArray<Float_t> gen_phi{reader, "GenPart_phi"};
    TTreeReaderArray<Int_t> gen_pdgId{reader, "GenPart_pdgId"};
    TTreeReaderArray<Int_t> gen_status{reader, "GenPart_status"};

    // L1Tk
    TTreeReaderValue<Int_t> n_L1TkMu{reader, "nL1TkMu"};
    TTreeReaderArray<Float_t> L1TkMu_eta{reader, "L1TkMu_eta"};
    TTreeReaderArray<Float_t> L1TkMu_phi{reader, "L1TkMu_phi"};

    // L3 Muons from L1Tk (IO) RED
    TTreeReaderValue<Int_t> n_L3_IO{reader, "nl3_mu_IO"};
    TTreeReaderArray<Float_t> L3_IO_eta{reader, "l3_mu_IO_eta"};
    TTreeReaderArray<Float_t> L3_IO_phi{reader, "l3_mu_IO_phi"};

    // L2 Muons from L1Tk (updated at vertex) PURPLE
    TTreeReaderValue<Int_t> n_L2{reader, "nl2_mu"};
    TTreeReaderArray<Float_t> L2_mu_eta{reader, "l2_mu_eta"};
    TTreeReaderArray<Float_t> L2_mu_phi{reader, "l2_mu_phi"};

    // L3 muons OI CYAN
    TTreeReaderValue<Int_t> n_L3_OI{reader, "nl3_mu_OI"};
    TTreeReaderArray<Float_t> L3_mu_OI_eta{reader, "l3_mu_OI_eta"};
    TTreeReaderArray<Float_t> L3_mu_OI_phi{reader, "l3_mu_OI_phi"};

    // L3 muons merged (before ID) YELLOW
    TTreeReaderValue<Int_t> n_L3_no_ID{reader, "nl3_mu_merged"};
    TTreeReaderArray<Float_t> L3_mu_no_id_eta{reader, "l3_mu_merged_eta"};
    TTreeReaderArray<Float_t> L3_mu_no_id_phi{reader, "l3_mu_merged_phi"};

    // L3 muons ID GREEN
    TTreeReaderValue<Int_t> n_L3_ID{reader, "nl3_mu_ID"};
    TTreeReaderArray<Float_t> L3_mu_id_eta{reader, "l3_mu_ID_eta"};
    TTreeReaderArray<Float_t> L3_mu_id_phi{reader, "l3_mu_ID_phi"};

    while (reader.Next())
    {
        std::vector<int> index;

        // Fill control gen plots
        for (std::size_t i_gen = 0; i_gen < *n_gen; ++i_gen)
        {
            if (!(std::abs(gen_pdgId[i_gen]) == 13 && gen_status[i_gen] == 1))
            {
                continue;
            }
            histos["gen_mu_eta"]->Fill(gen_eta[i_gen]);
            histos["gen_mu_phi"]->Fill(gen_phi[i_gen]);
            index.push_back(i_gen);
        }

        // Plots and Delta R for validators
        for (auto i_gen : index)
        {
            // L1Tk
            for (std::size_t i_mu = 0; i_mu < *n_L1TkMu; ++i_mu)
            {
                histos["L1Tk_mu_eta"]->Fill(L1TkMu_eta[i_mu]);
                histos["L1Tk_mu_phi"]->Fill(L1TkMu_phi[i_mu]);
                float delta_eta = gen_eta[i_gen] - L1TkMu_eta[i_mu];
                float delta_phi = gen_phi[i_gen] - L1TkMu_phi[i_mu];
                float delta_R = std::sqrt(delta_eta * delta_eta + delta_phi * delta_phi);
                histos["delta_R_L1Tk"]->Fill(delta_R);
            }
            // L3 IO
            for (std::size_t i_mu = 0; i_mu < *n_L3_IO; ++i_mu)
            {
                histos["L3_mu_IO_eta"]->Fill(L3_IO_eta[i_mu]);
                histos["L3_mu_IO_phi"]->Fill(L3_IO_phi[i_mu]);
                float delta_eta = gen_eta[i_gen] - L3_IO_eta[i_mu];
                float delta_phi = gen_phi[i_gen] - L3_IO_phi[i_mu];
                float delta_R = std::sqrt(delta_eta * delta_eta + delta_phi * delta_phi);
                histos["delta_R_L3_mu_IO"]->Fill(delta_R);
            }
            // L2 Muons from L1Tk
            for (std::size_t i_mu = 0; i_mu < *n_L2; ++i_mu)
            {
                histos["L2_mu_eta"]->Fill(L2_mu_eta[i_mu]);
                histos["L2_mu_phi"]->Fill(L2_mu_phi[i_mu]);
                float delta_eta = gen_eta[i_gen] - L2_mu_eta[i_mu];
                float delta_phi = gen_phi[i_gen] - L2_mu_phi[i_mu];
                float delta_R = std::sqrt(delta_eta * delta_eta + delta_phi * delta_phi);
                histos["delta_R_L2_mu"]->Fill(delta_R);
            }
            // L3 Muons OI
            for (std::size_t i_mu = 0; i_mu < *n_L3_OI; ++i_mu)
            {
                histos["L3_mu_OI_eta"]->Fill(L3_mu_OI_eta[i_mu]);
                histos["L3_mu_OI_phi"]->Fill(L3_mu_OI_phi[i_mu]);
                float delta_eta = gen_eta[i_gen] - L3_mu_OI_eta[i_mu];
                float delta_phi = gen_phi[i_gen] - L3_mu_OI_phi[i_mu];
                float delta_R = std::sqrt(delta_eta * delta_eta + delta_phi * delta_phi);
                histos["delta_R_L3_mu_OI"]->Fill(delta_R);
            }
            // L3 Muons no ID
            for (std::size_t i_mu = 0; i_mu < *n_L3_no_ID; ++i_mu)
            {
                histos["L3_mu_no_ID_eta"]->Fill(L3_mu_no_id_eta[i_mu]);
                histos["L3_mu_no_ID_phi"]->Fill(L3_mu_no_id_phi[i_mu]);
                float delta_eta = gen_eta[i_gen] - L3_mu_no_id_eta[i_mu];
                float delta_phi = gen_phi[i_gen] - L3_mu_no_id_phi[i_mu];
                float delta_R = std::sqrt(delta_eta * delta_eta + delta_phi * delta_phi);
                histos["delta_R_L3_mu_no_ID"]->Fill(delta_R);
            }
            // L3 Muons ID
            for (std::size_t i_mu = 0; i_mu < *n_L3_ID; ++i_mu)
            {
                histos["L3_mu_ID_eta"]->Fill(L3_mu_id_eta[i_mu]);
                histos["L3_mu_ID_phi"]->Fill(L3_mu_id_phi[i_mu]);
                float delta_eta = gen_eta[i_gen] - L3_mu_id_eta[i_mu];
                float delta_phi = gen_phi[i_gen] - L3_mu_id_phi[i_mu];
                float delta_R = std::sqrt(delta_eta * delta_eta + delta_phi * delta_phi);
                histos["delta_R_L3_mu_ID"]->Fill(delta_R);
            }
        }
    }

    std::filesystem::create_directory(RESULTS_FOLDER);

    auto dt_color = TColor::GetColorTransparent(kOrange - 2, 0.5);

    for (const auto &[name, histo] : histos)
    {
        TCanvas c{name.c_str(), name.c_str(), 3000, 1500};
        c.SetGrid();

        histo->SetFillColor(dt_color);
        histo->Draw();
        histo->SetMinimum(0.0);

        c.SaveAs(Form("%s/%s.png", RESULTS_FOLDER.c_str(), name.c_str()));
        // c.SaveAs(Form("%s/%s.pdf", RESULTS_FOLDER.c_str(), name.c_str()));
    }
}