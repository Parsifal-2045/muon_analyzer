#include <iostream>
#include <string>

#include <TFile.h>
#include <TH2.h>

int main()
{
    std::string path = "DQMData/Run 1/HLT/Run summary/Muon/MuonTrack/";
    std::vector<std::string> positions(5);
    positions[0] = "hltIter2Phase2L3FromL1TkMuonMerged";
    positions[1] = "hltL2MuonsFromL1TkMuon_UpdAtVtx";
    positions[2] = "hltL3MuonsPhase2L3OI";
    positions[3] = "hltPhase2L3MuonMerged";
    positions[4] = "hltPhase2L3MuonTrks";

    std::vector<std::string> names(5);
    names[0] = "L3 tracks from L1Tk";
    names[1] = "L2 standalone muons from L1Tk";
    names[2] = "L3 OI global muons";
    names[3] = "L3 tracks no ID";
    names[4] = "L3 tracks Muon ID";

    TFile *input = new TFile("../plots_no_purity_cut.root", "READ");

    std::vector<TH2F *> histos_2d(5);
    std::vector<TH1D *> projections_x(5);
    std::vector<TH1D *> projections_y(5);

    for (int i = 0; i != 5; i++)
    {
        std::string histo = path + positions[i] + "/PurityVsQuality";
        histos_2d[i] = (TH2F *)input->Get(histo.c_str());
        std::string title = "Purity vs Quality " + names[i] + " no purity cut";
        histos_2d[i]->SetTitle(title.c_str());
        std::string px_name = "quality " + names[i];
        projections_x[i] = histos_2d[i]->ProjectionX(px_name.c_str());
        std::string py_name = "purity " + names[i];
        projections_y[i] = histos_2d[i]->ProjectionY(py_name.c_str());
        std::string quality_title = "Quality " + names[i] + " no purity cut";
        projections_x[i]->SetTitle(quality_title.c_str());
        projections_x[i]->GetXaxis()->SetTitle("Quality");
        projections_x[i]->GetYaxis()->SetTitle("Occurrences");
        projections_x[i]->Scale(1./projections_x[i]->Integral());
        std::string purity_title = "Purity " + names[i] + " no purity cut";
        projections_y[i]->SetTitle(purity_title.c_str());
        projections_y[i]->GetXaxis()->SetTitle("Purity");
        projections_y[i]->GetYaxis()->SetTitle("Occurrences");
        projections_y[i]->Scale(1./projections_y[i]->Integral());

    }
    TFile *output = new TFile("no_purity_cut_projections.root", "RECREATE");
    for (int i = 0; i != 5; i++)
    {
        projections_x[i]->Write();
        projections_y[i]->Write();
    }
    input->Close();
    output->Close();
    delete input;
    delete output;
}