# Muon analysis
Plotter and configurations for preliminary muon HLT performance analysis for Phase 2 of the CMS experiment.

## Contents
All root files, plots and timing information can be found on [this website](https://lferragi.web.cern.ch/) (requires CERN login). Unless stated otherwise, results are obtained from CMSSW 14_0_0 with D98 geometry and Phase2 conditions. 

### Single Mu 200PU sample
9k events of Single Mu with 200PU from eos

- [purity_cut_0.75](/purity_cut_0.75/): Preliminary analysis 
- [purity_cut_0.25](/purity_cut_0.25/): Check effects of reduced purity cut in the associators for muons and muon tracks
- [no_purity_cut](/no_purity_cut/): Check effects of reduced purity cut in the associators for muons and muon tracks
- [purity_vs_quality](/purity_vs_quality/): Results of purity cut analysis: no major changes when relaxing the purity cuts

### ZMM
18k events of ZMM with 200PU from eos

- [ZMM](/ZMM/): Features of reconstructed muons 
- [ntuples](/ntuples/): Analysis of Fake and Reco tracks to understand where the contamination comes from
