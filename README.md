# Muon analysis
Plotter and configurations for preliminary muon HLT performance analysis for Phase 2 of the CMS experiment.
All root files, plots and timing information can be found on [this website](https://lferragi.web.cern.ch/) (requires CERN login). Unless stated otherwise

## Thesis plots
Contains plots and comparisons obtained from a ZMM RelVal based on CMSSW_14_1_0_pre5 (results updated august 2024)

## Other contents
- results in [14_0_X](/14_0_X/) are obtained from CMSSW 14_0_0 with D98 geometry and Phase-2 conditions
- results in [14_1_0_pre3](/14_1_0_pre3/) are obtained from 14_1_0_pre3 with D110 geometry with Phase-2 conditions.
- results in [14_1_0_pre5](/14_1_0_pre5/) are obtained from CMSSW_14_1_0_pre5 with D110 geometry and in Phase-2 conditions

Most up to date results: 14_1_0_pre5
They include:
- [seed_validator_comparison](/14_1_0_pre5/seed_validator_comparison/)
- [seeds_and_id](/14_1_0_pre5/seeds_and_id/)
- [new_validation](/14_1_0_pre5/new_validation/)
- [summary_plots](/14_1_0_pre5/summary_plots/)

### Deprecated 14_1_0_pre3

#### [Single Mu no PU sample](/14_1_0_pre3/SingleMu_phase2/)
SingleMu RelVal for matching testing 
- [dR_matching](/14_1_0_pre3/SingleMu_phase2/dR_matching/): matching logic associates stubs with segments that are within a dR window of 0.1 in both DTs and CSCs
- [dPhi+dTheta_matching](/14_1_0_pre3/SingleMu_phase2/dPhi+nHits_matching/): matching logic associates stubs with segments checking in a dPhi window first (0.1), if multiple segments are found in the same window the number of hits in phi (DT) or the total number of hits (CSC) is checked, if there still are multiple matches a window in dTheta is opened (???), finally the hit multiplicity in Theta is checked (DT-only). Logic implemented for barrel and endcap, overlap still relies on dR matching
- [dR_matching](/14_1_0_pre3/SingleMu_phase2/full_new_matching/): extension of new matching logic to the overlap region

#### [ZMM current](/14_1_0_pre3/ZMM_current/)
Current performance on ZMM events with 200PU in phase2 conditions

#### [ZMM Phase2](/14_1_0_pre3/ZMM_phase2/)
Various performance measurements on ZMM events with 200PU in phase2 conditions

#### [PU segments complexity](/14_1_0_pre3/PU_segments_complexity/)
Number of DT/CSC segments per event in ZMM events with 200PU

### Deprecated 14_0_X

#### Single Mu 200PU sample
9k events of Single Mu with 200PU from eos

- [purity_cut_0.75](/14_0_X/purity_cut_0.75/): Preliminary analysis 
- [purity_cut_0.25](/14_0_X/purity_cut_0.25/): Check effects of reduced purity cut in the associators for muons and muon tracks
- [no_purity_cut](/14_0_X/no_purity_cut/): Check effects of reduced purity cut in the associators for muons and muon tracks
- [purity_vs_quality](/14_0_X/purity_vs_quality/): Results of purity cut analysis: no major changes when relaxing the purity cuts

#### ZMM
18k events of ZMM with 200PU from eos

- [ZMM](/14_0_X/ZMM_new_associators/): Features of reconstructed muons 
