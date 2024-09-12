from ROOT import *
import math

class color:
   PURPLE = '\033[1;35;48m'
   CYAN = '\033[1;36;48m'
   BOLD = '\033[1;37;48m'
   BLUE = '\033[1;34;48m'
   GREEN = '\033[1;32;48m'
   YELLOW = '\033[1;33;48m'
   RED = '\033[1;31;48m'
   BLACK = '\033[1;30;48m'
   UNDERLINE = '\033[4;37;48m'
   END = '\033[1;37;0m'


gStyle.SetOptTitle(0)
gStyle.SetPaintTextFormat("1.3f")
gStyle.SetHistMinimumZero()



fileCurrent = TFile.Open('current/data/current_plots_thesis_ZMM_15k.root')
fileIOFirst = TFile.Open('IO_first/data/IO_first_plots_thesis_ZMM_15k.root')
fileOIFirst =TFile.Open('OI_first/data/OI_first_plots_thesis_ZMM_15k.root')


if not (fileCurrent and fileIOFirst and fileOIFirst) :
    print('[ERROR] File not found')
    sys.exit()

inputFiles = [fileCurrent, fileIOFirst, fileOIFirst]

l2MuFolder = 'DQMData/Run 1/HLT/Run summary/Muon/MuonTrack/hltL2MuonsFromL1TkMuon_UpdAtVtx'

l3IOFolders = ['DQMData/Run 1/HLT/Run summary/Muon/MuonTrack/hltIter2Phase2L3FromL1TkMuonMerged', 
                'DQMData/Run 1/HLT/Run summary/Muon/MuonTrack/phase2L3FilteredObjects_L3IOTrksFiltered',
                'DQMData/Run 1/HLT/Run summary/Muon/MuonTrack/hltIter2Phase2L3FromL1TkMuonMerged']

l3OIFolders = ['DQMData/Run 1/HLT/Run summary/Muon/MuonTrack/hltPhase2L3OIMuonTrackSelectionHighPurity',
               'DQMData/Run 1/HLT/Run summary/Muon/MuonTrack/hltPhase2L3OIMuonTrackSelectionHighPurity',
               'DQMData/Run 1/HLT/Run summary/Muon/MuonTrack/phase2L3FilteredObjects_L3OITrksFiltered']

l3MergedFolder = 'DQMData/Run 1/HLT/Run summary/Muon/MuonTrack/hltPhase2L3MuonMerged'

l3IDFolder = 'DQMData/Run 1/HLT/Run summary/Muon/MuonTrack/hltL3MuonIdTrks'


for i in range(0, len(inputFiles)):
    if i == 0:
        print('--------------------------------------------------------------------------------------------------------')
        print(color.RED + "CURRENT"+ color.END)
        print('--------------------------------------------------------------------------------------------------------')
    elif i == 1:
        print('--------------------------------------------------------------------------------------------------------')
        print(color.RED + "IO FIRST" + color.END)
        print('--------------------------------------------------------------------------------------------------------')
    elif i == 2:
        print('--------------------------------------------------------------------------------------------------------')
        print(color.RED + "OI FIRST" + color.END)
        print('--------------------------------------------------------------------------------------------------------')
    # L2 Muons updated at vertex efficiency
    print(color.CYAN + "L2 Muons" + color.END)
    inputDir = inputFiles[i].GetDirectory(l2MuFolder)
    hSimToRecoBarrel = inputDir.Get('num_assoSimToReco_phi_barrel').Clone()
    hSimBarrel = inputDir.Get('num_simul_phi_barrel').Clone()
    print("Barrel L2 Muon efficiency: ",
          hSimToRecoBarrel.Integral() / hSimBarrel.Integral(), 
          " +/- ", 1/math.sqrt(hSimToRecoBarrel.Integral()))
    hSimToRecoOverlap = inputDir.Get('num_assoSimToReco_phi_overlap').Clone()
    hSimOverlap = inputDir.Get('num_simul_phi_overlap').Clone()
    print("Overlap L2 Muon efficiency: ",
          hSimToRecoOverlap.Integral() / hSimOverlap.Integral(), 
          " +/- ", 1/math.sqrt(hSimToRecoOverlap.Integral()))
    hSimToRecoEndcap = inputDir.Get('num_assoSimToReco_phi_endcap').Clone()
    hSimEndcap = inputDir.Get('num_simul_phi_endcap').Clone()
    print("Endcap L2 Muon efficiency: ",
          hSimToRecoEndcap.Integral() / hSimEndcap.Integral(), 
          " +/- ", 1/math.sqrt(hSimToRecoEndcap.Integral()))
    # L3 IO Tracks efficiency
    print(color.CYAN + "L3 IO Tracks" + color.END)
    inputDir = inputFiles[i].GetDirectory(l3IOFolders[i])
    hSimToRecoBarrel = inputDir.Get('num_assoSimToReco_phi_barrel').Clone()
    hSimBarrel = inputDir.Get('num_simul_phi_barrel').Clone()
    print("Barrel L3 IO Tracks efficiency: ",
          hSimToRecoBarrel.Integral() / hSimBarrel.Integral(), 
          " +/- ", 1/math.sqrt(hSimToRecoBarrel.Integral()))
    hSimToRecoOverlap = inputDir.Get('num_assoSimToReco_phi_overlap').Clone()
    hSimOverlap = inputDir.Get('num_simul_phi_overlap').Clone()
    print("Overlap L3 IO Tracks efficiency: ",
          hSimToRecoOverlap.Integral() / hSimOverlap.Integral(), 
          " +/- ", 1/math.sqrt(hSimToRecoOverlap.Integral()))
    hSimToRecoEndcap = inputDir.Get('num_assoSimToReco_phi_endcap').Clone()
    hSimEndcap = inputDir.Get('num_simul_phi_endcap').Clone()
    print("Endcap L3 IO Tracks efficiency: ",
          hSimToRecoEndcap.Integral() / hSimEndcap.Integral(), 
          " +/- ", 1/math.sqrt(hSimToRecoEndcap.Integral()))
    # L3 OI Tracks efficiency
    print(color.CYAN + "L3 OI Tracks" + color.END)
    inputDir = inputFiles[i].GetDirectory(l3OIFolders[i])
    hSimToRecoBarrel = inputDir.Get('num_assoSimToReco_phi_barrel').Clone()
    hSimBarrel = inputDir.Get('num_simul_phi_barrel').Clone()
    print("Barrel L3 OI Tracks efficiency: ",
          hSimToRecoBarrel.Integral() / hSimBarrel.Integral(), 
          " +/- ", 1/math.sqrt(hSimToRecoBarrel.Integral()))
    hSimToRecoOverlap = inputDir.Get('num_assoSimToReco_phi_overlap').Clone()
    hSimOverlap = inputDir.Get('num_simul_phi_overlap').Clone()
    print("Overlap L3 OI Tracks efficiency: ",
          hSimToRecoOverlap.Integral() / hSimOverlap.Integral(), 
          " +/- ", 1/math.sqrt(hSimToRecoOverlap.Integral()))
    hSimToRecoEndcap = inputDir.Get('num_assoSimToReco_phi_endcap').Clone()
    hSimEndcap = inputDir.Get('num_simul_phi_endcap').Clone()
    print("Endcap L3 OI Tracks efficiency: ",
          hSimToRecoEndcap.Integral() / hSimEndcap.Integral(), 
          " +/- ", 1/math.sqrt(hSimToRecoEndcap.Integral()))
    # L3 Tracks merged
    print(color.CYAN + "L3 Tracks merged" + color.END)
    inputDir = inputFiles[i].GetDirectory(l3MergedFolder)
    hSimToRecoBarrel = inputDir.Get('num_assoSimToReco_phi_barrel').Clone()
    hSimBarrel = inputDir.Get('num_simul_phi_barrel').Clone()
    print("Barrel L3 Tracks merged efficiency: ",
          hSimToRecoBarrel.Integral() / hSimBarrel.Integral(), 
          " +/- ", 1/math.sqrt(hSimToRecoBarrel.Integral()))
    hSimToRecoOverlap = inputDir.Get('num_assoSimToReco_phi_overlap').Clone()
    hSimOverlap = inputDir.Get('num_simul_phi_overlap').Clone()
    print("Overlap L3 Tracks merged efficiency: ",
          hSimToRecoOverlap.Integral() / hSimOverlap.Integral(), 
          " +/- ", 1/math.sqrt(hSimToRecoOverlap.Integral()))
    hSimToRecoEndcap = inputDir.Get('num_assoSimToReco_phi_endcap').Clone()
    hSimEndcap = inputDir.Get('num_simul_phi_endcap').Clone()
    print("Endcap L3 Tracks merged efficiency: ",
          hSimToRecoEndcap.Integral() / hSimEndcap.Integral(), 
          " +/- ", 1/math.sqrt(hSimToRecoEndcap.Integral()))
    # Muon ID
    print(color.CYAN + "Muon ID" + color.END)
    inputDir = inputFiles[i].GetDirectory(l3IDFolder)
    hSimToRecoBarrel = inputDir.Get('num_assoSimToReco_phi_barrel').Clone()
    hSimBarrel = inputDir.Get('num_simul_phi_barrel').Clone()
    print("Barrel Muon ID efficiency: ",
          hSimToRecoBarrel.Integral() / hSimBarrel.Integral(), 
          " +/- ", 1/math.sqrt(hSimToRecoBarrel.Integral()))
    hSimToRecoOverlap = inputDir.Get('num_assoSimToReco_phi_overlap').Clone()
    hSimOverlap = inputDir.Get('num_simul_phi_overlap').Clone()
    print("Overlap Muon ID efficiency: ",
          hSimToRecoOverlap.Integral() / hSimOverlap.Integral(), 
          " +/- ", 1/math.sqrt(hSimToRecoOverlap.Integral()))
    hSimToRecoEndcap = inputDir.Get('num_assoSimToReco_phi_endcap').Clone()
    hSimEndcap = inputDir.Get('num_simul_phi_endcap').Clone()
    print("Endcap Muon ID efficiency: ",
          hSimToRecoEndcap.Integral() / hSimEndcap.Integral(), 
          " +/- ", 1/math.sqrt(hSimToRecoEndcap.Integral()))




#inputFolder = 'DQMData/Run 1/HLT/Run summary/Muon/MuonTrack/hltIter2Phase2L3FromL1TkMuonMerged'
#inputDir = fileCurrent.GetDirectory(inputFolder)
#
#for key in inputDir.GetListOfKeys():
#    print(key.GetName())
#
#
#hSimToRecoBarrel = inputDir.Get('num_assoSimToReco_pT_barrel').Clone()
#hSimBarrel = inputDir.Get('num_simul_pT_barrel').Clone()
#
#print("Barrel efficiency: ",
#       hSimToRecoBarrel.Integral() / hSimBarrel.Integral(), 
#       " +/- ", 1/math.sqrt(hSimToRecoBarrel.Integral()))