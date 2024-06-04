import os
import sys
import ROOT
from ROOT import gStyle, gPad
from array import array
import pandas as pd
import numpy as np

gStyle.SetOptTitle(1)
gStyle.SetOptStat(1110)

dt = pd.read_csv("/home/luca/muon_analyzer/14_1_0_pre3/segments_complexity/PU_segments.csv")["dt"]
csc = pd.read_csv("/home/luca/muon_analyzer/14_1_0_pre3/segments_complexity/PU_segments.csv")["csc"]

dt = array('i', dt)
csc = array('i',csc)

dtHisto = ROOT.TH1I("dt","Number of DT segments per ZMM PU 200 event; # DT segments; Occurrences",50, 0, 50)
cscHisto = ROOT.TH1I("csc", "Number of CSC segments per ZMM PU 200 event; # CSC segments; Occurrences", 200, 0, 200)

for value in dt:
    dtHisto.Fill(value)

for value in csc:
    cscHisto.Fill(value)

canvas = ROOT.TCanvas("canvas", "Number of DT/CSC segments per event", 3000, 2000)
canvas.Divide(1,2)
canvas.cd(1)
dtHisto.Draw()
canvas.cd(2)
cscHisto.Draw()

canvas.SaveAs("/home/luca/muon_analyzer/14_1_0_pre3/segments_complexity/PU_segments_complexity.png")
canvas.SaveAs("/home/luca/muon_analyzer/14_1_0_pre3/segments_complexity/PU_segments_complexity.pdf")


# wait for input to keep the GUI (which lives on a ROOT event dispatcher) alive
#if __name__ == '__main__':
#    rep = ''
#    while not rep in ['q', 'Q']:
#        # Check if we are in Python 2 or 3
#        if sys.version_info[0] > 2:
#            rep = input('enter "q" to quit: ')
#        else:
#            rep = raw_input('enter "q" to quit: ')
#        if 1 < len(rep):
#            rep = rep[0]