# This code is based on the code provided by Daiki Sekihata in his repository
# https://github.com/dsekihat/photon_sw_run3 for photon analysis in ALICE Run 3


import re, os
import numpy as np
import datetime
import math
import ROOT
import ctypes
from ROOT import TFile, TDirectory, THashList, TH1F, TH2F, TCanvas, TLegend, TPaletteAxis, TPaveText, TPython, TMath, TF1, TLine, TPython, TArc
from ROOT import gStyle, gROOT, gSystem, gPad
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kViolet
gStyle.SetOptStat(0);
gStyle.SetOptTitle(0);
from FomattingMaterialBudget import FrameSettings, ALICEtext, RatioLegendSettings, FrameSettingsRatio, FrameSettings2D, ALICEtext2D

#________________________________________________
def make_common_style(g1,marker,size,color,width=1,fill=0):
    g1.SetMarkerStyle(marker);
    g1.SetMarkerColor(color);
    g1.SetMarkerSize(size);
    g1.SetLineColor(color);
    g1.SetLineWidth(width);
    g1.SetFillColor(color);
    g1.SetFillStyle(fill);

#________________________________________________
def draw_AP_with_cuts(filename_data, filename_mc, cutname, suffix, folder, date):
    rootfile_data = TFile.Open(filename_data, "READ");
    rootdir_data = rootfile_data.Get("pcm-qc");
    list_v0_data = rootdir_data.Get("V0");
    list_ev_data = rootdir_data.Get("Event");
    list_cut_data = list_v0_data.FindObject(cutname);

    h1ev_data = list_ev_data.FindObject("hCollisionCounter").Clone("h1ev");
    nev_data = h1ev_data.GetBinContent(4);

    h1nch_data = list_ev_data.FindObject("hMultNTracksPV");
    nch_data = h1nch_data.GetMean();

    h2_data = list_cut_data.FindObject("hAPplot");
    h2_data.Sumw2();
    h2_data.SetXTitle("#it{#alpha} = (p_{L}^{-}-p_{L}^{+})/(p_{L}^{-}+p_{L}^{+})")
    h2_data.GetZaxis().SetTitleOffset(1.9);
    h2_data.SetDirectory(0);
    h2_data.GetXaxis().SetRangeUser(-100,100);
    h2_data.GetYaxis().SetRangeUser(-100,100);
    ROOT.SetOwnership(h2_data, False);

    c1 = TCanvas("c0","c0",0,0,800,800);
    c1.Divide(1,2,1e-5,1e-5); 
    c1.SetMargin(0.13,0.13,0.13,0.13);
    c1.SetTicks(1,1);

    frame1 = c1.DrawFrame(-1, 0, 1, 0.25);
    FrameSettings2D(frame1)

    gPad.SetLogz()
    h2_data.Draw("COLZ");

    ALICEtext2D("thesis")

    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    outname = os.path.join(folder, "{0}_material_budget_AP_plot_with_cuts_{1}_{2}.pdf".format(date, cutname, suffix))
    c1.SaveAs(outname);