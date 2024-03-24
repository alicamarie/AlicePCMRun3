# This code is based on the code provided by Daiki Sekihata in his repository
# https://github.com/dsekihat/photon_sw_run3 for photon analysis in ALICE Run 3


import re, os
import numpy as np
import datetime
import math
import ROOT
import ctypes
from ROOT import TFile, TDirectory, THashList, TH1F, TH2F, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TLine, TPython
from ROOT import gStyle, gROOT, gSystem, gPad
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan
gStyle.SetOptStat(0);
gStyle.SetOptTitle(0);
from FomattingMaterialBudget import FrameSettings, ALICEtext, FrameSettingsRatio, RatioLegendSettings, FrameSettings2D, ALICEtext2D

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
def draw_material_RZ(filename_data, filename_mc, suffix, cut, description, folder, date):
    r_bins = [0, 14, 30, 42, 58, 69, 90]#, 180]
    arr_rxy = r_bins
    rootfile_data = TFile.Open(filename_data, "READ");
    rootfile_mc   = TFile.Open(filename_mc  , "READ");

    rootdir_mc_gen  = rootfile_mc.Get("material-budget-mc")
    list_gen        = rootdir_mc_gen.Get("Generated");
    list_ev_gen     = rootdir_mc_gen.Get("Event");
    list_ev_mc_gen  = list_ev_gen.FindObject("PCMPCM");
    rootdir_mc_rec  = rootfile_mc.Get("material-budget-mc");
    list_v0_mc_rec  = rootdir_mc_rec.Get("V0");
    list_ev_rec     = rootdir_mc_rec.Get("Event")
    list_ev_mc_rec  = list_ev_rec.FindObject("PCMPCM");
    list_cut_mc_rec = list_v0_mc_rec.FindObject("qc");

    h1nch_mc_gen    = list_ev_mc_gen.FindObject("hMultNTracksPV").Clone("h1mult");
    nev_gen         = h1nch_mc_gen.GetEntries();
    nch_gen         = h1nch_mc_gen.GetMean();

    h1nch_mc_rec    = list_ev_mc_rec.FindObject("hMultNTracksPV");
    nch_rec         = h1nch_mc_rec.GetMean();
    nev_rec         = h1nch_mc_rec.GetEntries();

    h2rz = list_gen.FindObject("hPhotonRZ").Clone("h2rz");

    if cut >= 100:
        h2rz.GetYaxis().SetRangeUser(0,90);
    else:
        h2rz.GetYaxis().SetRangeUser(0,cut);
    h2rz.GetXaxis().SetRangeUser(-1.*cut,cut)
    h2rz.GetZaxis().SetTitleOffset(1.9);
#style   
    make_common_style(h2rz, 20, 1.0, kBlue+1, 1, 0);
    ROOT.SetOwnership(h2rz, False);

#canvas plotting
    c1 = TCanvas("c1","c1",0,0,800,800);
    c1.Divide(1,2,1e-5,1e-5); 
    c1.SetTicks(1,1);
    gPad.SetLogz()
    p1 = c1
    p1.SetMargin(0.13,0.2,0.13,0.13);

    if cut >= 100:
        frame1 = p1.DrawFrame(-1.*cut, 0, cut, 100);
    else:
        frame1 = p1.DrawFrame(-1.*cut, 0, cut, cut);
    # frame1 = p1.DrawFrame( 0,-250, 100, 250);
    frame1.GetYaxis().SetTitle("#it{R}_{xy} (cm)");
    frame1.GetXaxis().SetTitle("#{z} (cm)");
    FrameSettings2D(frame1)

    gPad.SetLogz()
    h2rz.Draw("COLZ SAME")
    txt = TPaveText(0.0,0.95,1.0,0.92,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(22);#centered,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.04);
    txt.AddText("R_{xy} vs. Z, M.C. gen. #gamma (LHC23d1k)")
    txt.Draw();
    ROOT.SetOwnership(txt,False);
    ALICEtext2D("thesis")
    
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs(os.path.join(folder,"{0}_material_budget_RZ_cut_{1}_{2}_{3}.png".format(date, cut, description, suffix)));