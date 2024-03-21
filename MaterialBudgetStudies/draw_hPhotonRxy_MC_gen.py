# This code is based on the code provided by Daiki Sekihata in his repository
# https://github.com/dsekihat/photon_sw_run3 for photon analysis in ALICE Run 3


import re
import numpy as np
import datetime
import math
import ROOT
import ctypes
import os, sys, shutil
import argparse
import yaml
from ROOT import TFile, TDirectory, THashList, TH1F, TH2F, TCanvas, TLegend, TPaletteAxis, TPaveText, TPython, TMath, TF1, TLine, TPython
from ROOT import gStyle, gROOT, gSystem, gPad, TArc
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan
gStyle.SetOptStat(0);
gStyle.SetOptTitle(0);
gStyle.SetPaintTextFormat(".1f")
from FomattingMaterialBudget import FrameSettings, ALICEtext, RatioLegendSettings, FrameSettingsRatio, FrameSettings2D, ALICEtext2D

#________________________________________________

period_data     = "LHC22f";
period_mc       = "LHC23d1k";
filename_mc = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_155278_LHC23d1k.root"
# # filename_mc = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_LHC23d1k_125889.root"
# filename_data = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_148193_LHC22f_apass4.root"
# # filename_data = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_LHC22f4_new_125184.root"

filename_data = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_155756_LHC22f_pass4.root"
cutname = "analysis_wo_mee";

r_bins = [0, 14, 30, 42, 58, 69, 90, 180]
arr_rxy = r_bins

rootfile_data = TFile.Open(filename_data, "READ");
rootfile_mc   = TFile.Open(filename_mc  , "READ");

rootdir_mc_gen  = rootfile_mc.Get("material-budget-mc")
list_gen        = rootdir_mc_gen.Get("Generated");
list_ev_gen     = rootdir_mc_gen.Get("Event")
rootdir_mc_rec  = rootfile_mc.Get("pcm-qc-mc");
list_generated = rootdir_mc_rec.Get("Generated")
list_v0_mc_rec  = rootdir_mc_rec.Get("V0");
list_ev_mc_rec  = rootdir_mc_rec.Get("Event");
list_ev_mc_gen  = list_ev_gen.FindObject("PCMPCM");
list_cut_mc_rec = list_v0_mc_rec.FindObject("qc");
rootdir_mc_pcmqc    = rootfile_mc.Get("pcm-qc-mc");
list_gen_pcmqc        = rootdir_mc_pcmqc.Get("Generated");
list_ev_mc_pcm  = rootdir_mc_pcmqc.Get("Event");
h1nch_mc_gen    = list_ev_mc_pcm.FindObject("hMultNTracksPV").Clone("h1mult");

# h1nch_mc_gen    = list_ev_mc_gen.FindObject("hMultNTracksPV").Clone("h1mult");
nev_gen         = h1nch_mc_gen.GetEntries();
nch_gen         = h1nch_mc_gen.GetMean();

h1nch_mc_rec    = h1nch_mc_gen #list_ev_mc_rec.FindObject("hMultNTracksPV");
nch_rec         = h1nch_mc_rec.GetMean();
nev_rec         = h1nch_mc_rec.GetEntries();

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

def get_bin_edges(axis):
    list_edges = [];
    for i in range(0, axis.GetNbins()+1):
        list_edges.append(axis.GetBinLowEdge(i+1));
    return list_edges;


def plotting_Rxy(histogram, xtitle, ytitle, log , mctype, suffix, circle, outname):
    c1 = TCanvas("c0","c0",0,0,800,800);
    c1.Divide(1,2,1e-5,1e-5); 
    c1.SetMargin(0.13,0.2,0.13,0.13);
    c1.SetTicks(1,1);

    frame1 = c1.DrawFrame(-10, -10, 10, 10);
    FrameSettings2D(frame1)

    if log == "log":
        gPad.SetLogz()
    histogram.Draw("COLZ");

    txt = TPaveText(0.0,0.95,1.0,0.92,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(22);#centered,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.05);
    txt.AddText("M.C. {0}. #gamma (LHC23d1k)".format(mctype))
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    ALICEtext2D("simulation")

    if circle == "circle":

        r = (42+58)/2
        # phi0 = 0.7853981633974482
        # phi1 = 2.399827721492203
        # phi2 = 3.8833575856873837
        # phi3 =  5.541420375081996
        # above were the maxima for the eta vs phi plot, see presentation from August 2023

        # those values are from the 1d plots and are the peaks in the regions between phi around 0/2pi and phi around pi
        phi0 = 1.13
        phi1 = 2.01
        phi2 = 4.28
        phi3 = 5.13

        x0 = r * TMath.Cos(phi0)
        y0 = r * TMath.Sin(phi0)

        x1 = r * TMath.Cos(phi1)
        y1 = r * TMath.Sin(phi1)

        x2 = r * TMath.Cos(phi2)
        y2 = r * TMath.Sin(phi2)

        x3 = r * TMath.Cos(phi3)
        y3 = r * TMath.Sin(phi3)


        for i in range(4):
            x = vars()['x{}'.format(i)]
            y = vars()['y{}'.format(i)]
            radius = (58-42)/2
            circle = TArc(x, y, radius)  # (x, y, radius)
            circle.SetLineColor(ROOT.kRed)
            circle.SetLineWidth(2) 
            circle.SetFillStyle(0);
            circle.Draw("SAME L")
            ROOT.SetOwnership(circle,False)

        circle = TArc(0, 0, 42)  # (x, y, radius)
        circle.SetLineColor(ROOT.kRed)
        circle.SetLineWidth(4) 
        circle.SetFillStyle(0);
        #circle.Draw("SAME L")
        ROOT.SetOwnership(circle,False)

        circle = TArc(0, 0, 58)  # (x, y, radius)
        circle.SetLineColor(ROOT.kRed)
        circle.SetLineWidth(4) 
        circle.SetFillStyle(0);
        #circle.Draw("SAME L")
        ROOT.SetOwnership(circle,False)

        # circle = TArc(0, 0, 42, 0, 90)  # (x, y, radius)
        # circle.SetLineColor(ROOT.kRed)
        # circle.SetLineWidth(4) 
        # circle.SetFillStyle(0);
        # circle.Draw("SAME L")
        # ROOT.SetOwnership(circle,False)

        # circle = TArc(0, 0, 58, 0, 90)  # (x, y, radius)
        # circle.SetLineColor(ROOT.kRed)
        # circle.SetLineWidth(4) 
        # circle.SetFillStyle(0);
        # circle.Draw("SAME L")
        # ROOT.SetOwnership(circle,False)

    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    date = datetime.date.today().strftime("%Y%m%d");

    c1.SaveAs(outname);

    c1.Close();   

def draw_hPhotonRxy_mc_gen(suffix, log, mctype, zoom, circle, folder, date):
    hPhotonRxy = list_generated.FindObject("hPhotonRxy").Clone("hPhotonRxy_gen");
    hPhotonRxy.SetDirectory(0);

    hPhotonRxy.Scale(1/nev_gen);
    hPhotonRxy.Scale(1/nch_gen);
    hPhotonRxy.SetZTitle("#frac{1}{<#it{N}_{ch}^{PV}>} #frac{1}{#it{N}_{ev}} #it{N}_{#gamma}");
    hPhotonRxy.GetZaxis().SetTitleOffset(1.9);

    # # For zooming in structure
    # hPhotonRxy.GetXaxis().SetRangeUser(0, zoom);
    # hPhotonRxy.GetYaxis().SetRangeUser(0, zoom);

    xtitle = "conversion point #it{#varphi} (rad.)"
    ytitle = "#frac{1}{<#it{N}_{ch}^{PV}>} #frac{1}{#it{N}_{ev}} #frac{d^{3}#it{N}_{#gamma}}{d#it{r}_{xy} d#it{#eta} d#it{#varphi}} (cm #upoint rad.)^{#minus1}"
    # date = datetime.date.today().strftime("%Y%m%d");
    outname = os.path.join(folder, "{0}_material_budget_hPhotonRxy_range{1}_{2}_{3}_{4}_{5}.png".format(date, zoom, circle, mctype, log, suffix))
    plotting_Rxy(hPhotonRxy, xtitle , ytitle , log , "gen", suffix, circle, outname)
    
def draw_hPhotonRxy_mc_rec(suffix, log, mctype, zoom, circle, folder, date):
    hPhotonRxy = list_cut_mc_rec.FindObject("hGammaRxy").Clone("hGammaRxy_rec");
    hPhotonRxy.SetDirectory(0);

    hPhotonRxy.Scale(1/nev_rec);
    hPhotonRxy.Scale(1/nch_rec);
    hPhotonRxy.SetZTitle("#frac{1}{<#it{N}_{ch}^{PV}>} #frac{1}{#it{N}_{ev}} #it{N}_{#gamma}");
    hPhotonRxy.GetZaxis().SetTitleOffset(1.9);

    # # For zooming in structure
    # hPhotonRxy.GetXaxis().SetRangeUser(0, zoom);
    # hPhotonRxy.GetYaxis().SetRangeUser(0, zoom);

    xtitle = "conversion point #it{#varphi} (rad.)"
    ytitle = "#frac{1}{<#it{N}_{ch}^{PV}>} #frac{1}{#it{N}_{ev}} #frac{d^{3}#it{N}_{#gamma}}{d#it{r}_{xy} d#it{#eta} d#it{#varphi}} (cm #upoint rad.)^{#minus1}"
    # date = datetime.date.today().strftime("%Y%m%d");
    outname = os.path.join(folder, "{0}_material_budget_hGammaRxy_rec_range{1}_{2}_{3}_{4}_{5}.png".format(date, zoom, circle, mctype, log, suffix))
    plotting_Rxy(hPhotonRxy, xtitle , ytitle , log , "rec", suffix, circle, outname)    
# if __name__ == "__main__":
#     period_data     = "LHC22f";
#     period_mc       = "LHC23d1k";
#     filename_mc     = "/Users/alicamarieenderich/AnalysisResults_LHC23d1k_125889.root"
#     filename_data   = "/Users/alicamarieenderich/AnalysisResults_LHC22f4_new_125184.root"
#     cutname = "analysis_wo_mee";
#     suffix = "" 
#     date = "this_thesis" #datetime.date.today().strftime("%Y%m%d");
#     folder = "/Users/alicamarieenderich/{0}_material_budget_plots/".format(date);  
#     os.makedirs(folder, exist_ok=True);

#     for log in ["log", ""]:
#         for mctype in ["gen"]:
#             for zoom in [100]:
#                 for circle in ["circle", ""]:
#                     draw_hPhotonRxy_mc_gen(suffix, log, mctype , zoom, circle, folder, date);