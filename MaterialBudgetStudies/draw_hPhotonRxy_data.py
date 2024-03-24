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
def draw_hPhotonRxy_data(filename_data, filename_mc, cutname, suffix, folder, date):
    rootfile_data = TFile.Open(filename_data, "READ");
    rootdir_data = rootfile_data.Get("pcm-qc");
    list_v0_data = rootdir_data.Get("V0");
    list_ev_data = rootdir_data.Get("Event");
    list_cut_data = list_v0_data.FindObject(cutname);

    h1ev_data = list_ev_data.FindObject("hCollisionCounter").Clone("h1ev");
    nev_data = h1ev_data.GetBinContent(4);

    h1nch_data = list_ev_data.FindObject("hMultNTracksPV");
    nch_data = h1nch_data.GetMean();

    h2_data = list_cut_data.FindObject("hGammaRxy")#("hGammaRxy_recalc");
    h2_data.Sumw2();
    h2_data.SetZTitle("#frac{1}{<N_{ch}^{PV}>} #frac{1}{N_{ev}} N_{#gamma} ");
    h2_data.GetZaxis().SetTitleOffset(1.9);
    h2_data.SetDirectory(0);
    h2_data.GetXaxis().SetRangeUser(-100,100); # previous: -100
    h2_data.GetYaxis().SetRangeUser(-100,100); # previous: -100
    h2_data.GetXaxis().SetTitleOffset(1.9)
    ROOT.SetOwnership(h2_data, False);

    h2_data.Scale(1/nev_data);
    h2_data.Scale(1/nch_data);

    c1 = TCanvas("c0","c0",0,0,800,800);
    c1.Divide(1,2,1e-5,1e-5); 
    c1.SetMargin(0.13,0.2,0.13,0.13);
    c1.SetTicks(1,1);

    frame1 = c1.DrawFrame(-100, -100, 100, 100);
    frame1.GetXaxis().SetTitle("conversion point #it{#varphi} (rad.)");
    frame1.GetYaxis().SetTitle("#frac{1}{<#it{N}_{ch}^{PV}>} #frac{1}{#it{N}_{ev}} #it{N}_{#gamma} ");
    FrameSettings2D(frame1)

    gPad.SetLogz()
    h2_data.Draw("COLZ");

    txt = TPaveText(0.0,0.95,1.0,0.92,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(22);#centered,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.05);
    txt.AddText("Data #gamma (LHC22f)")
    txt.Draw();
    ROOT.SetOwnership(txt,False);
    ALICEtext2D("thesis")

    r = (42+58)/2
    phi0 = 0.7853981633974482
    phi1 = 2.399827721492203
    phi2 = 3.8833575856873837
    phi3 =  5.541420375081996
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

    # circle = TArc(0, 0, 42)  # (x, y, radius)
    # circle.SetLineColor(ROOT.kRed)
    # circle.SetLineWidth(4) 
    # circle.SetFillStyle(0);
    # circle.Draw("SAME L")
    # ROOT.SetOwnership(circle,False)

    # circle = TArc(0, 0, 58)  # (x, y, radius)
    # circle.SetLineColor(ROOT.kRed)
    # circle.SetLineWidth(4) 
    # circle.SetFillStyle(0);
    # circle.Draw("SAME L")
    # ROOT.SetOwnership(circle,False)

    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    outname = os.path.join(folder, "{0}_material_budget_hPhotonRxy_data_circle_{1}_{2}.png".format(date, cutname, suffix))
    c1.SaveAs(outname);