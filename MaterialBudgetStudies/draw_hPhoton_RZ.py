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
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kPink
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

def RZ_line_cut(Rconv, eta, Z0):
    return (Rconv+Z0)/TMath.Tan(2*TMath.ATan(TMath.Exp(-1*eta)))

#________________________________________________
def draw_material_RZ(filename_data, filename_mc, suffix, cut, description,structures, folder, date):
    r_bins = [0, 14, 30, 42, 58, 69, 90]#, 180]
    arr_rxy = r_bins
    rootfile_data = TFile.Open(filename_data, "READ");
    rootfile_mc   = TFile.Open(filename_mc  , "READ");

    rootdir_mc_gen  = rootfile_mc.Get("pcm-qc-mc")
    list_gen        = rootdir_mc_gen.Get("Generated");
    list_ev_mc_gen     = rootdir_mc_gen.Get("Event");
    # list_ev_mc_gen  = list_ev_gen.FindObject("PCMPCM");
    rootdir_mc_rec  = rootfile_mc.Get("pcm-qc-mc");
    list_v0_mc_rec  = rootdir_mc_rec.Get("V0");
    list_ev_mc_rec     = rootdir_mc_rec.Get("Event")
    # list_ev_mc_rec  = list_ev_rec.FindObject("PCMPCM");

    list_cut_mc_rec = list_v0_mc_rec.FindObject("qc");

    h1nch_mc_gen    = list_ev_mc_gen.FindObject("hMultNTracksPV").Clone("h1mult");
    nev_gen         = h1nch_mc_gen.GetEntries();
    nch_gen         = h1nch_mc_gen.GetMean();

    h1nch_mc_rec    = list_ev_mc_rec.FindObject("hMultNTracksPV");
    nch_rec         = h1nch_mc_rec.GetMean();
    nev_rec         = h1nch_mc_rec.GetEntries();


    h2rz = list_gen.FindObject("hPhotonRZ").Clone("h2rz");

    if cut >= 100:
        h2rz.GetYaxis().SetRangeUser(0,100);
    else:
        h2rz.GetYaxis().SetRangeUser(0,cut);
    h2rz.GetXaxis().SetRangeUser(-1.*cut,cut)
    h2rz.GetZaxis().SetTitle("#frac{1}{<#it{N}_{ch}^{PV}>} #frac{1}{#it{N}_{ev}} #it{N}_{#gamma} ")
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
        frame1 = p1.DrawFrame(-1.*cut, 0, cut, 90);
    else:
        frame1 = p1.DrawFrame(-1.*cut, 0, cut, cut);
    # frame1 = p1.DrawFrame( 0,-250, 100, 250);
    frame1.GetYaxis().SetTitle("#it{R}_{xy} (cm)");
    frame1.GetXaxis().SetTitle("#it{z} (cm)");
    FrameSettings2D(frame1)

    gPad.SetLogz()
    # rotated_hist.Draw("COLZ SAME")
    h2rz.Draw("COLZ SAME")
    txt = TPaveText(0.0,0.95,1.0,0.92,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(22);#centered,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.04);
    txt.AddText("R_{xy} vs. z, M.C. gen. #gamma (LHC23d1k)")
    txt.Draw();
    ROOT.SetOwnership(txt,False);
    ALICEtext2D("simulation")



    if structures == True:
        R1 =1e-1* (26.7+22.7)/2
        line1 = TLine(RZ_line_cut(R1, 0.9, 7)+3,R1, RZ_line_cut(-R1, 0.9, -7)-3, R1);
        line1.SetLineColor(kRed);
        line1.SetLineStyle(1);
        line1.SetLineWidth(2);
        line1.Draw("");
        ROOT.SetOwnership(line1,False);
    
        R1 =1e-1* (30.1+34.6)/2
        line1 = TLine(RZ_line_cut(R1, 0.9, 7)+3,R1, RZ_line_cut(-R1, 0.9, -7)-3, R1);
        line1.SetLineColor(kRed);
        line1.SetLineStyle(1);
        line1.SetLineWidth(2);
        line1.Draw("");
        ROOT.SetOwnership(line1,False);

        R1 =1e-1* (37.8+42.1)/2
        line1 = TLine(RZ_line_cut(R1, 0.9, 7)+3,R1, RZ_line_cut(-R1, 0.9, -7)-3, R1);
        line1.SetLineColor(kRed);
        line1.SetLineStyle(1);
        line1.SetLineWidth(2);
        line1.Draw("");
        ROOT.SetOwnership(line1,False);
    
        R1 =1e-1* (194.4 + 197.7)/2
        line1 = TLine(RZ_line_cut(R1, 0.9, 7)+3,R1, RZ_line_cut(-R1, 0.9, -7)-3, R1);
        line1.SetLineColor(kRed);
        line1.SetLineStyle(1);
        line1.SetLineWidth(2);
        line1.Draw("");
        ROOT.SetOwnership(line1,False);
    
        R1 =1e-1* (243.0 + 247.0)/2
        line1 = TLine(RZ_line_cut(R1, 0.9, 7)+3,R1, RZ_line_cut(-R1, 0.9, -7)-3, R1);
        line1.SetLineColor(kRed);
        line1.SetLineStyle(1);
        line1.SetLineWidth(2);
        line1.Draw("");
        ROOT.SetOwnership(line1,False);
    
        R1 =1e-1* (342.3 + 345.4)/2
        line1 = TLine(RZ_line_cut(R1, 0.9, 7)+3,R1, RZ_line_cut(-R1, 0.9, -7)-3, R1);
        line1.SetLineColor(kRed);
        line1.SetLineStyle(1);
        line1.SetLineWidth(2);
        line1.Draw("");
        ROOT.SetOwnership(line1,False);
    
        R1 =1e-1* (391.8 + 394.9)/2
        line3 = TLine(RZ_line_cut(R1, 0.9, 7)+3,R1, RZ_line_cut(-R1, 0.9, -7)-3, R1);
        line3.SetLineColor(kRed);
        line3.SetLineStyle(1);
        line3.SetLineWidth(2);
        line3.Draw("");
        ROOT.SetOwnership(line3,False);
    
        R1 =1e-1* (449+461)/2
        line1 = TLine(RZ_line_cut(R1, 0.9, 7)+3,R1, RZ_line_cut(-R1, 0.9, -7)-3, R1);
        line1.SetLineColor(kBlack);
        line1.SetLineStyle(1);
        line1.SetLineWidth(3);
        line1.Draw("");
        ROOT.SetOwnership(line1,False);
    
        R1 =1e-1* (496 + 508)/2
        line1 = TLine(RZ_line_cut(R1, 0.9, 7)+3,R1, RZ_line_cut(-R1, 0.9, -7)-3, R1);
        line1.SetLineColor(kBlack);
        line1.SetLineStyle(1);
        line1.SetLineWidth(3);
        line1.Draw("");
        ROOT.SetOwnership(line1,False);
    
        R1 =1e-1* (540 + 550)/2
        line2 = TLine(RZ_line_cut(R1, 0.9, 7)+3,R1, RZ_line_cut(-R1, 0.9, -7)-3, R1);
        line2.SetLineColor(kBlack);
        line2.SetLineStyle(1);
        line2.SetLineWidth(3);
        line2.Draw("");
        ROOT.SetOwnership(line2,False);

        R1 =60.6
        line4 = TLine(RZ_line_cut(R1, 0.9, 7)+3,R1, RZ_line_cut(-R1, 0.9, -7)-3, R1);
        line4.SetLineColor(kMagenta);
        line4.SetLineStyle(1);
        line4.SetLineWidth(3);
        line4.Draw("");
        ROOT.SetOwnership(line4,False);

        R1 =78.8
        line5 = TLine(RZ_line_cut(R1, 0.9, 7)+3,R1, RZ_line_cut(-R1, 0.9, -7)-3, R1);
        line5.SetLineColor(kPink+5);
        line5.SetLineStyle(1);
        line5.SetLineWidth(3);
        line5.Draw("");
        ROOT.SetOwnership(line5,False);

        leg = TLegend(0.15,0.15,0.3,0.25);
        leg.SetBorderSize(0);
        leg.SetFillColor(kWhite);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.02);
        leg.AddEntry(line3   ,"Layers of the ITS", "L");
        #leg.AddEntry(h1ratio1 ,"M.C. gen / M.C. rec.","LP");
        leg.AddEntry(line2  , "support structure ITS & MFT", "L")
        leg.AddEntry(line4, "TPC inner containment vessel", "L")
        leg.AddEntry(line5, "TPC inner field cage vessel", "L")
        #if generated == True:
         #   leg.AddEntry(h1ratio2 ,"Data / M.C. gen.","LP");
        leg.Draw("");
        ROOT.SetOwnership(leg,False);


    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs(os.path.join(folder,"{0}_material_budget_RZ_cut_{1}_{2}_{3}.png".format(date, cut, description, suffix)));

# if __name__ == "__main__":
#     cutname = "qc"
#     period_mc = "LHC23d1k";
#     period_data = "LHC22f"
#     suffix = "AnyTrack";
#     filename_mc = "/Users/alicamarieenderich/AnalysisResults_LHC23d1k_125889.root"
#     filename_data = "/Users/alicamarieenderich/AnalysisResults_LHC22f4_new_125184.root"
#     date = "this_thesis" #datetime.date.today().strftime("%Y%m%d");
#     folder = "/Users/alicamarieenderich/{0}_material_budget_plots/".format(date);  
#     os.makedirs(folder, exist_ok=True);

#     draw_material_RZ(filename_data, filename_mc, suffix, 40, "wire", folder, date);
#     draw_material_RZ(filename_data, filename_mc, suffix, 250, "complete_comparison", folder, date);