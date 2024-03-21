# This code is based on the code provided by Daiki Sekihata in his repository
# https://github.com/dsekihat/photon_sw_run3 for photon analysis in ALICE Run 3


import re, os
import numpy as np
import datetime
import math
import ROOT
import ctypes
from ROOT import TBox, TFile, TDirectory, THashList, TH1F, TH2F, TH2Poly, TCanvas, TLegend, TPave, TPaletteAxis, TPaveText, TPaveLabel, TPython, TMath, TF1, TLine, TPython
from ROOT import gStyle, gROOT, gSystem, gRandom, gPad
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan
gStyle.SetOptStat(0);
gStyle.SetOptTitle(0);
from FomattingMaterialBudget import FrameSettings, ALICEtext, FrameSettingsRatio, RatioLegendSettings, FrameSettings2Dtwoplots, ALICEtext2Dtwoplots

#________________________________________________
def make_common_style(g1,marker,size,color,width=1,fill=0):
    g1.SetMarkerStyle(marker);
    g1.SetMarkerColor(color);
    g1.SetMarkerSize(size);
    g1.SetLineColor(color);
    g1.SetLineWidth(width);
    g1.SetFillColor(color);
    g1.SetFillStyle(fill);

def find_phi_for_eta(hist, eta_value, output_file):
    eta_bin = hist.GetYaxis().FindBin(eta_value)
    x_bins = hist.GetNbinsX()

    with open(output_file, 'w') as file:
        for x_bin in range(1, x_bins + 1):
            phi = hist.GetXaxis().GetBinCenter(x_bin)
            nev = hist.GetBinContent(x_bin, eta_bin)

            if nev > 4e-7:
                output_line = f"phi: {phi}, nev: {nev}\n"
                file.write(output_line)

#________________________________________________

def draw_material_z_vs_phi(folder, filename_data, filename_mc, cutname, period, rid, suffix, log, date):
    rootfile_data = TFile.Open(filename_data, "READ");
    rootfile_mc   = TFile.Open(filename_mc  , "READ");

    list_data = rootfile_data.Get(cutname);
    list_mc = rootfile_mc.Get(cutname);

    r_bins = [0, 14, 30, 42, 58, 69, 90, 180]
    
    h2data_complete = list_data.FindObject("h2etaphi_r{0:d}".format(rid));
    h2mc_complete = list_mc.FindObject("h2etaphi_r{0:d}".format(rid));
    h2data_complete.SetDirectory(0);
    h2mc_complete.SetDirectory(0);

#style   
    make_common_style(h2data_complete, 20, 1.0, kBlue+1, 1, 0);
    make_common_style(h2mc_complete  , 20, 1.0, kRed+1, 1, 0);
    ROOT.SetOwnership(h2data_complete, False);
    ROOT.SetOwnership(h2mc_complete, False);

#canvas plotting
    c1 = TCanvas("c1","c1",0,0,800,800);
    c1.Divide(1,2,1e-5,1e-5); 

    p1 = c1.cd(1);
    p1.SetMargin(0.1,0.13,0.13,0.2);
    p1.SetTicks(1,1);

    frame1 = p1.DrawFrame(0, -2, TMath.TwoPi(), 2);
    frame1.GetXaxis().SetTitle("#it{#varphi} (rad.)");
    frame1.GetYaxis().SetTitle("pseudorapidity #it{#eta}");
    FrameSettings2Dtwoplots(frame1)

    if log == "log":
        gPad.SetLogz()

    h2data_complete.Draw("COLZ SAME")
    max_z_value_data = h2data_complete.GetMaximum()

    txt = TPaveText(0.0,0.9,1.0,0.92,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(22);#centered,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.05);
    txt.AddText("#it{{#eta}} vs. #it{{#varphi}} in {0:3.2f} cm < #it{{r}}_{{xy}} < {1:3.2f} cm".format(r_bins[rid], r_bins[rid+1]))
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    txt = TPaveText(0.1,0.5,0.4,1.15,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.05);
    txt.AddText("Data #gamma candidates (LHC22f)")
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    ALICEtext2Dtwoplots("thesis",1)

    p2 = c1.cd(2);
    p2.SetMargin(0.1,0.13,0.22,0.05);
    p2.SetTicks(1,1);

    frame2 = p2.DrawFrame(0, -2, TMath.TwoPi(), 2);
    frame2.GetXaxis().SetTitle("#it{#varphi} (rad.)");
    frame2.GetYaxis().SetTitle("pseudorapidity #it{#eta}");
    FrameSettings2Dtwoplots(frame2)
    ROOT.SetOwnership(frame2,False);

    if log == "log":
        gPad.SetLogz()

# changed here to see difference
        #!!!!!!!!
    # h2mc_complete.SetMaximum(max_z_value_data)
    h2mc_complete.Draw("COLZ SAME")

    txt = TPaveText(0.1,0.8,0.4,1.15,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.05);
    txt.AddText("M.C. rec. #gamma (LHC23d1k)")
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    ALICEtext2Dtwoplots("thesis",2)

    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);

    if log == "log":
        outname = os.path.join(folder, "{0}_pp_13.6TeV_{1}_material_budget_eta_vs_phi_r{2}_log_{3}_{4}.png".format(date, period, rid, cutname, suffix))
        c1.SaveAs(outname);

    elif log == "":
        outname = os.path.join(folder, "{0}_pp_13.6TeV_{1}_material_budget_eta_vs_phi_r{2}_{3}_{4}.png".format(date, period, rid, cutname, suffix))
        c1.SaveAs(outname);

    if rid == 3:
        find_phi_for_eta(h2data_complete,0, os.path.join(folder, "{0}_2D_Eta_vs_Phi_phi_values_for_r3.txt".format(date)));

# if __name__ == "__main__":
#     cutname = "qc"
#     period_mc = "LHC23d1k";
#     period_data = "LHC22f"
#     suffix = "AnyTrack";
#     filename_data = "/Users/alicamarieenderich/20231028_material_budget_plots/20231028_material_budget_eta_vs_phi_data_pp_13.6TeV_LHC22fAnyTrack.root" 
#     filename_mc = "/Users/alicamarieenderich/20231028_material_budget_plots/20231028_material_budget_eta_vs_phi_mc_pp_13.6TeV_LHC22fAnyTrack.root" 
#     date ="this_thesis" #datetime.date.today().strftime("%Y%m%d");
#     folder = "/Users/alicamarieenderich/{0}_material_budget_plots/".format(date);  
#     os.makedirs(folder, exist_ok=True);
#     for r in range(7):
#         for log in ["log", ""]:
#                 draw_material_z_vs_phi(folder, filename_data, filename_mc, cutname, period_data, r, suffix, log, date);