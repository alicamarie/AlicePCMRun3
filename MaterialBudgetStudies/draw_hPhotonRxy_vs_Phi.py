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
from ROOT import gStyle, gROOT, gSystem, gPad
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan
gStyle.SetOptStat(0);
gStyle.SetOptTitle(0);
gStyle.SetPaintTextFormat(".1f")
from FomattingMaterialBudget import FrameSettings, ALICEtext, RatioLegendSettings, FrameSettingsRatio, FrameSettings2D, ALICEtext2D

period_data     = "LHC22f";
period_mc       = "LHC23d1k";
filename_mc = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_155278_LHC23d1k.root"
# filename_mc = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_LHC23d1k_125889.root"
filename_data = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_148193_LHC22f_apass4.root"
# filename_data = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_LHC22f4_new_125184.root"
cutname         = "analysis_wo_mee";
# date = datetime.date.today().strftime("%Y%m%d");
# folder = "/Users/alicamarieenderich/{0}_material_budget_plots/".format(date);  
# os.makedirs(folder, exist_ok=True);
# outname         = os.path.join(folder,"{0}_PhotonPhivsRxy.root".format(date));
# outfile         = TFile(outname, "RECREATE");

# filename_plots = outname
# file_plots = TFile.Open(filename_plots, "READ");
# file_plots.ls()

r_bins          = [0, 14, 30, 42, 58, 69, 90]#, 180]
arr_rxy         = r_bins

rootfile_data = TFile.Open(filename_data, "READ");
rootfile_mc   = TFile.Open(filename_mc  , "READ");

rootdir_mc_gen  = rootfile_mc.Get("material-budget-mc")
list_gen        = rootdir_mc_gen.Get("Generated");
list_rec        = rootdir_mc_gen.Get("V0")
list_ev_gen     = rootdir_mc_gen.Get("Event");
list_ev_mc_gen  = list_ev_gen.FindObject("PCMPCM");
rootdir_mc_rec  = rootfile_mc.Get("pcm-qc-mc");
list_generated = rootdir_mc_rec.Get("Generated")
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

#________________________________________________

def plotting_Rxy(histogram, xtitle, ytitle, log , mctype, suffix, folder, date):
    c1 = TCanvas("c0","c0",0,0,800,800);
    c1.Divide(1,2,1e-5,1e-5); 
    c1.SetMargin(0.13,0.2,0.13,0.13);
    c1.SetTicks(1,1);

    frame1 = c1.DrawFrame(-10, -10, 10, 10);
    frame1.GetXaxis().SetTitle(xtitle);
    frame1.GetYaxis().SetTitle(ytitle);
    frame1.GetXaxis().SetTitleSize(0.05);
    frame1.GetYaxis().SetTitleSize(0.05);
    frame1.GetXaxis().SetTitleOffset(1.0);
    frame1.GetYaxis().SetTitleOffset(1.3);
    frame1.GetXaxis().SetLabelSize(0.05);
    frame1.GetYaxis().SetLabelSize(0.05);
    frame1.GetYaxis().SetMaxDigits(2);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    ROOT.SetOwnership(frame1,False);

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

    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    # date = datetime.date.today().strftime("%Y%m%d");
    c1.SaveAs(os.path.join(folder,"{0}_material_budget_hPhotonRxy_{1}_{2}_{3}.pdf".format(date, mctype, log, suffix)));
    c1.Close();

#________________________________________________

def plotting_Phi_vs_Rxy(histogram, xtitle, ytitle, log , mctype, suffix, folder, date):

    c1 = TCanvas("c0","c0",0,0,800,800);
    c1.Divide(1,2,1e-5,1e-5); 
    c1.SetMargin(0.13,0.2,0.13,0.13);
    c1.SetTicks(1,1);

    frame1 = c1.DrawFrame(-10, -10, 10, 10);
    frame1.GetXaxis().SetTitle(xtitle);
    frame1.GetYaxis().SetTitle(ytitle);
    FrameSettings2D(frame1)

    if log == "log":
        gPad.SetLogz()
      
    histogram.Draw("COLZ");

    for r in range(1, len(r_bins)-1):
        line1 = TLine(0,r_bins[r],TMath.TwoPi(),r_bins[r]);
        line1.SetLineColor(kRed);
        line1.SetLineStyle(1);
        line1.SetLineWidth(2);
        line1.Draw("");
        ROOT.SetOwnership(line1,False);

    txt = TPaveText(0.0,0.95,1.0,0.92,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(22);#centered,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.05);
    txt.AddText("M.C. {0}. #gamma (LHC23d1k)".format(mctype));
    txt.Draw();
    ROOT.SetOwnership(txt,False);
    ALICEtext2D("simulation")

    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    # date = datetime.date.today().strftime("%Y%m%d");
    c1.SaveAs(os.path.join(folder,"{0}_material_budget_hPhotonPhivsRxy_{1}_{2}_{3}.png".format(date, mctype, log, suffix)));
    c1.Close();

#________________________________________________
def plotting_Phi_for_R(histogram, xtitle, ytitle, log , mctype, suffix, ir, folder, date):

        make_common_style(histogram, 20, 1.0, kRed+1, 1, 0);
        ROOT.SetOwnership(histogram, False);

        ymax = histogram.GetMaximum() * 1.7;
        ymin = histogram.GetMinimum() * -0.2;

        c1 = TCanvas("c0","c0",0,0,800,800);
        c1.Divide(1,2,1e-5,1e-5); 
        c1.SetMargin(0.13,0.2,0.13,0.13);
        c1.SetTicks(1,1);

        frame1 = c1.DrawFrame(0, ymin, TMath.TwoPi(), ymax);
        frame1.SetTitle("conversion point #it{{#varphi}} in {0:3.2f} < #it{{R}}_{{xy}} cm < {1:3.2f} cm".format(r_bins[ir], r_bins[ir+1]));
        frame1.GetXaxis().SetTitle(xtitle);
        frame1.GetYaxis().SetTitle(ytitle);
        FrameSettings2D(frame1)

        histogram.Draw("E0Hsame");

        txt = TPaveText(0.0,0.95,1.0,0.92,"NDC");
        txt.SetFillColor(kWhite);
        txt.SetFillStyle(0);
        txt.SetBorderSize(0);
        txt.SetTextAlign(22);#centered,left
        txt.SetTextFont(42);#helvetica
        txt.SetTextSize(0.03);
        txt.AddText("conversion point #it{{#varphi}} in {0:3.2f} cm < #it{{R}}_{{xy}} < {1:3.2f} cm".format(r_bins[ir], r_bins[ir+1]))
        txt.Draw();
        ROOT.SetOwnership(txt,False);
        ALICEtext2D("simulation")

        leg = TLegend(0.25,0.6,0.7,0.9);
        leg.SetBorderSize(0);
        leg.SetFillColor(kWhite);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.05);
        leg.AddEntry(histogram, "M.C. {0}. #gamma (LHC23d1k)".format(mctype),"LP");
        leg.Draw("");
        ROOT.SetOwnership(leg,False);
        # date = datetime.date.today().strftime("%Y%m%d");
        c1.Modified();
        c1.Update();
        c1.SaveAs(os.path.join(folder,"{0}_material_budget_hPhoton_vs_phi_r{1}_{2}_{3}.png".format(date, ir, mctype, suffix)));
        ROOT.SetOwnership(c1,False);
        c1.Close();
#________________________________________________

def draw_hPhotonRxy_vs_Phi_mc_gen(suffix, log, outfile, gen, folder, date):
   
    hPhotonRxy = list_generated.FindObject("hPhotonRxy").Clone("hPhotonRxy_gen");
    hPhotonRxy.SetDirectory(0);

    hPhotonRxy.Scale(1/nev_gen);
    hPhotonRxy.Scale(1/nch_gen);
    hPhotonRxy.SetZTitle("#frac{1}{<#it{N}_{ch}^{PV}>} #frac{1}{#it{N}_{ev}} #it{N}_{#gamma}");
    hPhotonRxy.GetZaxis().SetTitleOffset(1.9);

    outfile.WriteTObject(hPhotonRxy);

    plotting_Rxy(hPhotonRxy, "conversion point #it{#varphi} (rad.)", "#frac{1}{<#it{N}_{ch}^{PV}>} #frac{1}{#it{N}_{ev}} #frac{d^{3}#it{N}_{#gamma}}{d#it{R}_{xy} d#it{#eta} d#it{#varphi}} (cm #upoint rad.)^{#minus1}", log , gen, suffix, folder, date)

    hPhotonPhivsRxy = list_generated.FindObject("hPhotonPhivsRxy").Clone("hPhotonPhivsRxy_{0}".format(gen));
    hPhotonPhivsRxy.SetDirectory(0);
    hPhotonPhivsRxy.RebinX(5)
    hPhotonPhivsRxy.Scale(1., "width")
    hPhotonPhivsRxy.Scale(1/nev_gen);
    hPhotonPhivsRxy.Scale(1/nch_gen);
    hPhotonPhivsRxy.SetZTitle("#frac{1}{<N_{ch}^{PV}>} #frac{1}{N_{ev}} #frac{d^{2}N_{#gamma}}{d#it{R}_{xy} d#it{#varphi}} ");
    hPhotonPhivsRxy.GetZaxis().SetTitleOffset(1.9);
    outfile.WriteTObject(hPhotonPhivsRxy);

    plotting_Phi_vs_Rxy(hPhotonPhivsRxy, "conversion point #it{#varphi} (rad.)", "#frac{1}{<#it{N}_{ch}^{PV}>} #frac{1}{#it{N}_{ev}} #frac{d^{2}N_{#gamma}}{d#it{R}_{xy} d#it{#varphi}} (cm #upoint rad.)^{#minus1}", log , gen, suffix, folder, date)

    #loop over r_bins
    for ir in range(0, len(r_bins)-1):
        r1 = r_bins[ir];
        r2 = r_bins[ir+1];
        dr = r2 - r1;
        bin_r1 = hPhotonPhivsRxy.GetYaxis().FindBin(r1 + 1e-6);
        bin_r2 = hPhotonPhivsRxy.GetYaxis().FindBin(r2 - 1e-6);

        hPhoton_ir = hPhotonPhivsRxy.ProjectionX("hPhotonPhi_r{0:d}_gen".format(ir), bin_r1, bin_r2, "");
        hPhoton_ir.SetTitle("conversion point #it{{#varphi}} in {0:3.2f} < #it{{R}}_{{xy}} < {1:3.2f} cm".format(r1, r2));
        hPhoton_ir.SetXTitle("#varphi (rad.)");
        hPhoton_ir.SetYTitle("#frac{1}{<N_{ch}^{PV}>} #frac{1}{N_{ev}} #frac{d^{2}N_{#gamma}}{d#it{R}_{xy} d#it{#varphi}}");
        hPhoton_ir.Scale(1./dr);
        #hPhoton_ir.Scale(1, "width");
        outfile.WriteTObject(hPhoton_ir);
        ROOT.SetOwnership(hPhoton_ir,False);
        # plotting_Phi_for_R(hPhoton_ir, "#varphi (rad.)","#frac{1}{<N_{ch}^{PV}>} #frac{1}{N_{ev}} #frac{d^{2}N_{#gamma}}{d#it{R}_{xy} d#it{#varphi}}", log , gen, suffix, ir)
#________________________________________________

def draw_hPhotonRxy_vs_Phi_mc_rec(suffix, log, outfile, gen, folder, date):
    #rootdir_mat_bud_rec = rootdir_mc_gen.Get("VO")
    # list_rec.ls()
    mat_bud_rec_qc = list_rec.FindObject("qc")
    hs_mc_rec = mat_bud_rec_qc.FindObject("hs_conv_point").Clone("hs_rec");

    hPhotonRxy = list_cut_mc_rec.FindObject("hGammaRxy_recalc");
    if hPhotonRxy:
        hPhotonRxy = list_cut_mc_rec.FindObject("hGammaRxy_recalc").Clone("hPhotonRxy_rec");
        hPhotonRxy.SetDirectory(0);
        
        hPhotonRxy.Scale(1/nev_rec);
        hPhotonRxy.Scale(1/nch_rec);
        hPhotonRxy.SetZTitle("#frac{1}{<#it{N}_{ch}^{PV}>} #frac{1}{#it{N}_{ev}} #it{N}_{#gamma}");
        hPhotonRxy.GetZaxis().SetTitleOffset(1.9);
        outfile.WriteTObject(hPhotonRxy);

        plotting_Rxy(hPhotonRxy, "conversion point #it{#varphi} (rad.)", "#frac{1}{<#it{N}_{ch}^{PV}>} #frac{1}{#it{N}_{ev}} #frac{d^{3}#it{N}_{#gamma}}{d#it{R}_{xy} d#it{#eta} d#it{#varphi}} (cm #upoint rad.)^{#minus1}", log , gen, suffix, folder, date)

    #hs = list_cut_mc_rec.FindObject("hs_conv_point").Clone("h_conv_point_rec");   
    #outfile.WriteTObject(hs);


    hPhotonPhivsRxy = hs_mc_rec.Projection(1, 2, "");
    #hPhotonPhivsRxy = list_cut_mc_rec.FindObject("hGammaRxy_recalc").Clone("h_conv_point_rec"); 
    hPhotonPhivsRxy.SetName("hPhotonPhivsRxy_rec");
    hPhotonPhivsRxy.SetTitle("conversion point #it{{R}}_{{xy}} vs. #it{{#varphi}}");
    hPhotonPhivsRxy.SetXTitle("#varphi (rad.)");
    hPhotonPhivsRxy.SetYTitle("R_{xy} (cm)");
    hPhotonPhivsRxy.SetZTitle("#frac{1}{<N_{ch}^{PV}>} #frac{1}{N_{ev}} #frac{d^{2}N_{#gamma}}{d#it{R}_{xy} d#it{#varphi}} ");
    hPhotonPhivsRxy.GetZaxis().SetTitleOffset(1.9);
   # hPhotonPhivsRxy.Scale(1, "width");
    hPhotonPhivsRxy.Scale(1./nev_rec);
    hPhotonPhivsRxy.Scale(1./nch_rec);
    hPhotonPhivsRxy.Sumw2();
    outfile.WriteTObject(hPhotonPhivsRxy);
    hPhotonPhivsRxy.GetYaxis().SetRangeUser(0,90);

    plotting_Phi_vs_Rxy(hPhotonPhivsRxy, "conversion point #it{#varphi} (rad.)", "#frac{1}{<#it{N}_{ch}^{PV}>} #frac{1}{#it{N}_{ev}} #frac{d^{2}N_{#gamma}}{d#it{R}_{xy} d#it{#varphi}} (cm #upoint rad.)^{#minus1}", log , gen, suffix, folder, date)

    #loop over r_bins
    # for ir in range(0, len(r_bins)-1):
    #     r1 = r_bins[ir];
    #     r2 = r_bins[ir+1];
    #     dr = r2 - r1;
    #     bin_r1 = hPhotonPhivsRxy.GetYaxis().FindBin(r1 + 1e-6);
    #     bin_r2 = hPhotonPhivsRxy.GetYaxis().FindBin(r2 - 1e-6);

    #     hPhoton_ir = hPhotonPhivsRxy.ProjectionX("hPhotonPhi_r{0:d}_rec".format(ir), bin_r1, bin_r2, "");
    #     hPhoton_ir.SetTitle("conversion point #it{{#varphi}} in {0:3.2f} cm < #it{{R}}_{{xy}} < {1:3.2f} cm".format(r1, r2));
    #     hPhoton_ir.SetXTitle("#varphi (rad.)");
    #     hPhoton_ir.SetYTitle("#frac{1}{<N_{ch}^{PV}>} #frac{1}{N_{ev}} #frac{d^{2}N_{#gamma}}{d#it{R}_{xy} d#it{#varphi}}");
    #     hPhoton_ir.Scale(1./dr);
    #     #hPhoton_ir.Scale(1, "width");
    #     outfile.WriteTObject(hPhoton_ir);
    #     ROOT.SetOwnership(hPhoton_ir,False);

        # plotting_Phi_for_R(hPhoton_ir, "#varphi (rad.)", "#frac{1}{<N_{ch}^{PV}>} #frac{1}{N_{ev}} #frac{d^{2}N_{#gamma}}{d#it{R}_{xy} d#it{#varphi}}", log , gen, suffix, ir)

def plot_ratio(i, date):

    # date = datetime.date.today().strftime("%Y%m%d");
    folder = "/Users/alicamarieenderich/{0}_material_budget_plots/".format(date);  
    os.makedirs(folder, exist_ok=True);
    outname         = os.path.join(folder,"{0}_PhotonPhivsRxy.root".format(date));
    outfile         = TFile(outname, "RECREATE");

    filename_plots = outname
    file_plots = TFile.Open(filename_plots, "READ");
    file_plots.ls()


    h1_rec_name = "hPhotonPhi_r{}_rec".format(i)
    h1rec = file_plots.Get(h1_rec_name).Clone("h1rec");
    h1_gen_name = "hPhotonPhi_r{}_gen".format(i)
    h1gen = file_plots.Get(h1_gen_name).Clone("h1gen");
    make_common_style(h1rec, 20, 1.0, kBlue+1, 1, 0);
    make_common_style(h1gen  , 20, 1.0, kGreen+2, 1, 0);
    ROOT.SetOwnership(h1rec, False);
    ROOT.SetOwnership(h1gen, False);

    ymax = max(h1rec.GetMaximum() , h1gen.GetMaximum()) * 1.7;
    ymin = max(h1rec.GetMinimum() , h1gen.GetMinimum()) * -0.2;

    c1 = TCanvas("c0","c0",0,0,800,800);
    c1.Divide(1,2,1e-3,1e-3);
    p1 = c1.cd(1);
    p1.SetPad(0,0.35,1,1);
    p1.SetMargin(0.15,0.02,0.0,0.1);
    p1.SetTicks(1,1);

    frame1 = p1.DrawFrame(0, ymin, TMath.TwoPi(), ymax);
    frame1.SetTitle("conversion point #it{{#varphi}} in {0:3.2f} cm < #it{{R}}_{{xy}} < {1:3.2f} cm".format(r_bins[i], r_bins[i+1]));
    frame1.GetXaxis().SetTitle("conversion point #it{#varphi} (rad.)");
    frame1.GetYaxis().SetTitle("#frac{1}{<#it{N}_{ch}^{PV}>} #frac{1}{#it{N}_{ev}} #frac{d^{2}#it{N}_{#gamma}}{d#it{R}_{xy} d#it{#varphi}} (cm #upoint rad.)^{#minus1}");
    frame1.GetXaxis().SetTitleSize(0.05);
    frame1.GetYaxis().SetTitleSize(0.05);
    frame1.GetXaxis().SetTitleOffset(1.0);
    frame1.GetYaxis().SetTitleOffset(1.3);
    frame1.GetXaxis().SetLabelSize(0.05);
    frame1.GetYaxis().SetLabelSize(0.05);
    frame1.GetYaxis().SetMaxDigits(3);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    ROOT.SetOwnership(frame1,False);

    h1rec.Draw("E0Hsame");
    h1gen.Draw("E0Hsame");

    txt = TPaveText(0.0,0.95,1.0,0.92,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(22);#centered,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.05);
    txt.AddText("conversion point #it{{#varphi}} in {0:3.2f} cm < #it{{R}}_{{xy}} < {1:3.2f} cm".format(r_bins[i], r_bins[i+1]))
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    # txt = TPaveText(0.2,0.72,0.4,0.92,"NDC");
    # txt.SetFillColor(kWhite);
    # txt.SetFillStyle(0);
    # txt.SetBorderSize(0);
    # txt.SetTextAlign(12);#middle,left
    # txt.SetTextFont(42);#helvetica
    # txt.SetTextSize(0.05);
    # txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    # txt.Draw();
    # ROOT.SetOwnership(txt,False);

    txt = TPaveText(0.7,0.8,0.77,0.85,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(32);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.02);
    txt.AddText("ALICE this thesis");
    txt.AddText("pp at #sqrt{#it{s}}  = 13.6 TeV")
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    leg = TLegend(0.2,0.6,0.4,0.7);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.05);
    leg.AddEntry(h1rec, "M.C. rec. #gamma (LHC23d1k)","LP");
    leg.AddEntry(h1gen  , "M.C. gen. #gamma (LHC23d1k)","LP");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    p2 = c1.cd(2);
    p2.SetPad(0,0,1,0.35);
    p2.SetMargin(0.15,0.02,0.22,0.0);
    p2.SetTicks(1,1);

    frame2 = p2.DrawFrame(0,0.,TMath.TwoPi(),1.1);
    frame2.GetXaxis().SetTitle("#it{#varphi} (rad.)");
    frame2.GetYaxis().SetTitle("#frac{rec}{gen}");
    frame2.GetXaxis().SetTitleSize(0.10);
    frame2.GetYaxis().SetTitleSize(0.10);
    frame2.GetXaxis().SetTitleOffset(1.0);
    frame2.GetYaxis().SetTitleOffset(0.6);
    frame2.GetXaxis().SetLabelSize(0.1);
    frame2.GetYaxis().SetLabelSize(0.1);
    frame2.GetYaxis().CenterTitle(True);
    ROOT.SetOwnership(frame2,False);

    h1ratio = h1rec.Clone("h1ratio");
    h1ratio.Reset();
    h1ratio.Sumw2();
    h1ratio.Divide(h1rec, h1gen, 1., 1., "G");
    h1ratio.Draw("E0same");
    h1ratio.SetDirectory(0);
    ROOT.SetOwnership(h1ratio,False);

    line1 = TLine(0,1,TMath.TwoPi(),1);
    line1.SetLineColor(kBlack);
    line1.SetLineStyle(1);
    line1.SetLineWidth(2);
    line1.Draw("");
    ROOT.SetOwnership(line1,False);

    # date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs(os.path.join(folder,"{0}_material_budget_rec_gen_r{1}.pdf".format(date, i)));
    c1.Close();

# #________________________________________________
if __name__ == "__main__":
    period_data     = "LHC22f";
    period_mc       = "LHC23d1k";
    filename_mc = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_155278_LHC23d1k.root"
    # # filename_mc = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_LHC23d1k_125889.root"
    # filename_data = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_148193_LHC22f_apass4.root"
    # # filename_data = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_LHC22f4_new_125184.root"

    filename_data = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_155756_LHC22f_pass4.root"
    cutname         = "analysis_wo_mee";
    date = "this_thesis" #datetime.date.today().strftime("%Y%m%d");
    folder = "/Users/alicamarieenderich/{0}_material_budget_plots/".format(date);  
    os.makedirs(folder, exist_ok=True);
    outname         = os.path.join(folder,"{0}_PhotonPhivsRxy.root".format(date));
    outfile         = TFile(outname, "RECREATE");
    suffix = ""

    draw_hPhotonRxy_vs_Phi_mc_gen(suffix, "log", outfile, "gen", folder, date);
    # draw_hPhotonRxy_vs_Phi_mc_rec(suffix, "log", outfile, "rec", folder, date); 

    # for r in range(6):
    #     plot_ratio(r)