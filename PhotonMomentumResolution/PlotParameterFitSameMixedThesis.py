# This code is based on the code provided by Daiki Sekihata in his repository
# https://github.com/dsekihat/photon_sw_run3 for photon analysis in ALICE Run 3

import numpy as np
import yaml
import datetime
import ROOT
from ROOT import TFile, TDirectory, THashList, TH1F, TH2F, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TGaxis, TPad, TLatex, TBox, TString, TLine
from ROOT import gStyle, gROOT, gSystem
from ROOT import TFile, TDirectory, THashList, TF1, TH1F, TCanvas, TH1
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kFullCircle, kOrange, kGray, kTRUE
import re
import numpy as np
import datetime
import math
import ctypes
import os, sys, shutil

from HistoFormatting import FrameSettings, CanvasSettings, PadSettings, DrawHisto, SetTitle, SetStyleTLatex

def make_common_style(g1,marker,size,color,width=1,fill=0):
    g1.SetMarkerStyle(marker);
    g1.SetMarkerColor(color);
    g1.SetMarkerSize(size);
    g1.SetLineColor(color);
    g1.SetLineWidth(width);
    g1.SetFillColorAlpha(color, 0.65);
    g1.SetFillStyle(fill);

def draw_comparison_parameters(list_data0, index, fit_variable, date, period_str):
        ymin = [item.GetMinimum() for item in list_data0]
        ymax= [item.GetMaximum() for item in list_data0]
        yMin = 0.;
        yMax = 0.;
        for i in range(1):
            if isinstance(list_data0[0], TH1):
                for j in range ( list_data0[i].GetXaxis().FindBin(0.125), list_data0[i].GetXaxis().FindBin(0.138)):
                    if list_data0[i].GetBinContent(j) < yMin:
                        yMin = list_data0[i].GetBinContent(i);
                    if list_data0[i].GetBinContent(j) > yMax:
                        yMax = list_data0[i].GetBinContent(j);

    #canvas plotting
        c1 = TCanvas("c0","c0",0,0,800,800);
        c1.Divide(1,1,1e-3,1e-3);
        p1 = c1.cd(1);
        p1.SetPad(0,0., 1,1);
        p1.SetMargin(0.13,0.03,0.13,0.13);
        p1.SetTicks(1,1);
        if list_data0[0].GetTitle() == "raw yield":
            p1.SetLogy()

        frame1 = p1.DrawFrame(0, -2., 0.3, 1.2*yMax);
        frame1.GetXaxis().SetTitle("#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})"); 
        frame1.GetYaxis().SetTitle("Counts");
        frame1.GetXaxis().SetTitleSize(0.045);
        frame1.GetYaxis().SetTitleSize(0.045);
        frame1.GetXaxis().SetTitleOffset(1.2);
        frame1.GetYaxis().SetTitleOffset(1.2);
        frame1.GetXaxis().SetLabelSize(0.045);
        frame1.GetYaxis().SetLabelSize(0.045);
        frame1.GetYaxis().SetMaxDigits(3);
        frame1.GetXaxis().SetLabelOffset(0.01);
        frame1.GetYaxis().SetLabelOffset(0.01);
        frame1.GetYaxis().SetNdivisions(5, 5, 0, kTRUE)
        ROOT.SetOwnership(frame1,False);

        color = [kRed+1, kBlue+1, kGreen+2, kMagenta+2, kCyan+1]
        DrawHisto(list_data0[0], kGreen+2, "Hsame", 0.9, 24);       
        DrawHisto(list_data0[1], kGray+2, "Hsame", 0.9, 24);
        DrawHisto(list_data0[2], kBlue+1, "e,p,same", 0.9);
        list_data0[3].SetLineColor(kRed+1);
        list_data0[3].SetLineWidth(2);
        list_data0[3].SetFillColor(kRed+1);
        list_data0[3].SetFillStyle(0);
        list_data0[3].DrawCopy("same");
        if list_data0[0].GetTitle() == "mean":
                line2 = TLine(0,134.97,12,134.97);
                line2.SetLineColor(kGray +1);
                line2.SetLineStyle(2);
                line2.SetLineWidth(2);
                line2.Draw("");
                ROOT.SetOwnership(line2,False); 

        txt = TPaveText(0.,0.88,1.,0.995,"NDC");
        txt.SetFillColor(kWhite);
        txt.SetFillStyle(0);
        txt.SetBorderSize(0);
        txt.SetTextAlign(22);#middle,left
        txt.SetTextFont(42);#helvetica
        txt.SetTextSize(0.035);
        txt.AddText("Invariant mass distribution of reconstructed photon pairs");
        txt.AddText("0.4 GeV/c < #it{p}_{T}  < 0.8 GeV/c");
        txt.Draw();
        ROOT.SetOwnership(txt,False);

        txt = TPaveText(0.58,0.78,0.95,0.85,"NDC");
        txt.SetFillColor(kWhite);
        txt.SetFillStyle(0);
        txt.SetBorderSize(0);
        txt.SetTextAlign(32);#middle,left
        txt.SetTextFont(42);#helvetica
        txt.SetTextSize(0.03);
        txt.AddText("this thesis");
        txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
        txt.Draw();
        ROOT.SetOwnership(txt,False);

        leg = TLegend(0.52,0.63,0.95,0.77);
        leg.SetBorderSize(0);
        leg.SetFillColor(kWhite);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.03);
        leg.SetTextAlign(32);
        leg.SetTextFont(42);#helvetica
        cut = ["same evt. #it{M}_{#gamma #gamma}", "mixed evt. #it{M}_{#gamma #gamma} (scaled)", "same evt. #it{M}_{#gamma #gamma} (Signal-BG)", 'asymmetric Gaussian fit']
        for i in range(len(list_data0)-1):
            leg.AddEntry(list_data0[i], "{0}".format(cut[i]),"LP");
        leg.AddEntry(list_data0[len(list_data0)-1], "{0}".format(cut[len(list_data0)-1]), "L");
        leg.Draw("");
        ROOT.SetOwnership(leg,False);

        leg = TLegend(0.15,0.78,0.40,0.85);
        leg.SetBorderSize(0);
        leg.SetFillColor(kWhite);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.03);
        leg.SetTextAlign(12);
        leg.SetTextFont(42);#helvetica
        leg.AddEntry(list_data0[2], "{}".format(fit_variable),"LP");
        leg.Draw("");
        ROOT.SetOwnership(leg,False);

        c1.Modified();
        c1.Update();
        ROOT.SetOwnership(c1,False);
        if fit_variable == "fwhm/2.36":
             fit_variable = "fwhm"
        filepath = os.path.join(folder, "{0}_Parameter_Comparison_{1}{2}.pdf".format(date, period_str, fit_variable));    
        c1.SaveAs(filepath);

def run(filename_data0,config,type,folder, period_str):
    print(sys._getframe().f_code.co_name);
    arr_pt = np.array(config["common"]["pt_bin"],dtype=float);
    print("pT binning = ",arr_pt);
    print("type = ",type);
    print("input = ",filename_data0);
    rootfile_data0 = TFile.Open(filename_data0,"READ");
    meson = config["common"]["meson"];

    list_fit_func = config[type]["fit_func"];

    list_fit_min = config["common"]["fit_min"];
    list_fit_max = config["common"]["fit_max"];
    if len(list_fit_min) != len(list_fit_max):
        return;

    list_integral_min = config["common"]["integral_min"];
    list_integral_max = config["common"]["integral_max"];
    if len(list_integral_min) != len(list_integral_max):
        return;

    list_yield_min = config["common"]["yield_min"];
    list_yield_max = config["common"]["yield_max"];
    if len(list_yield_min) != len(list_yield_max):
        return;

    nsys = len(config[type]['subsystems']);
    print("nsys",nsys); 

    def loop_data(rootfile): 
            for isys in range(0,nsys):
                ssname = config[type]['subsystems'][isys]['name']; #subsystem name
                print("plot subsystem", ssname);
                cuts = config[type]["subsystems"][isys]['cuts']
                cutnames = [cut['name'] for cut in cuts]
                print("cutnames", cutnames); 
                nc = len(cutnames);
                nfit = len(list_fit_func)
                same_list = [];
                mixed_list = [];
                histo_list = [];
                fit_list = [];
                cut_list = [];
                list_ss   = rootfile.Get(ssname);
                print("nc", nc)
                for ic in range(0,nc):
                    cutname = cutnames[ic];
                    print(cutname)
                    list_ss_cut = list_ss.FindObject(cutname);
                    print("cutname ", cutname)        
                    cut_list.append(cutname)     
                    for ifit in range(0, nfit): 
                        fitname = "gausexplinear";
                        list_fitname  = list_ss_cut.FindObject(fitname);
                        print("fitname ", fitname)
                        list_plot = list_fitname.FindObject("fit_0.04_0.20_GeVc2")
                        histo_list.append(list_plot.FindObject("h1mgg_pt0"));
                        fit_list.append (list_plot.FindObject("f1total_pt0"));
                        mixed_list.append(list_plot.FindObject("h1mgg_mix_scaled_pt0"));
                        same_list.append(list_plot.FindObject("h1mgg_same_pt0"))          
            return same_list, mixed_list, histo_list, fit_list, cut_list

    same_list_data, mixed_list_data, histo_list_data, fit_list_data, cut_list  = loop_data(rootfile_data0)

    list_data0 = [same_list_data, mixed_list_data, histo_list_data, fit_list_data]
    list_data0_transposed = [list(t) for t in zip(*list_data0)]
    for i in range(4):
        draw_comparison_parameters(list_data0_transposed[i], i, cut_list[i], date, period_str)

if __name__ == "__main__":
    period_array = ["LHC22f"]
    filename_array = ["/Users/alicamarieenderich/AnalysisResults/AnalysisResults_124838_LHC22f_pass4.root"]
    type_array = ["data"]
    config_array = ["/Users/alicamarieenderich/202312_invariant_mass_code/configs/config_pp_13.6TeV_pi0_LHC22f.yml"]
    period_str = ""
    period_str += period_array[0]
    period_str += "_"
    period = period_array[0]
    filename = filename_array[0]
    cutname = "qc"
    suffix = "AnyTrack";
    type = type_array[0]
    config_file = config_array[0]
    with open(config_file, "r", encoding="utf-8") as config_yml:
        config = yaml.safe_load(config_yml)
    # Date or prefix "this_thesis"
    date = "this_thesis" #datetime.date.today().strftime("%Y%m%d");
    folder = "/Users/alicamarieenderich/{0}_{1}_invariant_mass_plots/".format(date, period);  
    os.makedirs(folder, exist_ok=True);
    # input files from fitting process
    filename_data0 = "/Users/alicamarieenderich/this_thesis_LHC22f_invariant_mass_plots/this_thesis_LHC22f_pi0_data_ptspectrum_pp_13.6TeV_LHC22f_AnyTrack.root"
    filename_data1 = "/Users/alicamarieenderich/this_thesis_LHC23zc_invariant_mass_plots/this_thesis_LHC23zc_pi0_data_ptspectrum_pp_13.6TeV_LHC22f_AnyTrack.root"
    filename_data2 = "/Users/alicamarieenderich/this_thesis_LHC22o_minBias_invariant_mass_plots/this_thesis_LHC22o_minBias_pi0_data_ptspectrum_pp_13.6TeV_LHC22f_AnyTrack.root"

    run(filename_data0, config, type, folder, period_str)