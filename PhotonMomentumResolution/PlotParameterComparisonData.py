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
gStyle.SetOptStat(0);
gStyle.SetOptTitle(0);
gStyle.SetPadTickX(1);
gStyle.SetPadTickY(1);
gStyle.SetErrorX(0)
gStyle.SetEndErrorSize(5)
ROOT.TGaxis.SetMaxDigits(2)

from HistoFormatting import FrameSettings, CanvasSettings, PadSettings, DrawHisto, SetTitle, SetStyleTLatex

def make_common_style(g1,marker,size,color,width=1,fill=0):
    g1.SetMarkerStyle(marker);
    g1.SetMarkerColor(color);
    g1.SetMarkerSize(size);
    g1.SetLineColor(color);
    g1.SetLineWidth(width);
    g1.SetFillColorAlpha(color, 0.65);
    g1.SetFillStyle(fill);

def draw_comparison_parameters(list_data0, list_data1, list_data2, index, fit_variable, date, period_str):
        ymin = [125,  1*1e-3, 1e-8]
        ymax = [145, 25*1e-3, 1e-2]
        for i in range(len(list_data0)):
            print("MAX", list_data0[i].GetMaximum(),list_data1[i].GetMaximum(), list_data2[i].GetMaximum())
    #canvas plotting
        c1 = TCanvas("c0","c0",0,0,800,800);
        c1.Divide(1,2,1e-3,1e-3);
        p1 = c1.cd(1);
        p1.SetPad(0,0.3,1,1);
        p1.SetMargin(0.15,0.02,0.,0.15);
        p1.SetTicks(1,1);
        if list_data0[0].GetTitle() == "raw yield":
            p1.SetLogy()

        frame1 = p1.DrawFrame(0, ymin[index], 10., ymax[index]);
        frame1.GetXaxis().SetTitle("#it{p}_{T} (GeV/c)");
        title_list = ["mean (MeV/#it{c}^{2})", "fwhm/2.36 (GeV/#it{c}^{2})", "#frac{1}{#it{N}_{ev}}#frac{d#it{N}}{d#it{p}_{T}} (GeV/#it{c})^{-1}"]
        if list_data0[0].GetTitle() == "raw yield":
             yTitle = title_list[2]
        if list_data0[0].GetTitle() == "fwhm/2.36":
             yTitle = title_list[1]
        if list_data0[0].GetTitle() == "mean":
             yTitle = title_list[0]
        frame1.GetYaxis().SetTitle(yTitle);
        frame1.GetXaxis().SetTitleSize(0.045);
        frame1.GetYaxis().SetTitleSize(0.045);
        frame1.GetXaxis().SetTitleOffset(1.0);
        frame1.GetYaxis().SetTitleOffset(1.52);
        frame1.GetXaxis().SetLabelSize(0.045);
        frame1.GetYaxis().SetLabelSize(0.045);
        frame1.GetYaxis().SetMaxDigits(3);
        frame1.GetXaxis().SetLabelOffset(0.01);
        frame1.GetYaxis().SetLabelOffset(0.01);
        frame1.GetYaxis().SetNdivisions(5, 5, 0, kTRUE)
        ROOT.SetOwnership(frame1,False);

        color = [kRed+1, kBlue+1, kGreen+2, kMagenta+2, kCyan+1]
        for i in range(len(list_data0)):
                print(i, list_data0[i].GetTitle(), len(list_data0))
                make_common_style(list_data0[i], 20, 1.5, color[0], 1, 0);
                make_common_style(list_data1[i], 21, 1.5, color[1], 1, 0); #0.9
                make_common_style(list_data2[i], 22, 1.5, color[2], 1, 0);

                if list_data0[0].GetTitle() == "mean":
                    list_data0[i].Scale(1000.)
                    list_data1[i].Scale(1000.)
                    list_data2[i].Scale(1000.)
                list_data0[i].Draw("E0same");
                list_data1[i].Draw("E0same");
                list_data2[i].Draw("E0same");

        if list_data0[0].GetTitle() == "mean":
                line2 = TLine(0,134.97,12,134.97);
                line2.SetLineColor(kGray +1);
                line2.SetLineStyle(2);
                line2.SetLineWidth(2);
                line2.Draw("");
                ROOT.SetOwnership(line2,False); 

        txt = TPaveText(0.,0.9,1.,0.95,"NDC");
        txt.SetFillColor(kWhite);
        txt.SetFillStyle(0);
        txt.SetBorderSize(0);
        txt.SetTextAlign(22);#middle,left
        txt.SetTextFont(42);#helvetica
        txt.SetTextSize(0.045);
        txt.AddText("{0} for the cut 'qc' in multiple datasets".format(fit_variable));
        txt.Draw();
        ROOT.SetOwnership(txt,False);

        txt = TPaveText(0.50,0.65,0.95,0.82,"NDC");
        txt.SetFillColor(kWhite);
        txt.SetFillStyle(0);
        txt.SetBorderSize(0);
        txt.SetTextAlign(32);#middle,left
        txt.SetTextFont(42);#helvetica
        txt.SetTextSize(0.045);#0.03
        txt.AddText("ALICE this thesis");
        txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
        txt.AddText("Data #gamma candidates for")
        txt.Draw();
        ROOT.SetOwnership(txt,False);

        leg = TLegend(0.60,0.50,0.94 ,0.65);
        leg.SetBorderSize(0);
        leg.SetFillColor(kWhite);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.045); #0.03
        leg.SetTextAlign(12);
        leg.SetTextFont(42);#helvetica
        # leg.AddEntry(list_data0[0], "Data #gamma candidates (LHC22f)","LP");
        # leg.AddEntry(list_data1[0] , "Data #gamma candidates (LHC23zc)","LP");
        # leg.AddEntry(list_data2[0] , "Data #gamma candidates (LHC22o min Bias)","LP");
        leg.AddEntry(list_data0[0], "LHC22f","LP");
        leg.AddEntry(list_data1[0] , "LHC23zc","LP");
        leg.AddEntry(list_data2[0] , "LHC22o min Bias","LP");
        leg.Draw("");
        ROOT.SetOwnership(leg,False);

        p2 = c1.cd(2);
        p2.SetPad(0,0,1,0.3);
        p2.SetMargin(0.15,0.02,0.22,0.0);
        p2.SetTicks(1,1);
        cut_in_ratio = 0.5
        ymin, ymax = 0.95, 1.05
        if index ==2:
            ymin, ymax = 8*1e-3, 2*1e0
            p2.SetLogy()
        if index ==1:
            ymin, ymax = 0.5,2.0

        frame2 = p2.DrawFrame(0,ymin,10.,ymax);
        frame2.GetXaxis().SetTitle("#it{p}_{T} (GeV/c)");
        frame2.GetYaxis().SetTitle("#frac{dataset}{LHC22f}");
        frame2.GetYaxis().SetNdivisions(5, 5, 0, kTRUE)
        frame2.GetXaxis().SetTitleSize(0.10);
        frame2.GetYaxis().SetTitleSize(0.10);
        frame2.GetXaxis().SetTitleOffset(1.0);
        frame2.GetYaxis().SetTitleOffset(0.7);
        frame2.GetXaxis().SetLabelSize(0.10);
        frame2.GetYaxis().SetLabelSize(0.10);
        frame2.GetYaxis().CenterTitle(True);
        frame2.GetXaxis().SetLabelOffset(0.01);
        frame2.GetYaxis().SetLabelOffset(0.01);
        ROOT.SetOwnership(frame2,False);

        color = [kRed+1, kBlue+1, kGreen+2, kMagenta+2, kCyan+1]
        for i in range(len(list_data0)):
                h1ratio = list_data1[i].Clone("h1ratio");
                h1ratio.Reset();
                h1ratio.Sumw2();
                h1ratio.Divide(list_data1[i], list_data0[i], 1., 1., "G");
                make_common_style(h1ratio, 21, 1.5, color[1], 1, 0); #0.9
                h1ratio.Draw("E0same");
                h1ratio.SetDirectory(0);
                ROOT.SetOwnership(h1ratio,False);
        for i in range(len(list_data0)):
                h1ratio = list_data2[i].Clone("h1ratio");
                h1ratio.Reset();
                h1ratio.Sumw2();
                h1ratio.Divide(list_data2[i], list_data0[i], 1., 1., "G");
                make_common_style(h1ratio, 22, 1.5, color[2], 1, 0); #0.9
                h1ratio.Draw("E0same");
                h1ratio.SetDirectory(0);
                ROOT.SetOwnership(h1ratio,False);
        line1 = TLine(0.,1,10.,1);
        line1.SetLineColor(kBlack);
        line1.SetLineStyle(1);
        line1.SetLineWidth(2);
        line1.Draw("");
        ROOT.SetOwnership(line1,False);

        c1.Modified();
        c1.Update();
        ROOT.SetOwnership(c1,False);
        if fit_variable == "fwhm/2.36":
             fit_variable = "fwhm"
        filepath = os.path.join(folder, "{0}_Parameter_Comparison_{1}{2}.pdf".format(date, period_str, fit_variable));    
        c1.SaveAs(filepath);

def run(filename_data0, filename_data1, filename_data2,config,type,folder, period_str):
    print(sys._getframe().f_code.co_name);
    arr_pt = np.array(config["common"]["pt_bin"],dtype=float);
    print("pT binning = ",arr_pt);
    print("type = ",type);
    print("input = ",filename_data0);
    rootfile_data0 = TFile.Open(filename_data0,"READ");
    rootfile_data1 = TFile.Open(filename_data1,"READ");
    rootfile_data2 = TFile.Open(filename_data2,"READ");
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
                list_parameters_comparison = [];
                yield_list = [];
                mean_list = [];
                fwhm_list = [];
                list_ss   = rootfile.Get(ssname);
                print("nc", nc)
                ic =1
                cutname = cutnames[ic];
                print(cutname)
                list_ss_cut = list_ss.FindObject(cutname);
                print("cutname ", cutname)             
                for ifit in range(0, nfit): 
                    fitname = "gausexplinear";
                    list_fitname  = list_ss_cut.FindObject(fitname);
                    print("fitname ", fitname)

                    list_plot = list_fitname.FindObject("fit_0.04_0.20_GeVc2")
                    yield_list.append(list_plot.FindObject("h1yield_param"))
                    mean_list.append(list_plot.FindObject("h1mean_param"))
                    fwhm_list.append(list_plot.FindObject("h1fwhm_param"))                        
            return mean_list, fwhm_list, yield_list  
            return mean_list, fwhm_list, yield_list  
    # def loop_mc(rootfile): 
    #         for isys in range(0,nsys):
    #             ssname = config[type]['subsystems'][isys]['name']; #subsystem name
    #             print("plot subsystem", ssname);
    #             cuts = config[type]["subsystems"][isys]['cuts']
    #             cutnames = [cut['name'] for cut in cuts]
    #             print("cutnames", cutnames); 
    #             nc = len(cutnames);
    #             nfit = len(list_fit_func)
    #             list_parameters_comparison = [];
    #             yield_list = [];
    #             list_ss   = rootfile.Get(ssname);
    #             for ic in range(0,nc):
    #                 cutname = cutnames[ic];
    #                 list_ss_cut = list_ss.FindObject(cutname);
    #                 print("cutname ", cutname)             
    #                 for ifit in range(0, nfit): 
    #                     fitname = "gausexp";
    #                     list_fitname  = list_ss_cut.FindObject(fitname);
    #                     print("fitname ", fitname)

    #                     list_plot = list_fitname.FindObject("fit_0.04_0.20_GeVc2")
    #                     yield_list.append(list_plot.FindObject("h1yield_param"))


    #                     parameter_list = [];
    #                     same_list = [];
    #                     mixed_list = [];
    #                     it = ROOT.TIter(list_plot)
    #                     obj = it.Next()
    #                     while obj:
    #                         objName = obj.GetName()
    #                         if "param" in objName:
    #                             parameter_list.append(obj)
    #                         obj = it.Next()
    #                     print("HERE")
    #                     for i in range(len(parameter_list)):
    #                         print(parameter_list[i].GetTitle())
    #     # pdf output of mass, amplitude and width for each fit and comparison of all cuts
    #                         parameter_comparison = [];
    #                         for i_parameter in range(5):
    #                             parameter_comparison.append(parameter_list[i_parameter]);
    #                         list_parameters_comparison.append(parameter_comparison);
    #         return list_parameters_comparison, yield_list     
            return mean_list, fwhm_list, yield_list   
    # def loop_mc(rootfile): 
    #         for isys in range(0,nsys):
    #             ssname = config[type]['subsystems'][isys]['name']; #subsystem name
    #             print("plot subsystem", ssname);
    #             cuts = config[type]["subsystems"][isys]['cuts']
    #             cutnames = [cut['name'] for cut in cuts]
    #             print("cutnames", cutnames); 
    #             nc = len(cutnames);
    #             nfit = len(list_fit_func)
    #             list_parameters_comparison = [];
    #             yield_list = [];
    #             list_ss   = rootfile.Get(ssname);
    #             for ic in range(0,nc):
    #                 cutname = cutnames[ic];
    #                 list_ss_cut = list_ss.FindObject(cutname);
    #                 print("cutname ", cutname)             
    #                 for ifit in range(0, nfit): 
    #                     fitname = "gausexp";
    #                     list_fitname  = list_ss_cut.FindObject(fitname);
    #                     print("fitname ", fitname)

    #                     list_plot = list_fitname.FindObject("fit_0.04_0.20_GeVc2")
    #                     yield_list.append(list_plot.FindObject("h1yield_param"))


    #                     parameter_list = [];
    #                     same_list = [];
    #                     mixed_list = [];
    #                     it = ROOT.TIter(list_plot)
    #                     obj = it.Next()
    #                     while obj:
    #                         objName = obj.GetName()
    #                         if "param" in objName:
    #                             parameter_list.append(obj)
    #                         obj = it.Next()
    #                     print("HERE")
    #                     for i in range(len(parameter_list)):
    #                         print(parameter_list[i].GetTitle())
    #     # pdf output of mass, amplitude and width for each fit and comparison of all cuts
    #                         parameter_comparison = [];
    #                         for i_parameter in range(5):
    #                             parameter_comparison.append(parameter_list[i_parameter]);
    #                         list_parameters_comparison.append(parameter_comparison);
    #         return list_parameters_comparison, yield_list     

    mean_list_data0,fwhm_list_data0, yield_list_data0 = loop_data(rootfile_data0)
    mean_list_data1,fwhm_list_data1, yield_list_data1 = loop_data(rootfile_data1)
    mean_list_data2,fwhm_list_data2, yield_list_data2 = loop_data(rootfile_data2)

    list_data0 = [mean_list_data0,fwhm_list_data0, yield_list_data0]
    list_data1 = [mean_list_data1,fwhm_list_data1, yield_list_data1]
    list_data2 = [mean_list_data2,fwhm_list_data2, yield_list_data2]
    print("length", len(list_data0), (len(mean_list_data0,)), (len(fwhm_list_data0)), (len(yield_list_data0)))
    print("length", len(list_data1), (len(mean_list_data1,)), (len(fwhm_list_data1)), (len(yield_list_data1)))
    print("length", len(list_data2), (len(mean_list_data2,)), (len(fwhm_list_data2)), (len(yield_list_data2)))
    for i in range(3):
        draw_comparison_parameters(list_data0[i], list_data1[i], list_data2[i], i, list_data0[i][0].GetTitle(), date, period_str)

period_array = ["LHC22f", "LHC23zc", "LHC22o_minBias"]

filename_array = ["/Users/alicamarieenderich/AnalysisResults/AnalysisResults_124838_LHC22f_pass4.root",
                  "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_134765_LHC23zc_pass1_relval_itstpcmap_1.root",
            "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_135860_LHC22o_pass4_minBias_medium.root"]

type_array = ["data", "data", "data"]

config_array = ["/Users/alicamarieenderich/202312_invariant_mass_code/configs/config_pp_13.6TeV_pi0_LHC22f.yml",
                "/Users/alicamarieenderich/202312_invariant_mass_code/configs/config_pp_13.6TeV_pi0_LHC23zc.yml",
                "/Users/alicamarieenderich/202312_invariant_mass_code/configs/config_pp_13.6TeV_pi0_LHC22o_medium.yml"]

period_str = ""
for i in range(len(period_array)):
    period_str += period_array[i]
    period_str += "_"
    print("period_str", period_str)

for i in range(len(period_array)):
    period = period_array[i]
    filename = filename_array[i]
    cutname = "qc"
    suffix = "AnyTrack";
    type = type_array[i]
    config_file = config_array[i]
    with open(config_file, "r", encoding="utf-8") as config_yml:
        config = yaml.safe_load(config_yml)
    # Date or prefix "this_thesis"
    date = "presentation" #"this_thesis" #datetime.date.today().strftime("%Y%m%d");
    folder = "/Users/alicamarieenderich/{0}_{1}_invariant_mass_plots/".format(date, period);  
    os.makedirs(folder, exist_ok=True);
    # input files from fitting process
    filename_data0 = "/Users/alicamarieenderich/this_thesis_LHC22f_invariant_mass_plots/this_thesis_LHC22f_pi0_data_ptspectrum_pp_13.6TeV_LHC22f_AnyTrack.root"
    filename_data1 = "/Users/alicamarieenderich/this_thesis_LHC23zc_invariant_mass_plots/this_thesis_LHC23zc_pi0_data_ptspectrum_pp_13.6TeV_LHC22f_AnyTrack.root"
    filename_data2 = "/Users/alicamarieenderich/this_thesis_LHC22o_minBias_invariant_mass_plots/this_thesis_LHC22o_minBias_pi0_data_ptspectrum_pp_13.6TeV_LHC22f_AnyTrack.root"
    # filename_mc = "/Users/alicamarieenderich/this_thesis_LHC23d1k_invariant_mass_plots/this_thesis_LHC23d1k_pi0_mc_ptspectrum_pp_13.6TeV_LHC22f_AnyTrack.root"
    print(config_file)
    print(period)
    print(filename)
    run(filename_data0, filename_data1, filename_data2,config,type,folder, period_str)
    print("input period_str", period_str)