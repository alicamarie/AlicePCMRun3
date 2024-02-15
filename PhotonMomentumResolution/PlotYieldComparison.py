# This code is based on the code provided by Daiki Sekihata in his repository
# https://github.com/dsekihat/photon_sw_run3 for photon analysis in ALICE Run 3


import numpy as np
import yaml
import datetime
import ROOT
from ROOT import TFile, TDirectory, THashList, TH1F, TH2F, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TGaxis, TPad, TLatex, TBox, TString, TLine
from ROOT import gStyle, gROOT, gSystem
from ROOT import TFile, TDirectory, THashList, TF1, TH1F, TCanvas, TH1
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kFullCircle, kOrange
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

def draw_comparison_yields(list_data, list_mc, index, cutname, date):
        # ymin = [6*1e-3, 125*1e-3,  10*1e-4]
        # ymax = [22*1e-3, 138*1e-3 , 110*1e-4]
        # for i in range(len(list_data)):
        #     print("MAX", list_data[i].GetMaximum(),list_mc[i].GetMaximum())
    #canvas plotting
        c1 = TCanvas("c0","c0",0,0,800,800);
        c1.Divide(1,2,1e-3,1e-3);
        p1 = c1.cd(1);
        p1.SetPad(0,0.3,1,1);
        p1.SetMargin(0.15,0.02,0.,0.15);
        p1.SetTicks(1,1);
        p1.SetLogy()

        ymin, ymax = 1e-9, 4*1e-3

        frame1 = p1.DrawFrame(0, ymin, 10., ymax);
        frame1.GetXaxis().SetTitle("#it{p}_{T} (GeV/c)");
        frame1.GetYaxis().SetTitle("fitparameter value");
        frame1.GetXaxis().SetTitleSize(0.045);
        frame1.GetYaxis().SetTitleSize(0.045);
        frame1.GetXaxis().SetTitleOffset(1.0);
        frame1.GetYaxis().SetTitleOffset(1.52);
        frame1.GetXaxis().SetLabelSize(0.045);
        frame1.GetYaxis().SetLabelSize(0.045);
        frame1.GetYaxis().SetMaxDigits(3);
        frame1.GetXaxis().SetLabelOffset(0.01);
        frame1.GetYaxis().SetLabelOffset(0.01);
        ROOT.SetOwnership(frame1,False);

        color = [kRed+1, kBlue+1, kGreen+2, kMagenta+2, kCyan+1]
        for i in range(len(list_data)):

                make_common_style(list_data[i], 20, 0.9, color[i], 1, 0);
                make_common_style(list_mc[i], 25, 0.9, color[i], 1, 0);
                #mc_list[i].SetLineStyle(2);
                list_data[i].Draw("E0same");
                list_mc[i].Draw("E0same");

        txt = TPaveText(0.,0.9,1.,0.95,"NDC");
        txt.SetFillColor(kWhite);
        txt.SetFillStyle(0);
        txt.SetBorderSize(0);
        txt.SetTextAlign(22);#middle,left
        txt.SetTextFont(42);#helvetica
        txt.SetTextSize(0.045);
        txt.AddText("{0} for different cuts in data and MC".format(cutname));
        txt.Draw();
        ROOT.SetOwnership(txt,False);

        txt = TPaveText(0.65,0.77,0.95,0.87,"NDC");
        txt.SetFillColor(kWhite);
        txt.SetFillStyle(0);
        txt.SetBorderSize(0);
        txt.SetTextAlign(32);#middle,left
        txt.SetTextFont(42);#helvetica
        txt.SetTextSize(0.02);
        txt.AddText("ALICE this thesis");
        txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
        txt.Draw();
        ROOT.SetOwnership(txt,False);

        leg = TLegend(0.65,0.62,0.95,0.77);
        leg.SetBorderSize(0);
        leg.SetFillColor(kWhite);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.03);
        leg.SetTextAlign(32);
        leg.SetTextFont(42);#helvetica
        cut = ['analysis_analysis', 'qc_qc', 'qc_ITSTPC_qc_ITSTPC', 'qc_TPConly_qc_TPConly']
        for i in range(len(list_data)):
            leg.AddEntry(list_data[i], "{0}".format(cut[i]),"LP");
        leg.Draw("");
        ROOT.SetOwnership(leg,False);

        leg = TLegend(0.15,0.72,0.40,0.82);
        leg.SetBorderSize(0);
        leg.SetFillColor(kWhite);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.03);
        leg.SetTextAlign(12);
        leg.SetTextFont(42);#helvetica
        leg.AddEntry(list_data[0], "Data #gamma candidates (LHC22f)","LP");
        leg.AddEntry(list_mc[0] , "M.C. rec. #gamma (LHC23d1k)","LP");
        leg.Draw("");
        ROOT.SetOwnership(leg,False);

        p2 = c1.cd(2);
        p2.SetPad(0,0,1,0.3);
        p2.SetMargin(0.15,0.02,0.22,0.0);
        p2.SetTicks(1,1);
        cut_in_ratio = 0.5


        frame2 = p2.DrawFrame(0,0.2,10.,1.6);
        frame2.GetXaxis().SetTitle("#it{p}_{T} (GeV/c)");
        frame2.GetYaxis().SetTitle("#frac{Data}{M.C.}");
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
        for i in range(len(list_data)):
                h1ratio = list_data[i].Clone("h1ratio");
                h1ratio.Reset();
                h1ratio.Sumw2();
                h1ratio.Divide(list_data[i], list_mc[i], 1., 1., "G");
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
        filepath = os.path.join(folder, "{0}_Parameter_Comparison_{1}.pdf".format(date, cutname));    
        c1.SaveAs(filepath);


def run(filename_data, filename_mc,config,type,folder):
    print(sys._getframe().f_code.co_name);
    arr_pt = np.array(config["common"]["pt_bin"],dtype=float);
    print("pT binning = ",arr_pt);
    print("type = ",type);
    print("input = ",filename_data);
    rootfile_data = TFile.Open(filename_data,"READ");
    rootfile_mc = TFile.Open(filename_mc,"READ");
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
    print(nsys); 

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
                list_ss   = rootfile.Get(ssname);
                for ic in range(0,nc):
                    cutname = cutnames[ic];
                    list_ss_cut = list_ss.FindObject(cutname);
                    print("cutname ", cutname)             
                    for ifit in range(0, nfit): 
                        fitname = "gausexplinear";
                        list_fitname  = list_ss_cut.FindObject(fitname);
                        print("fitname ", fitname)

                        list_plot = list_fitname.FindObject("fit_0.04_0.20_GeVc2")
                        yield_list.append(list_plot.FindObject("h1yield_param"))

                        parameter_list = [];
                        it = ROOT.TIter(list_plot)
                        obj = it.Next()
                        while obj:
                            objName = obj.GetName()
                            if "param" in objName:
                                parameter_list.append(obj)
                            obj = it.Next()
                        # print("HERE")
                        # for i in range(len(parameter_list)):
                        #     print(parameter_list[i].GetTitle())
        # pdf output of mass, amplitude and width for each fit and comparison of all cuts
                            parameter_comparison = [];
                            for i_parameter in range(5):
                                parameter_comparison.append(parameter_list[i_parameter]);
                            list_parameters_comparison.append(parameter_comparison);
            return list_parameters_comparison, yield_list  
    def loop_mc(rootfile): 
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
                list_ss   = rootfile.Get(ssname);
                for ic in range(0,nc):
                    cutname = cutnames[ic];
                    list_ss_cut = list_ss.FindObject(cutname);
                    print("cutname ", cutname)             
                    for ifit in range(0, nfit): 
                        fitname = "gausexp";
                        list_fitname  = list_ss_cut.FindObject(fitname);
                        print("fitname ", fitname)

                        list_plot = list_fitname.FindObject("fit_0.04_0.20_GeVc2")
                        yield_list.append(list_plot.FindObject("h1yield_param"))


                        parameter_list = [];
                        same_list = [];
                        mixed_list = [];
                        it = ROOT.TIter(list_plot)
                        obj = it.Next()
                        while obj:
                            objName = obj.GetName()
                            if "param" in objName:
                                parameter_list.append(obj)
                            obj = it.Next()
                        # print("HERE")
                        # for i in range(len(parameter_list)):
                        #     print(parameter_list[i].GetTitle())
        # pdf output of mass, amplitude and width for each fit and comparison of all cuts
                            parameter_comparison = [];
                            for i_parameter in range(5):
                                parameter_comparison.append(parameter_list[i_parameter]);
                            list_parameters_comparison.append(parameter_comparison);
            return list_parameters_comparison, yield_list     
    list_parameters_data, yield_list_data = loop_data(rootfile_data)
    list_parameters_mc, yield_list_mc = loop_mc(rootfile_mc)


    # print("START")
    # for i in range(len(list_parameters_data)):
    #     for j in range(len(list_parameters_data[1])):
    #         print(list_parameters_data[i][j].GetTitle())
    

    # list_data_mean = []
    # list_data_lambda = []
    # list_data_sigma = []
    # for i in range(len(list_parameters_data)):
    #     for j in range(len(list_parameters_data[1])):
    #         if list_parameters_data[i][j].GetTitle() == "mean":
    #             list_data_mean.append(list_parameters_data[i][j])
    #         if list_parameters_data[i][j].GetTitle() == "lambda":
    #             list_data_lambda.append(list_parameters_data[i][j])
    #         if list_parameters_data[i][j].GetTitle() == "sigma":
    #             list_data_sigma.append(list_parameters_data[i][j])


    # list_mc_mean = []
    # list_mc_lambda = []
    # list_mc_sigma = []
    # for i in range(len(list_parameters_mc)):
    #     for j in range(len(list_parameters_mc[1])):
    #         if list_parameters_mc[i][j].GetTitle() == "mean":
    #             list_mc_mean.append(list_parameters_mc[i][j])
    #         if list_parameters_mc[i][j].GetTitle() == "lambda":
    #             list_mc_lambda.append(list_parameters_mc[i][j])
    #         if list_parameters_mc[i][j].GetTitle() == "sigma":
    #             list_mc_sigma.append(list_parameters_mc[i][j])
    
    # print(list_data_lambda)
    # print(list_data_mean)
    # print(list_data_sigma)
    # print(list_mc_lambda)
    # print(list_mc_mean)
    # print(list_mc_sigma)

    # list_data = [list_data_lambda, list_data_mean, list_data_sigma]
    # list_mc = [list_mc_lambda, list_mc_mean, list_mc_sigma]

    print(yield_list_data)
    print(yield_list_mc)
    # outname_histo_linear = os.path.join(folder, "{0}_{1}_InvMass_Parameters_Comparison.pdf".format(date, period));
    for i in range(3):
        draw_comparison_yields(yield_list_data, yield_list_mc, i,yield_list_data[0].GetTitle(), date)

period_array = ["LHC22f", "LHC23d1k"]

filename_array = ["/Users/alicamarieenderich/AnalysisResults/AnalysisResults_124838_LHC22f_pass4.root",
            "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_124837_LHC23d1k.root",]

type_array = ["data", "mc"]

config_array = ["/Users/alicamarieenderich/202312_invariant_mass_code/configs/config_pp_13.6TeV_pi0_LHC22f.yml",
                "/Users/alicamarieenderich/202312_invariant_mass_code/configs/config_pp_13.6TeV_pi0_LHC23d1k.yml"]

for i in range(len(period_array)):
    period = period_array[i]
    filename = filename_array[i]
    cutname = "qc"
    suffix = "AnyTrack";
    type = type_array[i]
    config_file = config_array[i]
    with open(config_file, "r", encoding="utf-8") as config_yml:
        config = yaml.safe_load(config_yml)
    date = "this_thesis" #datetime.date.today().strftime("%Y%m%d");
    folder = "/Users/alicamarieenderich/{0}_{1}_invariant_mass_plots/".format(date, period);  
    os.makedirs(folder, exist_ok=True);
    filename_data = "/Users/alicamarieenderich/this_thesis_LHC22f_invariant_mass_plots/this_thesis_LHC22f_pi0_data_ptspectrum_pp_13.6TeV_LHC22f_AnyTrack.root"
    filename_mc = "/Users/alicamarieenderich/this_thesis_LHC23d1k_invariant_mass_plots/this_thesis_LHC23d1k_pi0_mc_ptspectrum_pp_13.6TeV_LHC22f_AnyTrack.root"

    run(filename_data, filename_mc,config,type,folder)