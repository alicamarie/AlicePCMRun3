# This code is based on the code provided by Daiki Sekihata in his repository
# https://github.com/dsekihat/photon_sw_run3 for photon analysis in ALICE Run 3


import os, sys, shutil
from array import array
import math
import argparse
import numpy as np
from ctypes import *
import datetime
import yaml
import ROOT
ROOT.gROOT.SetBatch(True);
from ROOT import TFile, THashList, TF1
from ROOT import TFile, TDirectory, THashList, TH1F, TH2F, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TLine, TPython
from ROOT import gStyle, gROOT, gSystem
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange
gStyle.SetOptStat(0);
gStyle.SetOptTitle(0);
from FomattingMaterialBudget import FrameSettings, FrameSettingsRatio, RatioLegendSettings, ALICEtext

#_____________________________________________________________________
def make_common_style(g1,marker,size,color,width=1,fill=0):
    g1.SetMarkerStyle(marker);
    g1.SetMarkerColor(color)
    g1.SetMarkerSize(size);
    g1.SetLineColorAlpha(color, 0.8)
    g1.SetLineWidth(width);
    g1.SetFillColor(color)
    g1.SetFillStyle(fill);

class calc_ratio:
    def __init__(self):
        print("default constructor is called");
    def __init__(self, config, filename_data, filename_mc, cutname, folder, period_data, period_mc, suffix):
        #print("filename_data = {0} , filename_mc = {1} , cutname = {2}".format(filename_data, filename_mc, cutname));
        self.rootfile_data = TFile.Open(filename_data, "READ");
        self.rootfile_mc   = TFile.Open(filename_mc  , "READ");
        self.suffix = suffix;
        self.folder = folder;
        self.cutname = cutname
        self.config = config

    def __del__(self):
        if self.rootfile_data.IsOpen():
            print("close input data root file.");
            self.rootfile_data.Close();
        if self.rootfile_mc.IsOpen():
            print("close input mc root file.");
            self.rootfile_mc.Close();

    def draw_ratio_weights_combined(self, filename_isospin, filename_wire, arr_pt_cuts, date):
        rootfile_isospin = TFile.Open(filename_isospin, "READ");
        rootfile_wire = TFile.Open(filename_wire, "READ");
        arr_pt_cuts = self.config["common"]["pt_cuts"];

        list_isospin = rootfile_isospin.Get("qc");
        list_wire = rootfile_wire.Get("qc")
        list_isospin.ls()
        list_wire.ls() 

        r_bins = [0, 14, 30, 42, 58, 69, 90]
        
        data_list_isospin = []
        mc_list_isospin = []            
        for ipt in range(len(arr_pt_cuts)):
            h1data_complete = list_isospin.FindObject("h1NGamma_data_pt{0}".format(ipt));
            h1mc_complete = list_isospin.FindObject("h1NGamma_mc_pt{0}".format(ipt));
            h1data_complete.SetDirectory(0);
            h1mc_complete.SetDirectory(0);
    
        #normalization
            h1data_complete.Sumw2()
            h1mc_complete.Sumw2()  

            data_list_isospin.append(h1data_complete)
            mc_list_isospin.append(h1mc_complete)

        data_list_wire = []
        mc_list_wire = []     
        for ipt in range(len(arr_pt_cuts)):
            h1data_complete = list_wire.FindObject("h1NGamma_data_pt{0}".format(ipt));
            h1mc_complete = list_wire.FindObject("h1NGamma_mc_pt{0}".format(ipt));
            h1data_complete.SetDirectory(0);
            h1mc_complete.SetDirectory(0);

        #normalization
            h1data_complete.Sumw2()
            h1mc_complete.Sumw2()  

            data_list_wire.append(h1data_complete)
            mc_list_wire.append(h1mc_complete)

        h1ratio_isospin = []
        for i in range(len(data_list_isospin)):
            h1ratio = mc_list_isospin[i].Clone("h1ratio");
            h1ratio.Reset();
            h1ratio.Sumw2();
            h1ratio.Divide(data_list_isospin[i], mc_list_isospin[i], 1., 1., "G");
            # h1ratio.Draw("E0same");
            h1ratio.SetDirectory(0);
            markerstlye_data = [21, 20, 22]
            markerstlye_mc = [25, 24, 26 ]
            size = [2.0,1.8, 1.6]
            color = [kRed+1,kGreen+2, kBlue+1, kCyan+2, kMagenta+2]
            make_common_style(h1ratio, markerstlye_data[i], size[i], color[i], 1, 0)
            h1ratio_isospin.append(h1ratio)
            ROOT.SetOwnership(h1ratio,False);

        h1ratio_wire = []
        for i in range(len(data_list_wire)):
            h1ratio = mc_list_wire[i].Clone("h1ratio");
            h1ratio.Reset();
            h1ratio.Sumw2();
            h1ratio.Divide(data_list_wire[i], mc_list_wire[i], 1., 1., "G");
            # h1ratio.Draw("E0same");
            h1ratio.SetDirectory(0);
            markerstlye_data = [21, 20, 22]
            markerstlye_mc = [25, 24, 26 ]
            size = [2.0,1.8, 1.6]
            color = [kRed+1,kGreen+2, kBlue+1, kCyan+2, kMagenta+2]
            make_common_style(h1ratio, markerstlye_mc[i], size[i], color[i], 1, 0)
            h1ratio_wire.append(h1ratio)
            ROOT.SetOwnership(h1ratio,False);
                
    #canvas plotting
        c1 = TCanvas("c0","c0",0,0,800,800);
        c1.Divide(1,2,1e-3,1e-3);
        p1 = c1.cd(1);
        p1.SetPad(0,0.3,1,1);
        p1.SetMargin(0.15,0.02,0.,0.15);
        p1.SetTicks(1,1);
        # p1.SetLogy();

        frame1 = p1.DrawFrame(0., 0, 90., 4)#(0., ymin, 90, ymax); #(0., 1e-20, 10., 1e-1);#
        frame1.GetXaxis().SetTitle("#it{r}_{xy} (cm)")
        frame1.GetYaxis().SetTitle("Calibration weights");
        FrameSettings(frame1)

        color = [kRed+1,kGreen+2, kBlue+1, kCyan+2, kMagenta+2]
        for i in range(len(h1ratio_isospin)):

            markerstlye_data = [21, 20, 22]
            markerstlye_mc = [25, 24, 26 ]
            size = [2.0,1.8, 1.6]
            make_common_style(h1ratio_isospin[i], markerstlye_data[i], size[i], color[i], 1, 0);
            make_common_style(h1ratio_wire[i]  , markerstlye_mc[i], size[i], color[i], 1, 0);
            h1ratio_isospin[i].Draw("E0same");
            h1ratio_wire[i].Draw("E0same");

        txt = TPaveText(0.,0.9,1.,0.95,"NDC");
        txt.SetFillColor(kWhite);
        txt.SetFillStyle(0);
        txt.SetBorderSize(0);
        txt.SetTextAlign(22);#middle,left
        txt.SetTextFont(42);#helvetica
        txt.SetTextSize(0.045);
        txt.AddText("Ratio of the calibration weights #it{#Omega}_{#it{i}} / #it{#omega}_{#it{i}} for different cuts in #it{p}_{T}");
        txt.Draw();
        ROOT.SetOwnership(txt,False);

        ALICEtext("thesis")
        
        leg = TLegend(0.25,0.72,0.40,0.82);
        leg.SetBorderSize(0);
        leg.SetFillColor(kWhite);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.03);
        leg.SetTextAlign(12);
        leg.SetTextFont(42);#helvetica
        leg.AddEntry(h1ratio_isospin[0], "Isospin calibration weights #it{#Omega}_{#it{i}}","LP");
        leg.AddEntry(h1ratio_wire[0] , "Wire calibration weights #it{#omega}_{#it{i}}","LP");
        leg.Draw("");
        ROOT.SetOwnership(leg,False);

        p2 = c1.cd(2);
        p2.SetPad(0,0,1,0.3);
        p2.SetMargin(0.15,0.02,0.22,0.0);
        p2.SetTicks(1,1);
        cut_in_ratio = 0.5


        frame2 = p2.DrawFrame(0.,0.2,90.,1.2);
        frame2.GetXaxis().SetTitle("#it{R}_{xy} (cm)");
        frame2.GetYaxis().SetTitle("#it{#Omega}_{#it{i}} / #it{#omega}_{#it{i}}");
        FrameSettingsRatio(frame2)

        minimum = []
        maximum = []
        for i in range(len(h1ratio_isospin)):
            h1ratio = h1ratio_isospin[i].Clone("h1ratio");
            h1ratio.Reset();
            h1ratio.Sumw2();
            h1ratio.Divide(h1ratio_isospin[i], h1ratio_wire[i], 1., 1., "G");
            h1ratio.Draw("E0same");
            h1ratio.SetDirectory(0);
            ROOT.SetOwnership(h1ratio,False);
            minimum.append(h1ratio.GetMinimum())
            maximum.append(h1ratio.GetMaximum())

        line1 = TLine(0.,1,90,1);
        line1.SetLineColor(kBlack);
        line1.SetLineStyle(1);
        line1.SetLineWidth(2);
        line1.Draw("");
        ROOT.SetOwnership(line1,False);

        txt = TPaveText(0.,0.9,1.,0.95,"NDC");
        txt.SetFillColor(kWhite);
        txt.SetFillStyle(0);
        txt.SetBorderSize(0);
        txt.SetTextAlign(22);#middle,left
        txt.SetTextFont(42);#helvetica
        txt.SetTextSize(0.075);
        txt.AddText("Ratio = {}".format(round(h1ratio.GetMaximum(),4)));
        txt.Draw();
        ROOT.SetOwnership(txt,False);

        c1.Modified();
        c1.Update();
        ROOT.SetOwnership(c1,False);
        filepath = os.path.join(self.folder, "{0}_NGamma_pt_Ratio_Combined_{2}.pdf".format(date, ipt, self.suffix));    
        c1.SaveAs(filepath);

# #________________________________________________
# if __name__:
#     cutname = "qc"
#     period_mc = "LHC23d1k";
#     period_data = "LHC22f"
#     suffix = "AnyTrack";
#     filename_data = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_147811_LHC22_pass4_lowIR.root"
#     filename_mc = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_147812_LHC23d1k.root"
#     config_file = "CalibrationWeights/config_pp_13.6TeV_LHC22f_material.yml"
#     with open(config_file, "r", encoding="utf-8") as config_yml:
#         config = yaml.safe_load(config_yml)
#     date = "this_thesis" # datetime.date.today().strftime("%Y%m%d");
#     folder = "/Users/alicamarieenderich/{0}_calibration_weights/".format(date);  
#     os.makedirs(folder, exist_ok=True);

#     for type in ["data", "mc"]:
#         if type == "data": 
#             file = filename_data;
#         elif type == "mc":
#             file = filename_mc;
#     NGamma = calc_ratio(config, filename_data, filename_mc, cutname, folder, period_data, period_mc, suffix);
#     arr_pt_cuts = config["common"]["pt_cuts"];

#     file_isospin = "/Users/alicamarieenderich/this_thesis_calibration_weights/this_thesis_calibration_weights_pp_13.6TeV_LHC22f_and_LHC23d1k_AnyTrack_new.root"
#     file_wire = "/Users/alicamarieenderich/this_thesis_calibration_weights/this_thesis_calibration_weights_wire_pp_13.6TeV_LHC22f_and_LHC23d1k_AnyTrack_new.root"

#     NGamma.draw_ratio_weights_combined(file_isospin, file_wire, arr_pt_cuts, date)