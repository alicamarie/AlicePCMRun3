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

class calc_Omega_isospin:
    def __init__(self):
        print("default constructor is called");
    def __init__(self, config, filename_data, filename_mc, cutname, folder, period_data, period_mc, suffix):
        self.config = config
        self.cutname = cutname;
        self.folder = folder;
        self.suffix = suffix;

        self.rootfile_data = TFile.Open(filename_data, "READ");
        self.rootfile_mc   = TFile.Open(filename_mc  , "READ");

        self.rootdir_mc_gen  = self.rootfile_mc.Get("material-budget-mc")
        self.list_gen        = self.rootdir_mc_gen.Get("Generated");
        self.list_ev_gen     = self.rootdir_mc_gen.Get("Event");
        self.list_ev_mc_gen  = self.list_ev_gen.FindObject("PCMDalitzEE");
        self.rootdir_mc_rec  = self.rootfile_mc.Get("material-budget-mc");
        self.list_v0_mc_rec  = self.rootdir_mc_rec.Get("V0");
        self.list_ev_rec     = self.rootdir_mc_rec.Get("Event")
        self.list_ev_mc_rec  = self.list_ev_rec.FindObject("PCMDalitzEE");
        self.list_cut_mc_rec = self.list_v0_mc_rec.FindObject(cutname);

        self.rootdir_mc_pcmqc    = self.rootfile_mc.Get("pcm-qc-mc");
        self.list_gen_pcmqc        = self.rootdir_mc_pcmqc.Get("Generated");
        self.list_ev_mc_pcm  = self.rootdir_mc_pcmqc.Get("Event");
        self.h1nch_mc_gen    = self.list_ev_mc_pcm.FindObject("hMultNTracksPV").Clone("h1mult");
        self.nev_gen         = self.h1nch_mc_gen.GetEntries();
        self.nch_gen         = self.h1nch_mc_gen.GetMean();

        self.h1nch_mc_rec    = self.list_ev_mc_pcm.FindObject("hMultNTracksPV");
        self.nch_rec         = self.h1nch_mc_rec.GetMean();
        self.nev_rec         = self.h1nch_mc_rec.GetEntries();
    
        self.rootdir_data    = self.rootfile_data.Get("material-budget");
        self.list_v0_data    = self.rootdir_data.Get("V0");
        self.list_ev_data_1  = self.rootdir_data.Get("Event");
        self.list_ev_data    = self.list_ev_data_1.FindObject("PCMDalitzEE");
        self.list_cut_data   = self.list_v0_data.FindObject(cutname);
        self.rootdir_data_pcmqc    = self.rootfile_data.Get("pcm-qc");
        self.list_ev_data_pcm  = self.rootdir_data_pcmqc.Get("Event");
        self.h1nch_data      = self.list_ev_data_pcm.FindObject("hMultNTracksPV");
        self.nev_data        = self.h1nch_data.GetEntries();
        self.nch_data        = self.h1nch_data.GetMean();
        
        self.period_data = period_data;
        self.period_mc = period_mc
        self.arr_rxy = np.array([0,1,2,3,4,5], dtype=float);

    def __del__(self):
        if self.rootfile_data.IsOpen():
            print("close input data root file.");
            self.rootfile_data.Close();
        if self.rootfile_mc.IsOpen():
            print("close input mc root file.");
            self.rootfile_mc.Close();

    #_____________________________________________________________________
    def calc_in_Rbins(self, arr_rxy, arr_pt_cuts):

        outlist = THashList();
        outlist.SetOwner(True);
        outlist.SetName("outlist");
        outlist.Add(self.h1nch_mc_gen);
        outlist.Add(self.h1nch_data)

        hs_mc = self.list_cut_mc_rec.FindObject("hs_conv_point").Clone("hs_rec");
        outlist.Add(hs_mc)
        h1_mc = hs_mc.Projection(2 ,"");
        h1_mc.SetDirectory(0);
        ROOT.SetOwnership(h1_mc, False);
        h1_mc.Sumw2();
        #h1_mc_rec.RebinX(2);
        h1_mc.Scale(1,"width");
        h1_mc.Scale(1/self.nev_rec);
        h1_mc.Scale(1/self.nch_rec);#nch

        hs_data = self.list_cut_data.FindObject("hs_conv_point").Clone("hs_data");
        outlist.Add(hs_data)
        h1_data = hs_data.Projection(2 ,"");
        h1_data.SetDirectory(0);
        ROOT.SetOwnership(h1_data, False);
        h1_data.Sumw2();
        h1_data.Scale(1,"width");
        h1_data.Scale(1/self.nev_data);
        h1_data.Scale(1/self.nch_data);

        self.arr_rxy = arr_rxy
        h2_mc = hs_mc.Projection(0,1, "")
        h2_mc.SetName("h2_rec")
        outlist.Add(h2_mc)
        h2_data = hs_data.Projection(0,1, "")
#        h2_data.Scale(1,"width")
        h2_data.SetName("h2_data")
        h2_data.Sumw2();
        outlist.Add(h2_data)
   
        for ipt in range(0, len(arr_pt_cuts)):
            pt_cut = arr_pt_cuts[ipt]
            bin_pt_end_mc = h2_mc.GetYaxis().GetNbins()
            bin_pt_cut_mc = h2_mc.GetYaxis().FindBin(pt_cut + 1e-6);
            print("pt:", pt_cut)
            rxy = [0., 14., 30., 42., 58., 69., 90.]
            h1NGamma_mc = TH1F("h1NGamma", "NGamma", len(rxy)-1, array('d', rxy))
            h1NGamma_mc.SetName("h1NGamma_mc_pt{0}".format(ipt))
            h1NGamma_mc.SetTitle("NGamma/Nch for #it{{p}}_{{T}} >= {0} GeV/c, mc".format(pt_cut))
            h1NGamma_mc.SetXTitle("#it{r}_{xy} (cm)")
            h1NGamma_mc.SetYTitle("#it{N}_{#gamma}/#it{N}_{ch}")

            for ir in range(0, len(self.arr_rxy)-1):
                r1 = arr_rxy[ir];
                r2 = arr_rxy[ir+1];
                dr = r2 - r1;
                bin_r1_mc = h2_mc.GetXaxis().FindBin(r1 + 1e-6);
                bin_r2_mc = h2_mc.GetXaxis().FindBin(r2 - 1e-6);
                ngamma_err = c_double(0.0)
                ngamma = h2_mc.IntegralAndError(bin_r1_mc, bin_r2_mc, bin_pt_cut_mc, bin_pt_end_mc, ngamma_err, "")#"width")
                print("isospin method", r1, "ngamma mc", ngamma, "+/-", ngamma_err.value)
                h1NGamma_mc.SetBinContent(ir+1, ngamma/self.nch_rec/self.nev_rec)
                h1NGamma_mc.SetBinError(ir+1, ngamma_err.value/self.nch_rec/self.nev_rec)
            outlist.Add(h1NGamma_mc)

            pt_cut = arr_pt_cuts[ipt]
            bin_pt_end_data = h2_data.GetYaxis().GetNbins()
            bin_pt_cut_data = h2_data.GetYaxis().FindBin(pt_cut);

            rxy = [0., 14., 30., 42., 58., 69., 90.]
            h1NGamma_data = TH1F("h1NGamma", "NGamma", len(rxy)-1, array('d', rxy))
            h1NGamma_data.SetName("h1NGamma_data_pt{0}".format(ipt))
            h1NGamma_data.SetTitle("NGamma/Nch for #it{{p}}_{{T}} >= {0} GeV/c, data".format(pt_cut))
            h1NGamma_data.SetXTitle("#it{r}_{xy} (cm)")
            h1NGamma_data.SetYTitle("#it{N}_{#gamma}/#it{N}_{ch}")

            for ir in range(0, len(self.arr_rxy)-1):
                r1 = arr_rxy[ir];
                r2 = arr_rxy[ir+1];
                dr = r2 - r1;
                bin_r1_data = h2_data.GetXaxis().FindBin(r1 + 1e-6);
                bin_r2_data = h2_data.GetXaxis().FindBin(r2 - 1e-6);
                ngamma_err = c_double(0.0)
                ngamma = h2_data.IntegralAndError(bin_r1_data, bin_r2_data, bin_pt_cut_data, bin_pt_end_data, ngamma_err, "")
                print("isospin method", r1, "ngamma data", ngamma, "+/-", ngamma_err.value)
                h1NGamma_data.SetBinContent(ir+1, ngamma/self.nch_data/self.nev_data)
                h1NGamma_data.SetBinError(ir+1, ngamma_err.value/self.nch_data/self.nev_data)

            outlist.Add(h1NGamma_data)

            h1Omega = h1NGamma_data.Clone("h1Omega_isopsin_pt{0}".format(pt_cut))
            h1Omega.Divide(h1NGamma_data, h1NGamma_mc, 1., 1., "G")
            h1Omega.SetTitle("#it{{#Omega}} for #it{{p}}_{{T}} >= {0} GeV/c, data".format(pt_cut))
            h1Omega.SetXTitle("#it{r}_{xy} (cm)")
            h1Omega.SetYTitle("#it{#Omega}")
            outlist.Add(h1Omega)   
        return outlist
    
    #_____________________________________________________________________
    def run(self, filename, type, date):
        print(sys._getframe().f_code.co_name);
        print("reading...", filename);

        taskname = "material-budget";
        if type == "mc":
            taskname = "material-budget-mc";

        rootfile = TFile.Open(filename, "READ");
        rootdire = rootfile.Get(taskname);
        list_v0 = rootdire.Get("V0");
        list_ev_1     = rootdire.Get("Event");
        list_ev = list_ev_1.FindObject("PCMPCM");

        outname = os.path.join(self.folder, "{0}_calibration_weights_{1}_{2}TeV_{3}_and_{4}_{5}_new.root".format(date, self.config["common"]["system"], self.config["common"]["energy"], self.config["common"]["period"], self.config["common"]["period_mc"] ,self.suffix));
        print("out file name = ", outname);
        outfile = TFile(outname, "RECREATE");

        arr_rxy = self.config["common"]["rxy_bin"];
        arr_eta = self.config["common"]["eta_bin"];
        arr_pt_cuts = self.config["common"]["pt_cuts"];
        cutnames =self. config[type]["subsystems"][0]['cutnames']
        print("cutnames", cutnames); 
        nc = len(cutnames);
        for ic in range(0,nc):
            cutname = cutnames[ic];
            print("analyzing...", cutname);
            list_v0_cut = list_v0.FindObject(cutname);
            outlist = self.calc_in_Rbins(arr_rxy, arr_pt_cuts);
            outlist.SetName(cutname);
            outlist.SetOwner(True);
            outfile.WriteTObject(outlist);
            outlist.Clear();
        outfile.Close();

        rootfile.Close();

    def draw_NGamma_isospin(self, filename, ipt, date):
        rootfile= TFile.Open(filename, "READ");

        arr_pt_cuts = self.config["common"]["pt_cuts"];

        list= rootfile.Get("qc");

        r_bins = [0, 14, 30, 42, 58, 69, 90]
        
        h1data_complete = list.FindObject("h1NGamma_data_pt{0}".format(ipt));
        h1mc_complete = list.FindObject("h1NGamma_mc_pt{0}".format(ipt));
        h1data_complete.SetDirectory(0);
        h1mc_complete.SetDirectory(0);
    
    #normalization
        h1data_complete.Sumw2()
        h1mc_complete.Sumw2()  

    #style   
        make_common_style(h1data_complete, 20, 1.0, kBlue+1, 1, 0);
        make_common_style(h1mc_complete  , 20, 1.0, kRed+1, 1, 0);
        ROOT.SetOwnership(h1data_complete, False);
        ROOT.SetOwnership(h1mc_complete, False);

        ymax = max(h1data_complete.GetMaximum() , h1mc_complete.GetMaximum()) * 3.5;
        ymin = max(h1data_complete.GetMinimum() , h1mc_complete.GetMinimum()) * -0.3;
        if ymin ==0:
            ymin = 1e-6

    #canvas plotting
        c1 = TCanvas("c0","c0",0,0,800,800);
        c1.Divide(1,2,1e-3,1e-3);
        p1 = c1.cd(1);
        p1.SetPad(0,0.3,1,1);
        p1.SetMargin(0.15,0.02,0.,0.15);
        p1.SetTicks(1,1);
       #p1.SetLogy();

        frame1 = p1.DrawFrame(0., ymin, 90, ymax); #(0., 1e-20, 10., 1e-1);#
        frame1.GetXaxis().SetTitle("#it{r}_{xy} (cm)")
        frame1.GetYaxis().SetTitle("#it{N}_{#gamma}/#it{N}_{ch}");
        FrameSettings(frame1)

        h1mc_complete.Draw("E0Hsame");
        h1data_complete.Draw("E0Hsame");

        txt = TPaveText(0.,0.9,1.,0.95,"NDC");
        txt.SetFillColor(kWhite);
        txt.SetFillStyle(0);
        txt.SetBorderSize(0);
        txt.SetTextAlign(22);#middle,left
        txt.SetTextFont(42);#helvetica
        txt.SetTextSize(0.045);
        txt.AddText("Calibration weight #it{{#Omega}}_{{#it{{i}}}} for #it{{p}}_{{T}} > {0} GeV/c".format(arr_pt_cuts[ipt]));
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
        leg.AddEntry(h1data_complete, "Data (LHC22f)","LP");
        leg.AddEntry(h1mc_complete  , "M.C. rec. (LHC23d1k)","LP");
        leg.Draw("");
        ROOT.SetOwnership(leg,False);

        p2 = c1.cd(2);
        p2.SetPad(0,0,1,0.3);
        p2.SetMargin(0.15,0.02,0.22,0.0);
        p2.SetTicks(1,1);
        cut_in_ratio = 0.5


        frame2 = p2.DrawFrame(0.,0.6,90.,1.2);
        frame2.GetXaxis().SetTitle("#it{r}_{xy} (cm)");
        frame2.GetYaxis().SetTitle("#it{#Omega}_{#it{i}}");
        FrameSettingsRatio(frame2)

        h1ratio = h1data_complete.Clone("h1ratio");
        h1ratio.Reset();
        h1ratio.Sumw2();
        h1ratio.Divide(h1data_complete, h1mc_complete, 1., 1., "G");
        h1ratio.Draw("E0same");
        h1ratio.SetDirectory(0);
        ROOT.SetOwnership(h1ratio,False);

        line1 = TLine(0.,1,90,1);
        line1.SetLineColor(kBlack);
        line1.SetLineStyle(1);
        line1.SetLineWidth(2);
        line1.Draw("");
        ROOT.SetOwnership(line1,False);

        c1.Modified();
        c1.Update();
        ROOT.SetOwnership(c1,False);
        filepath = os.path.join(self.folder, "{0}_NGamma_isospin_pt{1}_{2}.pdf".format(date, ipt, self.suffix));    
        c1.SaveAs(filepath);

    def draw_NGamma_isospin_combined(self, filename, arr_pt_cuts, date):
        rootfile= TFile.Open(filename, "READ");

        arr_pt_cuts = self.config["common"]["pt_cuts"];

        list= rootfile.Get("qc");

        r_bins = [0, 14, 30, 42, 58, 69, 90]
        
        data_list = []
        mc_list = []            
        for ipt in range(len(arr_pt_cuts)):
            h1data_complete = list.FindObject("h1NGamma_data_pt{0}".format(ipt));
            h1mc_complete = list.FindObject("h1NGamma_mc_pt{0}".format(ipt));
            h1data_complete.SetDirectory(0);
            h1mc_complete.SetDirectory(0);
    
        #normalization
            h1data_complete.Sumw2()
            h1mc_complete.Sumw2()  

        #style   
            markerstlye_data = [21, 20, 22]
            markerstlye_mc = [25, 24, 26 ]
            make_common_style(h1data_complete, markerstlye_data[ipt], 1.0, kBlue+1, 1, 0);
            make_common_style(h1mc_complete  , markerstlye_mc[ipt], 1.0, kRed+1, 1, 0);
            ROOT.SetOwnership(h1data_complete, False);
            ROOT.SetOwnership(h1mc_complete, False);

            data_list.append(h1data_complete)
            mc_list.append(h1mc_complete)

        ymax = max(h1data_complete.GetMaximum() , h1mc_complete.GetMaximum()) * 2.;
        ymin = max(h1data_complete.GetMinimum() , h1mc_complete.GetMinimum()) * -1.;

    #canvas plotting
        c1 = TCanvas("c0","c0",0,0,800,800);
        c1.Divide(1,2,1e-3,1e-3);
        p1 = c1.cd(1);
        p1.SetPad(0,0.3,1,1);
        p1.SetMargin(0.15,0.02,0.,0.15);
        p1.SetTicks(1,1);
        p1.SetLogy();

        frame1 = p1.DrawFrame(0., 1e-5, 90., 4*1e-2)#(0., ymin, 90, ymax); #(0., 1e-20, 10., 1e-1);#
        frame1.GetXaxis().SetTitle("#it{r}_{xy} (cm)")
        frame1.GetYaxis().SetTitle("#it{N}_{#gamma}/#it{N}_{ch}");
        FrameSettings(frame1)

        color = [kRed+1,kGreen+2, kBlue+1, kCyan+2, kMagenta+2]
        for i in range(len(data_list)):

            markerstlye_data = [21, 20, 22]
            markerstlye_mc = [25, 24, 26 ]
            size = [2.0,1.8, 1.6]
            make_common_style(data_list[i], markerstlye_data[i], size[i], color[i], 1, 0);
            make_common_style(mc_list[i]  , markerstlye_mc[i], size[i], color[i], 1, 0);
            #mc_list[i].SetLineStyle(2);
            data_list[i].Draw("E0same");
            mc_list[i].Draw("E0same");

        txt = TPaveText(0.,0.9,1.,0.95,"NDC");
        txt.SetFillColor(kWhite);
        txt.SetFillStyle(0);
        txt.SetBorderSize(0);
        txt.SetTextAlign(22);#middle,left
        txt.SetTextFont(42);#helvetica
        txt.SetTextSize(0.045);
        txt.AddText("Calibration weight #it{#Omega}_{#it{i}} for different cuts in #it{p}_{T}");
        txt.Draw();
        ROOT.SetOwnership(txt,False);

        ALICEtext("thesis")

        leg = TLegend(0.65,0.50,0.95,0.65);
        leg.SetBorderSize(0);
        leg.SetFillColor(kWhite);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.045);
        leg.SetTextAlign(32);
        leg.SetTextFont(42);#helvetica
        for i in range(len(data_list)):
            leg.AddEntry(data_list[i], "#it{{p}}_{{T}} > {0} GeV/c".format(arr_pt_cuts[i]),"LP");
            #leg.AddEntry(h1mc_complete  , "M.C. rec. #gamma (LHC23d1k)","LP");
        leg.Draw("");
        ROOT.SetOwnership(leg,False);

        leg = TLegend(0.25,0.72,0.40,0.82);
        leg.SetBorderSize(0);
        leg.SetFillColor(kWhite);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.03);
        leg.SetTextAlign(12);
        leg.SetTextFont(42);#helvetica
        leg.AddEntry(data_list[0], "Data (LHC22f)","LP");
        leg.AddEntry(mc_list[0] , "M.C. rec. (LHC23d1k)","LP");
        leg.Draw("");
        ROOT.SetOwnership(leg,False);

        p2 = c1.cd(2);
        p2.SetPad(0,0,1,0.3);
        p2.SetMargin(0.15,0.02,0.22,0.0);
        p2.SetTicks(1,1);
        cut_in_ratio = 0.5

        frame2 = p2.DrawFrame(0.,0.5,90.,1.5);
        frame2.GetXaxis().SetTitle("#it{R}_{xy} (cm)");
        frame2.GetYaxis().SetTitle("#it{#Omega}_{#it{i}}");
        FrameSettingsRatio(frame2)

        color = [kRed+1,kGreen+2, kBlue+1, kCyan+2, kMagenta+2]
        for i in range(len(data_list)):
            h1ratio = mc_list[i].Clone("h1ratio");
            h1ratio.Reset();
            h1ratio.Sumw2();
            h1ratio.Divide(data_list[i], mc_list[i], 1., 1., "G");
            h1ratio.Draw("E0same");
            h1ratio.SetDirectory(0);
            ROOT.SetOwnership(h1ratio,False);

        line1 = TLine(0.,1,90,1);
        line1.SetLineColor(kBlack);
        line1.SetLineStyle(1);
        line1.SetLineWidth(2);
        line1.Draw("");
        ROOT.SetOwnership(line1,False);

        c1.Modified();
        c1.Update();
        ROOT.SetOwnership(c1,False);
        filepath = os.path.join(self.folder, "{0}_NGamma_pt_isospin_combined_{2}.pdf".format(date, ipt, self.suffix));    
        c1.SaveAs(filepath);

# #________________________________________________
# if __name__:
#     cutname = "qc"
#     period_mc = "LHC23d1k";
#     period_data = "LHC22f"
#     suffix = "AnyTrack";
#     filename_data = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_147811_LHC22_pass4_lowIR.root"
#     filename_mc = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_147812_LHC23d1k.root"
#     config_file = "config_pp_13.6TeV_LHC22f_material.yml"
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
#     NGamma = calc_Omega_isospin(config, filename_data, filename_mc, cutname, folder, period_data, period_mc, suffix);
#     NGamma.run(file, type, date);

#     NGamma = calc_Omega_isospin(config, filename_data, filename_mc, cutname, folder, period_data, period_mc, suffix);
#     arr_pt_cuts = config["common"]["pt_cuts"];

#     file = "/Users/alicamarieenderich/this_thesis_calibration_weights/this_thesis_calibration_weights_pp_13.6TeV_LHC22f_and_LHC23d1k_AnyTrack_new.root"

#     for ipt in range(len(arr_pt_cuts)):
#             NGamma.draw_NGamma_isospin(file, ipt, date)
#     NGamma.draw_NGamma_isospin_combined(file, arr_pt_cuts, date)