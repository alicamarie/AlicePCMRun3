# This code is based on the code provided by Daiki Sekihata in his repository
# https://github.com/dsekihat/photon_sw_run3 for photon analysis in ALICE Run 3


import numpy as np
import datetime
import math
import ROOT
import os
from ROOT import TFile, TDirectory, THashList, TH1F, TH2F, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TLine, TPython
from ROOT import gStyle, gROOT, gSystem, gPad
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kTRUE
gStyle.SetOptStat(0);
gStyle.SetOptTitle(0);
from FomattingMaterialBudget import FrameSettings, ALICEtext, RatioLegendSettings, FrameSettingsRatio
#_____________________________________________________________________
def make_common_style(g1,marker,size,color,width=1,fill=0):
    g1.SetMarkerStyle(marker);
    g1.SetMarkerColor(color);
    g1.SetMarkerSize(size);
    g1.SetLineColor(color);
    g1.SetLineWidth(width);
    g1.SetFillColor(color);
    g1.SetFillStyle(fill);

class draw_1D_integrated:
    def __init__(self):
        print("default constructor is called");
    def __init__(self, filename_data, filename_mc, cutname, folder, period, suffix):
        # self.rootfile_data = TFile.Open(filename_data, "READ");
        # self.rootfile_mc   = TFile.Open(filename_mc  , "READ");

        # self.rootdir_mc_gen  = self.rootfile_mc.Get("pcm-qc-mc")
        # self.list_gen        = self.rootdir_mc_gen.Get("Generated");
        # self.list_ev_mc_gen     = self.rootdir_mc_gen.Get("Event");
        # # self.list_ev_mc_gen  = self.list_ev_gen.FindObject("PCMPCM");
        # self.rootdir_mc_rec  = self.rootfile_mc.Get("material-budget-mc");
        # self.list_v0_mc_rec  = self.rootdir_mc_rec.Get("V0");
        # self.list_ev_rec     = self.rootdir_mc_rec.Get("Event")
        # self.list_ev_mc_rec  = self.list_ev_rec.FindObject("PCMPCM");
        # self.list_cut_mc_rec = self.list_v0_mc_rec.FindObject(cutname);

        # self.h1nch_mc_gen    = self.list_ev_mc_gen.FindObject("hMultNTracksPV").Clone("h1mult");
        # self.nev_gen         = self.h1nch_mc_gen.GetEntries();
        # self.nch_gen         = self.h1nch_mc_gen.GetMean();

        # self.h1nch_mc_rec    = self.list_ev_mc_rec.FindObject("hMultNTracksPV");
        # self.nch_rec         = self.h1nch_mc_rec.GetMean();
        # self.nev_rec         = self.h1nch_mc_rec.GetEntries();
    
        # self.rootdir_data    = self.rootfile_data.Get("material-budget");
        # self.list_v0_data    = self.rootdir_data.Get("V0");
        # self.list_ev_data_1  = self.rootdir_data.Get("Event");
        # self.list_ev_data    = self.list_ev_data_1.FindObject("PCMDalitzEE");
        # self.list_cut_data   = self.list_v0_data.FindObject(cutname);

        # self.h1nch_data      = self.list_ev_data.FindObject("hMultNTracksPV");
        # self.nev_data        = self.h1nch_data.GetEntries();
        # self.nch_data        = self.h1nch_data.GetMean();
    
        self.cutname = cutname;
        self.folder = folder;
        self.period = period;
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
        # print("nev_gen  = {0:e}".format(self.nev_gen));
        # print("nch_gen  = {0:e}".format(self.nch_gen));  

        self.h1nch_mc_rec    = self.list_ev_mc_pcm.FindObject("hMultNTracksPV");
        self.nch_rec         = self.h1nch_mc_rec.GetMean();
        self.nev_rec         = self.h1nch_mc_rec.GetEntries();
        # print("nev_rec  = {0:e}".format(self.nev_rec));
        # print("nch_rec  = {0:e}".format(self.nch_rec));
    
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
        # print("nev_data = {0:e}".format(self.nev_data));
        # print("nch_data = {0:e}".format(self.nch_data)); 
    
        # print("period_data = {0} , period_mc = {1} , config = {2}, suffix = {3}".format(period_data,period_mc, config, suffix));
        # self.period_data = period_data;
        # self.period_mc = period_mc
        self.suffix = suffix;
        # with open(config, "r", encoding="utf-8") as config_yml:
        #     self.config = yaml.safe_load(config_yml)
        self.folder = folder;
        self.cutname = cutname
        # self.config = config
        self.arr_rxy = np.array([0,1,2,3,4,5], dtype=float);

    def __del__(self):
        if self.rootfile_data.IsOpen():
            print("close input data root file.");
            self.rootfile_data.Close();
        if self.rootfile_mc.IsOpen():
            print("close input mc root file.");
            self.rootfile_mc.Close();
    
    #_____________________________________________________________________
    def draw_material_phi(self, cuts, generated, date):

        hs_mc_rec = self.list_cut_mc_rec.FindObject("hs_conv_point").Clone("hs_rec");
        h1_mc_rec = hs_mc_rec.Projection(2 ,"");
        h1_mc_rec.SetDirectory(0);
        ROOT.SetOwnership(h1_mc_rec, False);
        h1_mc_rec.Sumw2();
        #h1_mc_rec.RebinX(2);
        h1_mc_rec.Scale(1,"width");
        #h1_mc_rec.Scale(1/dr);
        h1_mc_rec.Scale(1/self.nev_rec);
        h1_mc_rec.Scale(1/self.nch_rec);#nch
        make_common_style(h1_mc_rec, 20, 1.0, kRed+1, 1, 0);

        h2_mc_gen = self.list_gen_pcmqc.FindObject("hPhotonPhivsRxy");
        h2_mc_gen.Sumw2();
        h2_mc_gen.RebinX(2);
        h2_mc_gen.Scale(1,"width");
        h1_mc_gen = h2_mc_gen.ProjectionX("h1_mc_gen");

        h2_mc_gen.SetDirectory(0);
        h1_mc_gen.SetDirectory(0);
        ROOT.SetOwnership(h2_mc_gen,False);
        ROOT.SetOwnership(h1_mc_gen,False);

        h1_mc_gen.Scale(1/self.nev_gen);
        h1_mc_gen.Scale(1/self.nch_gen);
        make_common_style(h1_mc_gen, 20, 1.0, kGreen+2, 1, 0);

        hs_data = self.list_cut_data.FindObject("hs_conv_point").Clone("hs_data");
        h1_data = hs_data.Projection(2 ,"");
        h1_data.SetDirectory(0);
        ROOT.SetOwnership(h1_data, False);
        h1_data.Sumw2();
        #h1_data.RebinX(2);
        h1_data.Scale(1,"width");
        h1_data.Scale(1/self.nev_data);
        h1_data.Scale(1/self.nch_data);
        make_common_style(h1_data, 20, 1.0, kBlue+1, 1, 0);

    #_____________________________________________________________________
        if generated == True:
            ymax = max(h1_data.GetMaximum() , h1_mc_rec.GetMaximum(), h1_mc_gen.GetMaximum()) * 1.6;
            ymin = max(h1_data.GetMinimum() , h1_mc_rec.GetMinimum(), h1_mc_rec.GetMinimum()) * -0.1;
        else:
            ymax = max(h1_data.GetMaximum() , h1_mc_rec.GetMaximum()) * 1.6;
            ymin = max(h1_data.GetMinimum() , h1_mc_rec.GetMinimum()) * -0.1;
        
        c1 = TCanvas("c0","c0",0,0,800,800);
        c1.Divide(1,2,1e-3,1e-3);
        p1 = c1.cd(1);
        p1.SetPad(0,0.3,1,1);
        p1.SetMargin(0.15,0.02,0.,0.15);
        p1.SetTicks(1,1);
        # p1.SetLogy();
        if generated == True:
            ymin = 1e-3
            ymax = 5*1e-1
        frame1 = p1.DrawFrame(0, ymin, TMath.TwoPi(), ymax);
        frame1.GetXaxis().SetTitle("conversion point #it{#varphi} (rad.)");
        frame1.GetYaxis().SetTitle("#frac{1}{<#it{N}_{ch}^{PV}>} #frac{1}{#it{N}_{ev}} #frac{d#it{N}_{#gamma}}{d#it{#varphi}}(rad.)^{-1}");
        # frame1.GetYaxis().SetTitle("#frac{1}{<#it{N}_{ch}^{PV}>} #frac{1}{#it{N}_{ev}} #frac{d^{2}#it{N}_{#gamma}}{d#it{r}_{xy} d#it{#eta}} (cm)^{#minus1}");
        
        FrameSettings(frame1)
        if generated == True:
            gPad.SetLogy();       
        ROOT.SetOwnership(frame1,False);

        if generated == True:
            h1_mc_gen.Draw("E0h,same");
        h1_mc_rec.Draw("E0h,same");
        h1_data.Draw("E0h,same");
      
        txt = TPaveText(0.,0.9,1.,0.95,"NDC");
        txt.SetFillColor(kWhite);
        txt.SetFillStyle(0);
        txt.SetBorderSize(0);
        txt.SetTextAlign(22);#middle,left
        txt.SetTextFont(42);#helvetica
        txt.SetTextSize(0.045);
        txt.AddText("Reconstructed photons as a function of #it{#varphi}");
        txt.Draw();
        ROOT.SetOwnership(txt,False);       

        ALICEtext("thesis")

        leg = TLegend(0.17,0.72,0.35,0.82);
        leg.SetBorderSize(0);
        leg.SetFillColor(kWhite);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.045);
        leg.AddEntry(h1_data   ,"Data #gamma candidates (LHC22f pass4)","LP");
        if generated == True:
            leg.AddEntry(h1_mc_gen ,"M.C. gen. (LHC23d1k)","LP");
        leg.AddEntry(h1_mc_rec ,"M.C. rec. primary #gamma (LHC23d1k)","LP");
        leg.Draw("");
        ROOT.SetOwnership(leg,False);


        p2 = c1.cd(2);
        p2.SetPad(0,0,1,0.3);
        p2.SetMargin(0.15,0.02,0.22,0.0);
        p2.SetTicks(1,1);
        #p2.SetLogy(1);

        frame2 = p2.DrawFrame(0,0.,TMath.TwoPi(),2.);
        frame2.GetXaxis().SetTitle("#it{#varphi} (rad.)");
        frame2.GetYaxis().SetTitle("#frac{Data}{M.C.}");
        FrameSettingsRatio(frame2)

        line1 = TLine(0,1,TMath.TwoPi(),1);
        line1.SetLineColor(kBlack);
        line1.SetLineStyle(1);
        line1.SetLineWidth(2);
        line1.Draw("");
        ROOT.SetOwnership(line1,False);

            # 5% lines:
        line2 = TLine(0,1.05,TMath.TwoPi(),1.05);
        line2.SetLineColor(kBlack);
        line2.SetLineStyle(2);
        line2.SetLineWidth(2);
        line2.Draw("");
        ROOT.SetOwnership(line2,False); 

        line2 = TLine(0,0.95,TMath.TwoPi(),0.95);
        line2.SetLineColor(kBlack);
        line2.SetLineStyle(2);
        line2.SetLineWidth(2);
        line2.Draw("");
        ROOT.SetOwnership(line2,False); 

        h1ratio = h1_data.Clone("h1ratio");
        #make_common_style(h1ratio, 20, 1.0, kRed+1, 1, 0);
        h1ratio.Reset();
        h1ratio.Sumw2();
        h1ratio.Divide(h1_data, h1_mc_rec, 1., 1., "G");
        h1ratio.Draw("E0,same");
        h1ratio.SetDirectory(0);
        ROOT.SetOwnership(h1ratio,False);

        if generated == True:
            h1ratio2 = h1_data.Clone("h1ratio");
            make_common_style(h1ratio2, 20, 1.0, kGreen+2, 1, 0);
            h1ratio2.Reset();
            h1ratio2.Sumw2();
            h1ratio2.Divide(h1_data, h1_mc_gen, 1., 1., "G");
            h1ratio2.Draw("E0,same");
            h1ratio2.SetDirectory(0);
            ROOT.SetOwnership(h1ratio2,False);

        leg = RatioLegendSettings()
        leg.AddEntry(h1ratio   ,"Data / M.C. rec.","LP");
        #leg.AddEntry(h1ratio1 ,"M.C. gen / M.C. rec.","LP");
        leg.AddEntry(line2  , "ratio \pm 5%", "l")
        # if generated == True:
        #     leg.AddEntry(h1ratio2 ,"Data / M.C. gen.","LP");
        leg.Draw("");
        ROOT.SetOwnership(leg,False);

        # date = datetime.date.today().strftime("%Y%m%d");
        c1.Modified();
        c1.Update();
        ROOT.SetOwnership(c1,False);
        if cuts == True:
            self.suffix += "_with_cuts";
        if generated == True:
            self.suffix += "_with_generated";
        filepath = os.path.join(self.folder, "{0}_pp_13.6TeV_{1}_material_budget_vs_phi_{2}_{3}.pdf".format(date, self.period, self.cutname, self.suffix))
        c1.SaveAs(filepath);

        self.rootfile_data.Close();
        self.rootfile_mc  .Close();
        c1.Close();

        if cuts == True:
            self.suffix = self.suffix.replace("_with_cuts", "");
        if generated == True:
            self.suffix = self.suffix.replace("_with_generated", "");
    #_____________________________________________________________________
    def draw_material_rxy(self, cuts, generated, date):
        hs_mc_rec = self.list_cut_mc_rec.FindObject("hs_conv_point").Clone("hs_rec");
        h1_mc_rec = hs_mc_rec.Projection(1 ,"");
        h1_mc_rec.SetDirectory(0);
        ROOT.SetOwnership(h1_mc_rec, False);
        h1_mc_rec.Sumw2();
        #h1_mc_rec.RebinX(2);
        h1_mc_rec.Scale(1,"width");
        #h1_mc_rec.Scale(1/dr);
        h1_mc_rec.Scale(1/self.nev_rec);
        h1_mc_rec.Scale(1/self.nch_rec);#nch
        h1_mc_rec.GetXaxis().SetRangeUser(0.,90.)
        make_common_style(h1_mc_rec, 20, 1.0, kRed+1, 1, 0)

        h2_mc_gen = self.list_gen_pcmqc.FindObject("hPhotonRZ");
        h2_mc_gen.Sumw2();
        h1_mc_gen = h2_mc_gen.ProjectionY("h1_mc_gen");
        h2_mc_gen.SetDirectory(0);
        h1_mc_gen.SetDirectory(0);
        ROOT.SetOwnership(h2_mc_gen,False);
        ROOT.SetOwnership(h1_mc_gen,False);
        h1_mc_gen.GetXaxis().SetRangeUser(0.,90.)
        h1_mc_gen.RebinX(10);
        h1_mc_gen.Scale(1,"width");
        h1_mc_gen.Scale(1/self.nev_gen);
        h1_mc_gen.Scale(1/self.nch_gen);
        make_common_style(h1_mc_gen, 20, 1.0, kGreen+2, 1, 0);

        hs_data = self.list_cut_data.FindObject("hs_conv_point").Clone("hs_data");
        h1_data = hs_data.Projection(1 ,"");
        h1_data.SetDirectory(0);
        ROOT.SetOwnership(h1_data, False);
        h1_data.Sumw2();
        #h1_data.RebinX(1);
        h1_data.Scale(1,"width");
        h1_data.Scale(1/self.nev_data);
        h1_data.Scale(1/self.nch_data);
        make_common_style(h1_data, 20, 1.0, kBlue+1, 1, 0);    

    #_________________________________________________________________________________________    

        c1 = TCanvas("c0","c0",0,0,800,800);
        c1.Divide(1,2,1e-3,1e-3);
        p1 = c1.cd(1);
        p1.SetPad(0,0.3,1,1);
        p1.SetMargin(0.15,0.02,0.,0.15);
        p1.SetTicks(1,1);
        p1.SetLogy(1);

        frame1 = p1.DrawFrame(0,4e-8,90,1e-1);  #previously 1e-2 as ymax
        frame1.GetXaxis().SetTitle("conversion radius #it{R}_{xy} (cm)");
        frame1.GetYaxis().SetTitle("#frac{1}{<#it{N}_{ch}^{PV}>} #frac{1}{#it{N}_{ev}} #frac{d^{2}#it{N}_{#gamma}}{d#it{r}_{xy} d#it{#eta}} (cm)^{#minus1}");
        FrameSettings(frame1)

        if generated == True:
            h1_mc_gen.Draw("E0h,same");
        h1_mc_rec.Draw("E0h,same");
        h1_data.Draw("E0h,same");
    
        if cuts == True:
            line2 = TLine(42,1e-9,42,1e-1);
            line2.SetLineColor(kMagenta);
            line2.SetLineStyle(2);
            line2.SetLineWidth(2);
            line2.Draw("");
            ROOT.SetOwnership(line2,False);

            line3 = TLine(58,1e-9,58,1e-1);
            line3.SetLineColor(kMagenta);
            line3.SetLineStyle(2);
            line3.SetLineWidth(2);
            line3.Draw("");
            ROOT.SetOwnership(line3,False);

            line7 = TLine(14,1e-9,14,1e-1);
            line7.SetLineColor(kMagenta);
            line7.SetLineStyle(2);
            line7.SetLineWidth(2);
            line7.Draw("");
            ROOT.SetOwnership(line7,False);

            line8 = TLine(30,1e-9,30,1e-1);
            line8.SetLineColor(kMagenta);
            line8.SetLineStyle(2);
            line8.SetLineWidth(2);
            line8.Draw("");
            ROOT.SetOwnership(line8,False);

            line10 = TLine(69,1e-9,69,1e-1);
            line10.SetLineColor(kMagenta);
            line10.SetLineStyle(2);
            line10.SetLineWidth(2);
            line10.Draw("");
            ROOT.SetOwnership(line10,False);

        txt = TPaveText(0.,0.9,1.,0.95,"NDC");
        txt.SetFillColor(kWhite);
        txt.SetFillStyle(0);
        txt.SetBorderSize(0);
        txt.SetTextAlign(22);#middle,left
        txt.SetTextFont(42);#helvetica
        txt.SetTextSize(0.045);
        txt.AddText("Reconstructed photons as a function of #it{R_{xy}}, |#it{#eta}_{#gamma}| < 0.9 ");
        txt.Draw();
        ROOT.SetOwnership(txt,False);

        ALICEtext("thesis")

        leg = TLegend(0.17,0.72,0.35,0.82);
        leg.SetBorderSize(0);
        leg.SetFillColor(kWhite);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.045);
        leg.AddEntry(h1_data   ,"Data #gamma candidates (LHC22f pass4)","LP");
        if generated == True:
            leg.AddEntry(h1_mc_gen ,"M.C. gen. (LHC23d1k)","LP");
        leg.AddEntry(h1_mc_rec ,"M.C. rec. primary #gamma (LHC23d1k)","LP");
        leg.Draw("");
        ROOT.SetOwnership(leg,False);

        p2 = c1.cd(2);
        p2.SetPad(0,0,1,0.3);
        p2.SetMargin(0.15,0.02,0.22,0.0);
        p2.SetTicks(1,1);
        #p2.SetLogy(1);

        frame2 = p2.DrawFrame(0,0.,90,2.);
        frame2.GetXaxis().SetTitle("#it{R}_{xy} (cm)");
        frame2.GetYaxis().SetTitle("#frac{Data}{M.C.}");
        FrameSettingsRatio(frame2)

        line1 = TLine(0,1,90,1);
        line1.SetLineColor(kBlack);
        line1.SetLineStyle(1);
        line1.SetLineWidth(2);
        line1.Draw("");
        ROOT.SetOwnership(line1,False);

    # 5% lines:
        line2 = TLine(0,1.05,90,1.05);
        line2.SetLineColor(kBlack);
        line2.SetLineStyle(2);
        line2.SetLineWidth(2);
        line2.Draw("");
        ROOT.SetOwnership(line2,False); 

        line2 = TLine(0,0.95,90,0.95);
        line2.SetLineColor(kBlack);
        line2.SetLineStyle(2);
        line2.SetLineWidth(2);
        line2.Draw("");
        ROOT.SetOwnership(line2,False);  

        if cuts == True:
            line4 = TLine(42,0,42,4);
            line4.SetLineColor(kMagenta);
            line4.SetLineStyle(2);
            line4.SetLineWidth(2);
            line4.Draw("");
            ROOT.SetOwnership(line4,False);

            line5 = TLine(58,0,58,4);
            line5.SetLineColor(kMagenta);
            line5.SetLineStyle(2);
            line5.SetLineWidth(2);
            line5.Draw("");
            ROOT.SetOwnership(line5,False);

            line6 = TLine(14,0,14,4);
            line6.SetLineColor(kMagenta);
            line6.SetLineStyle(2);
            line6.SetLineWidth(2);
            line6.Draw("");
            ROOT.SetOwnership(line6,False);

            line9 = TLine(30,0,30,4);
            line9.SetLineColor(kMagenta);
            line9.SetLineStyle(2);
            line9.SetLineWidth(2);
            line9.Draw("");
            ROOT.SetOwnership(line9,False);

            line11 = TLine(69,0,69,4);
            line11.SetLineColor(kMagenta);
            line11.SetLineStyle(2);
            line11.SetLineWidth(2);
            line11.Draw("");
            ROOT.SetOwnership(line11,False);

        h1ratio = h1_data.Clone("h1ratio");
        #make_common_style(h1ratio, 20, 1.0, kRed+1, 1, 0);
        h1ratio.Reset();
        h1ratio.Sumw2();
        h1ratio.Divide(h1_data, h1_mc_rec, 1., 1., "G");
        h1ratio.Draw("E0,same");
        h1ratio.SetDirectory(0);
        ROOT.SetOwnership(h1ratio,False);

        h1ratio1 = h1_mc_gen.Clone("h1ratio1");
        make_common_style(h1ratio1, 20, 1.0, kBlack, 1, 0);
        h1ratio1.Reset();
        h1ratio1.Sumw2();
        h1ratio1.Divide(h1_mc_gen, h1_mc_rec, 1., 1., "G");
        #h1ratio1.Draw("E0,same");
        h1ratio1.SetDirectory(0);
        ROOT.SetOwnership(h1ratio1,False);

        if generated == True:
            h1ratio2 = h1_data.Clone("h1ratio");
            make_common_style(h1ratio2, 20, 1.0, kGreen+2, 1, 0);
            h1ratio2.Reset();
            h1ratio2.Sumw2();
            h1ratio2.Divide(h1_data, h1_mc_gen, 1., 1., "G");
            h1ratio2.Draw("E0,same");
            h1ratio2.SetDirectory(0);
            ROOT.SetOwnership(h1ratio2,False);
        
        leg = RatioLegendSettings()
        leg.AddEntry(h1ratio   ,"Data / M.C. rec.","LP");
        #leg.AddEntry(h1ratio1 ,"M.C. gen / M.C. rec.","LP");
        leg.AddEntry(line2  , "ratio \pm 5%", "l")
        # if generated == True:
        #     leg.AddEntry(h1ratio2 ,"Data / M.C. gen.","LP");
        leg.Draw("");
        ROOT.SetOwnership(leg,False);

        # date = datetime.date.today().strftime("%Y%m%d");
        c1.Modified();
        c1.Update();
        ROOT.SetOwnership(c1,False);
        if cuts == True:
            self.suffix += "_with_cuts"
        if generated == True:
            self.suffix += "_with_generated"
        filepath = os.path.join(self.folder, "{0}_pp_13.6TeV_{1}_material_budget_vs_rxy_{2}_{3}.pdf".format(date, self.period, self.cutname, self.suffix))
        c1.SaveAs(filepath);

        self.rootfile_data.Close();
        self.rootfile_mc  .Close();
        c1.Close();

        if cuts == True:
            self.suffix = self.suffix.replace("_with_cuts", "");
        if generated == True:
            self.suffix = self.suffix.replace("_with_generated", "");
    #_____________________________________________________________________
    # add generated when available! right now still missing

    def draw_material_eta(self, cuts, generated, date):

        hs_mc_rec = self.list_cut_mc_rec.FindObject("hs_conv_point").Clone("hs_rec");
        h1_mc_rec = hs_mc_rec.Projection(3 ,"");
        h1_mc_rec.SetDirectory(0);
        ROOT.SetOwnership(h1_mc_rec, False);
        h1_mc_rec.Sumw2();
        #h1_mc_rec.RebinX(2);
        h1_mc_rec.Scale(1,"width");
        #h1_mc_rec.Scale(1/dr);
        h1_mc_rec.Scale(1/self.nev_rec);
        h1_mc_rec.Scale(1/self.nch_rec);#nch
        make_common_style(h1_mc_rec, 20, 1.0, kRed+1, 1, 0)

        hs_data = self.list_cut_data.FindObject("hs_conv_point").Clone("hs_data");
        h1_data = hs_data.Projection(3 ,"");
        h1_data.SetDirectory(0);
        ROOT.SetOwnership(h1_data, False);
        h1_data.Sumw2();
        #h1_data.RebinX(1);
        h1_data.Scale(1,"width");
        h1_data.Scale(1/self.nev_data);
        h1_data.Scale(1/self.nch_data);
        make_common_style(h1_data, 20, 1.0, kBlue+1, 1, 0);   

    #_____________________________________________________________________
        ymax = max(h1_data.GetMaximum() , h1_mc_rec.GetMaximum()) * 1.6;
        ymin = max(h1_data.GetMinimum() , h1_mc_rec.GetMinimum()) * -0.1;

        c1 = TCanvas("c0","c0",0,0,800,800);
        c1.Divide(1,2,1e-3,1e-3);
        p1 = c1.cd(1);
        p1.SetPad(0,0.3,1,1);
        p1.SetMargin(0.15,0.02,0.,0.15);
        p1.SetTicks(1,1);
        #p1.SetLogy(1);

        frame1 = p1.DrawFrame(-1., -0.3*1e-3, 1., ymax);
        frame1.GetXaxis().SetTitle("conversion point #it{#varphi} (rad.)");
        frame1.GetYaxis().SetTitle("#frac{1}{<#it{N}_{ch}^{PV}>} #frac{1}{#it{N}_{ev}} #frac{d#it{N}_{#gamma}}{d#it{#eta}}");
        FrameSettings(frame1)
        #h1_mc_gen.Draw("E0h,same");
        h1_mc_rec.Draw("E0h,same");
        h1_data.Draw("E0h,same");

        if cuts == True:
            rectangle = ROOT.TBox(-1., -0.4*1e-3,-0.9, ymax)  # (x1, y1, x2, y2)
            rectangle.SetFillStyle(3544)  # Set the fill style to hatches
            rectangle.SetLineColor(ROOT.kMagenta)
            rectangle.SetLineWidth(1)
            rectangle.SetFillColor(ROOT.kMagenta) 
            rectangle.Draw()
            ROOT.SetOwnership(rectangle,False);

            rectangle = ROOT.TBox(1., -0.4*1e-3,0.9, ymax)  # (x1, y1, x2, y2)
            rectangle.SetFillStyle(3544)  # Set the fill style to hatches
            rectangle.SetLineColor(ROOT.kMagenta)
            rectangle.SetLineWidth(1)
            rectangle.SetFillColor(ROOT.kMagenta)
            rectangle.Draw()
            ROOT.SetOwnership(rectangle,False);

            line2 = TLine(0.5,-0.4*1e-3,0.5,ymax);
            line2.SetLineColor(kMagenta);
            line2.SetLineStyle(2);
            line2.SetLineWidth(2);
            line2.Draw("");
            ROOT.SetOwnership(line2,False);

            line2 = TLine(-0.5,-0.4*1e-3,-0.5,ymax);
            line2.SetLineColor(kMagenta);
            line2.SetLineStyle(2);
            line2.SetLineWidth(2);
            line2.Draw("");
            ROOT.SetOwnership(line2,False);

            line2 = TLine(0.9,-0.4*1e-3,0.9,ymax);
            line2.SetLineColor(kMagenta);
            line2.SetLineStyle(2);
            line2.SetLineWidth(2);
            line2.Draw("");
            ROOT.SetOwnership(line2,False);

            line2 = TLine(-0.9,-0.4*1e-3,-0.9,ymax);
            line2.SetLineColor(kMagenta);
            line2.SetLineStyle(2);
            line2.SetLineWidth(2);
            line2.Draw("");
            ROOT.SetOwnership(line2,False);

            line2 = TLine(0.0,-0.4*1e-3,0.0,ymax);
            line2.SetLineColor(kMagenta);
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
        txt.AddText("Reconstructed photons as a function of #it{#eta} ");
        txt.Draw();
        ROOT.SetOwnership(txt,False);

        ALICEtext("thesis")

        leg = TLegend(0.17,0.72,0.35,0.82);
        leg.SetBorderSize(0);
        leg.SetFillColor(kWhite);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.045);
        #leg.AddEntry(h1_mc_gen ,"M.C. gen. (LHC23d1k)","LP");
        leg.AddEntry(h1_data   ,"Data #gamma candidates (LHC22f pass4)","LP");
        leg.AddEntry(h1_mc_rec ,"M.C. rec. primary #gamma (LHC23d1k)","LP");

        leg.Draw("");
        ROOT.SetOwnership(leg,False);

        p2 = c1.cd(2);
        p2.SetPad(0,0,1,0.3);
        p2.SetMargin(0.15,0.02,0.22,0.0);
        p2.SetTicks(1,1);
        #p2.SetLogy(1);

        frame2 = p2.DrawFrame(-1.,0.,1.,2.);
        frame2.GetXaxis().SetTitle("#it{#eta}");
        frame2.GetYaxis().SetTitle("#frac{Data}{M.C.}");
        FrameSettingsRatio(frame2)

        line1 = TLine(-1.,1,1.,1);
        line1.SetLineColor(kBlack);
        line1.SetLineStyle(1);
        line1.SetLineWidth(2);
        line1.Draw("");
        ROOT.SetOwnership(line1,False);

        # 5% lines:
        line2 = TLine(-1.,1.05,1.,1.05);
        line2.SetLineColor(kBlack);
        line2.SetLineStyle(2);
        line2.SetLineWidth(2);
        line2.Draw("");
        ROOT.SetOwnership(line2,False); 

        line2 = TLine(-1.,0.95,1.,0.95);
        line2.SetLineColor(kBlack);
        line2.SetLineStyle(2);
        line2.SetLineWidth(2);
        line2.Draw("");
        ROOT.SetOwnership(line2,False); 

        if cuts == True :
            rectangle = ROOT.TBox(-1., 0.,-0.9, 3.2)  # (x1, y1, x2, y2)
            rectangle.SetFillStyle(3544)  # Set the fill style to hatches
            rectangle.SetLineColor(ROOT.kMagenta)
            rectangle.SetLineWidth(1)
            rectangle.SetFillColor(ROOT.kMagenta)  # Set the fill color (white)
            rectangle.Draw()
            ROOT.SetOwnership(rectangle,False);

            rectangle = ROOT.TBox(1., 0. ,0.9, 3.2)  # (x1, y1, x2, y2)
            rectangle.SetFillStyle(3544)  # Set the fill style to hatches
            rectangle.SetLineColor(ROOT.kMagenta)
            rectangle.SetLineWidth(1)
            rectangle.SetFillColor(ROOT.kMagenta)  # Set the fill color (white)
            rectangle.Draw()
            ROOT.SetOwnership(rectangle,False);

            line2 = TLine(0.5,0.,0.5,3.2);
            line2.SetLineColor(kMagenta);
            line2.SetLineStyle(2);
            line2.SetLineWidth(2);
            line2.Draw("");
            ROOT.SetOwnership(line2,False);

            line2 = TLine(-0.5,0.,-0.5,3.2);
            line2.SetLineColor(kMagenta);
            line2.SetLineStyle(2);
            line2.SetLineWidth(2);
            line2.Draw("");
            ROOT.SetOwnership(line2,False);

            line2 = TLine(0.0,0.,0.0,3.2);
            line2.SetLineColor(kMagenta);
            line2.SetLineStyle(2);
            line2.SetLineWidth(2);
            line2.Draw("");
            ROOT.SetOwnership(line2,False);

            line2 = TLine(0.9,0.,0.9,3.2);
            line2.SetLineColor(kMagenta);
            line2.SetLineStyle(2);
            line2.SetLineWidth(2);
            line2.Draw("");
            ROOT.SetOwnership(line2,False);

            line2 = TLine(-0.9,0.,-0.9,3.2);
            line2.SetLineColor(kMagenta);
            line2.SetLineStyle(2);
            line2.SetLineWidth(2);
            line2.Draw("");
            ROOT.SetOwnership(line2,False);
                

        h1ratio = h1_data.Clone("h1ratio");
        #make_common_style(h1ratio, 20, 1.0, kRed+1, 1, 0);
        h1ratio.Reset();
        h1ratio.Sumw2();
        h1ratio.Divide(h1_data, h1_mc_rec, 1., 1., "G");
        h1ratio.Draw("E0,same");
        h1ratio.SetDirectory(0);
        ROOT.SetOwnership(h1ratio,False);

        leg = RatioLegendSettings()
        leg.AddEntry(h1ratio   ,"Data / M.C. rec.","LP");
        #leg.AddEntry(h1ratio1 ,"M.C. gen / M.C. rec.","LP");
        leg.AddEntry(line2  , "ratio \pm 5%", "l")
        #if generated == True:
         #   leg.AddEntry(h1ratio2 ,"Data / M.C. gen.","LP");
        leg.Draw("");
        ROOT.SetOwnership(leg,False);

        # date = datetime.date.today().strftime("%Y%m%d");
        c1.Modified();
        c1.Update();
        ROOT.SetOwnership(c1,False);
        if cuts == True:
            self.suffix += "_with_cuts"
        if generated == True:
            self.suffix += "_with_generated"     
        filepath = os.path.join(self.folder, "{0}_pp_13.6TeV_{1}_material_budget_vs_eta_{2}_{3}.pdf".format(date, self.period, self.cutname, self.suffix))
        c1.SaveAs(filepath)

        self.rootfile_data.Close();
        self.rootfile_mc  .Close();
        c1.Close();
        
        if cuts == True:
            self.suffix = self.suffix.replace("_with_cuts", "");
        if generated == True:
            self.suffix = self.suffix.replace("_with_generated", "");

    def draw_material_nch(self, cuts, date):
        h1nch_mc = self.h1nch_mc_gen #self.list_ev_mc_gen.FindObject("hMultNTracksPV").Clone("h1mult");
        h1nch_mc.SetDirectory(0);
        ROOT.SetOwnership(h1nch_mc, False);
        h1nch_mc.Sumw2();
        #h1_mc_rec.RebinX(2);
        h1nch_mc.Scale(1,"width");
        make_common_style(h1nch_mc, 20, 1.0, kRed+1, 1, 0);

        h1nch_data = self.h1nch_data 
        h1nch_data.SetDirectory(0);
        ROOT.SetOwnership(h1nch_data, False);
        h1nch_data.Sumw2();
        #h1_data.RebinX(2);
        h1nch_data.Scale(1,"width");
        make_common_style(h1nch_data, 20, 1.0, kBlue+1, 1, 0);

    #_____________________________________________________________________
        ymax = max(h1nch_data.GetMaximum() , h1nch_mc.GetMaximum()) * 1.6;
        ymin = max(h1nch_data.GetMinimum() , h1nch_mc.GetMinimum()) * -0.1;

        c1 = TCanvas("c0","c0",0,0,800,800);
        c1.Divide(1,2,1e-3,1e-3);
        p1 = c1.cd(1);
        p1.SetPad(0,0.3,1,1);
        p1.SetMargin(0.15,0.02,0.,0.15);
        p1.SetTicks(1,1);
        p1.SetLogy()

        if ymin <= 0:
            ymin = 1.
        frame1 = p1.DrawFrame(0, ymin, 100, ymax);
        frame1.GetXaxis().SetTitle("#it{N}_{track} to PV");
        frame1.GetYaxis().SetTitle("#it{N}_{ch}");
        FrameSettings(frame1)
        h1nch_mc.Draw("E0h,same");
        h1nch_data.Draw("E0h,same");

        txt = TPaveText(0.,0.9,1.,0.95,"NDC");
        txt.SetFillColor(kWhite);
        txt.SetFillStyle(0);
        txt.SetBorderSize(0);
        txt.SetTextAlign(22);#middle,left
        txt.SetTextFont(42);#helvetica
        txt.SetTextSize(0.045);
        txt.AddText("Number of charged particles");
        txt.Draw();
        ROOT.SetOwnership(txt,False);

        ALICEtext("thesis")
        leg = TLegend(0.17,0.72,0.35,0.82);
        leg.SetBorderSize(0);
        leg.SetFillColor(kWhite);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.045);
        leg.AddEntry(h1nch_data   ,"Data (LHC22f pass4)","LP");
        leg.AddEntry(h1nch_mc ,"M.C. rec. (LHC23d1k)","LP");

        leg.Draw("");
        ROOT.SetOwnership(leg,False);


        p2 = c1.cd(2);
        p2.SetPad(0,0,1,0.3);
        p2.SetMargin(0.15,0.02,0.22,0.0);
        p2.SetTicks(1,1);
        #p2.SetLogy(1);

        frame2 = p2.DrawFrame(0,0.,100,2.);
        frame2.GetXaxis().SetTitle("#it{N}_{track} to PV");
        frame2.GetYaxis().SetTitle("#frac{Data}{M.C.}");
        FrameSettingsRatio(frame2)

        line1 = TLine(0,1,100,1);
        line1.SetLineColor(kBlack);
        line1.SetLineStyle(1);
        line1.SetLineWidth(2);
        line1.Draw("");
        ROOT.SetOwnership(line1,False);

            # 5% lines:
        line2 = TLine(0,1.05,100,1.05);
        line2.SetLineColor(kBlack);
        line2.SetLineStyle(2);
        line2.SetLineWidth(2);
        line2.Draw("");
        ROOT.SetOwnership(line2,False); 

        line2 = TLine(0,0.95,100,0.95);
        line2.SetLineColor(kBlack);
        line2.SetLineStyle(2);
        line2.SetLineWidth(2);
        line2.Draw("");
        ROOT.SetOwnership(line2,False); 

        h1ratio = h1nch_data.Clone("h1ratio");
       #make_common_style(h1ratio, 20, 1.0, kRed+1, 1, 0);
        h1ratio.Reset();
        h1ratio.Sumw2();
        h1ratio.Divide(h1nch_data, h1nch_mc, 1., 1., "G");
        h1ratio.Draw("E0,same");
        h1ratio.SetDirectory(0);
        ROOT.SetOwnership(h1ratio,False);

        leg = RatioLegendSettings()
        leg.AddEntry(h1ratio   ,"Data / M.C. rec.","LP");
        #leg.AddEntry(h1ratio1 ,"M.C. gen / M.C. rec.","LP");
        leg.AddEntry(line2  , "ratio \pm 5%", "l")
        leg.Draw("");
        ROOT.SetOwnership(leg,False);

        # date = datetime.date.today().strftime("%Y%m%d");
        c1.Modified();
        c1.Update();
        ROOT.SetOwnership(c1,False);
        filepath = os.path.join(self.folder, "{0}_pp_13.6TeV_{1}_material_budget_nch_{2}_{3}.pdf".format(date, self.period, self.cutname, self.suffix))
        c1.SaveAs(filepath);

        self.rootfile_data.Close();
        self.rootfile_mc  .Close();
        c1.Close();

    #_____________________________________________________________________

# if __name__ == "__main__":
#     cuts = True
#     generated = True
#     cutname = "qc"
#     period = "LHC23d1k";
#     suffix = "AnyTrack";
#     filename_mc = "/Users/alicamarieenderich/AnalysisResults_LHC23d1k_125889.root"
#     filename_data = "/Users/alicamarieenderich/AnalysisResults_LHC22f4_new_125184.root"
#     date = datetime.date.today().strftime("%Y%m%d");
#     folder = "/Users/alicamarieenderich/{0}_material_budget_plots/".format(date);  
#     os.makedirs(folder, exist_ok=True)
#     draw = draw_1D_integrated(filename_data, filename_mc, cutname, folder, period, suffix)
#     draw.draw_material_rxy(cuts, generated);
#     draw.draw_material_phi(cuts, generated);
#     draw.draw_material_eta(cuts, generated);