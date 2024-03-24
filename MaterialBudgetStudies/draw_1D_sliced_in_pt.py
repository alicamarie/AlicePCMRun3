# This code is based on the code provided by Daiki Sekihata in his repository
# https://github.com/dsekihat/photon_sw_run3 for photon analysis in ALICE Run 3


import os, sys, shutil
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
from ROOT import gStyle, gROOT, gSystem, kTRUE
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange
gStyle.SetOptStat(0);
gStyle.SetOptTitle(0);
from FomattingMaterialBudget import ALICEtext, FrameSettingsRatio, FrameSettings, RatioLegendSettings
#________________________________________________
def make_common_style(g1,marker,size,color,width=1,fill=0):
    g1.SetMarkerStyle(marker);
    g1.SetMarkerColor(color);
    g1.SetMarkerSize(size);
    g1.SetLineColor(color);
    g1.SetLineWidth(width);
    g1.SetFillColorAlpha(color, 0.65);
    g1.SetFillStyle(fill);

class draw_1D_sliced_pt:
    def __init__(self):
        print("default constructor is called");
    def __init__(self,  config, suffix, folder, period_data, period_mc, cutname):
        print("period_data = {0} , period_mc = {1} , config = {2}, suffix = {3}".format(period_data,period_mc, config, suffix));
        self.period_data = period_data;
        self.period_mc = period_mc
        self.suffix = suffix;
        self.folder = folder;
        self.cutname = cutname
        self.config = config

    def __del__(self):
        if self.rootfile.IsOpen():
            print("close input data root file.");
            self.rootfile.Close();
        if self.rootfile.IsOpen():
            print("close input mc root file.");
            self.rootfile.Close();
    
    #________________________________________________
    def get_bin_edges(self, axis):
        list_edges = [];
        for i in range(0, axis.GetNbins()+1):
            list_edges.append(axis.GetBinLowEdge(i+1));
        return list_edges;
    #________________________________________________
    def analyze(self, list_ev, list_v0cut, list_v0, arr_rxy, arr_eta, type):
        h1mult = list_ev.FindObject("hMultNTracksPV").Clone("h1mult");
        nev = h1mult.GetEntries();
        nch = h1mult.GetMean();
        hs = list_v0cut.FindObject("hs_conv_point").Clone("hs"); #Clone is called to avoid acess deleted object.
        outlist = THashList();
        outlist.SetOwner(True);
        outlist.SetName("outlist");
        outlist.Add(h1mult);
        outlist.Add(hs);

        arr_phi = np.array(self.get_bin_edges(hs.GetAxis(2)), dtype=float);
        hs.Sumw2();
        hs.Scale(1./nev);
        hs.Scale(1./nch); 

        list_wire = list_v0.FindObject("wwire_ib");
        hs_wire = list_wire.FindObject("hs_conv_point").Clone("hs_wire")
        hs_wire.Sumw2();
        hs_wire.Scale(1./nev);
        hs_wire.Scale(1./nch); 
        h1wire_pt = hs_wire.Projection(0, "");
        h1wire_pt.SetName("h1wire_pt");
        outlist.Add(h1wire_pt)
    #first, eta loop
        phi1 = 0;
        phi2 = len(arr_phi)-1;
        dphi = phi2 - phi1;
        h2rpt = hs.Projection(1, 0, ""); 
        h2rpt.SetName("h2rpt_phi");
        h2rpt.SetTitle("conversion point #it{{r}}_{{xy}} vs. #it{{p}}_{{T}} in {0:3.2f} < #it{{#varphi}} < {1:3.2f}".format(phi1, phi2));
        h2rpt.SetXTitle("#it{p}_{T} (GeV/c)");
        h2rpt.SetYTitle("r_{xy} (cm)");
        h2rpt.SetZTitle("#frac{1}{<N_{ch}^{PV}>} #frac{1}{N_{ev}} #frac{dN_{#gamma}}{d#it{{p}}_{{T}}}");
        outlist.Add(h2rpt);

        folder = "/Users/alicamarieenderich/this_thesis_material_budget_plots/"
        os.makedirs(folder, exist_ok=True);

        #second, rxy loop
        for ir in range(0, len(arr_rxy)-1):
            r1 = arr_rxy[ir];
            r2 = arr_rxy[ir+1];
            dr = r2 - r1;
            ptmax = h2rpt.GetXaxis().GetNbins()
            pt_bin_cut = h2rpt.GetXaxis().FindBin(5. -1e-6);
            h2rpt.GetXaxis().SetRange(0, pt_bin_cut)
            bin_r1 = h2rpt.GetYaxis().FindBin(r1 + 1e-6);
            bin_r2 = h2rpt.GetYaxis().FindBin(r2 - 1e-6);

            h1pt = h2rpt.ProjectionX("h1pt0_r{0:d}".format(ir), bin_r1, bin_r2, "");
            h1pt.SetTitle("transverse impulse #it{{p}}_{{T}} in {0:3.2f} < #it{{#varphi}} < {1:3.2f} , {2:3.2f} < #it{{r}}_{{xy}} < {3:3.2f} cm".format(phi1, phi2, r1, r2));
            h1pt.SetXTitle("#it{p}_{T}");
            h1pt.SetYTitle("#frac{1}{<N_{ch}^{PV}>} #frac{1}{N_{ev}} #frac{d^{3}N_{#gamma}}{dr_{xy} d#it{p}_{T} d#varphi} (cm #upoint rad.)^{-1}");
            h1pt.Scale(1, "width")
            pt_bin_cut_plus = h2rpt.GetXaxis().FindBin(5. +1e-6);
            h2rpt.GetXaxis().SetRange(0,ptmax)

            h2rpt.GetXaxis().SetRange(pt_bin_cut_plus, ptmax)
            bin_r1 = h2rpt.GetYaxis().FindBin(r1 + 1e-6);
            bin_r2 = h2rpt.GetYaxis().FindBin(r2 - 1e-6);

            h1pt1 = h2rpt.ProjectionX("h1pt1_r{0:d}".format(ir), bin_r1, bin_r2, "");
            h1pt1.SetTitle("transverse impulse #it{{p}}_{{T}} in {0:3.2f} < #it{{#varphi}} < {1:3.2f} , {2:3.2f} < #it{{r}}_{{xy}} < {3:3.2f} cm".format(phi1, phi2, r1, r2));
            h1pt1.SetXTitle("#it{p}_{T}");
            h1pt1.SetYTitle("#frac{1}{<N_{ch}^{PV}>} #frac{1}{N_{ev}} #frac{d^{3}N_{#gamma}}{dr_{xy} d#it{p}_{T} d#varphi} (cm #upoint rad.)^{-1}");
            h1pt1.Scale(1, "width")
            outlist.Add(h1pt);
            outlist.Add(h1pt1)

        hs.GetAxis(2).SetRange(0,0);
        return outlist;

    def run(self, filename, type, date):
        taskname = "material-budget";
        if type == "mc":
            taskname = "material-budget-mc";
        pcmname = "pcm-qc";
        if type == "mc":
            pcmname = "pcm-qc-mc";
        rootfile = TFile.Open(filename, "READ");
        rootdire_pcm = rootfile.Get(pcmname); 
        rootdire = rootfile.Get(taskname); 
        list_v0 = rootdire.Get("V0");
        list_ev     = rootdire_pcm.Get("Event");
        outname = os.path.join(self.folder, "{0}_material_budget_dPt_{1}_{2}_{3}TeV_{4}{5}.root".format(date, type, self.config["common"]["system"], self.config["common"]["energy"], self.config["common"]["period"], self.suffix));
        outfile = TFile(outname, "RECREATE");

        arr_rxy = self.config["common"]["rxy_bin"];
        arr_eta = self.config["common"]["eta_bin"];

        cutnames =self. config[type]["subsystems"][0]['cutnames']
        nc = len(cutnames);
        for ic in range(0,nc):
            cutname = cutnames[ic];
            list_v0_cut = list_v0.FindObject(cutname);
            outlist = self.analyze(list_ev, list_v0_cut, list_v0, arr_rxy, arr_eta, type);
            outlist.SetName(cutname);
            outlist.SetOwner(True);
            outfile.WriteTObject(outlist);
            outlist.Clear();
        outfile.Close();

        rootfile.Close();

    def draw_material_pt_combined(self, filename_data, filename_mc, arr_rxy, cutname, date):
        rootfile_data = TFile.Open(filename_data, "READ");
        rootfile_mc   = TFile.Open(filename_mc  , "READ");

        list_data = rootfile_data.Get(cutname);
        list_mc = rootfile_mc.Get(cutname);

        r_bins = [0, 14, 30, 42, 58, 69, 90]
        
        h1wire_data = list_data.FindObject("h1wire_pt")
        h1wire_mc = list_mc.FindObject("h1wire_pt")
        data_list = []
        mc_list = []
        for rid in range(len(r_bins)-1):
        #pt0
            h1data_0 = list_data.FindObject("h1pt0_r{0:d}".format(rid));
            h1mc_0 = list_mc.FindObject("h1pt0_r{0:d}".format(rid));
            h1data_0.SetDirectory(0);
            h1mc_0.SetDirectory(0);
        
        #normalization
            h1data_0.Sumw2()
            h1mc_0.Sumw2()  

        #style   
            make_common_style(h1data_0, 20, 0.6, kRed+1, 1, 0);
            make_common_style(h1mc_0  , 24, 0.6, kRed+1, 1, 0);
            #h1mc_complete.SetLineStyle(2);
            ROOT.SetOwnership(h1data_0, False);
            ROOT.SetOwnership(h1mc_0, False);

        #pt1
            h1data_1 = list_data.FindObject("h1pt1_r{0:d}".format(rid));
            h1mc_1 = list_mc.FindObject("h1pt1_r{0:d}".format(rid));
            h1data_1.SetDirectory(0);
            h1mc_1.SetDirectory(0);
        
        #normalization
            h1data_1.Sumw2()
            h1mc_1.Sumw2()  

        #style   
            make_common_style(h1data_1, 20, 0.6, kRed+1, 1, 0);
            make_common_style(h1mc_1  , 24, 0.6, kRed+1, 1, 0);
            #h1mc_complete.SetLineStyle(2);
            ROOT.SetOwnership(h1data_1, False);
            ROOT.SetOwnership(h1mc_1, False);
        
            data_list.append([h1data_0, h1data_1])
            mc_list.append([h1mc_0, h1mc_1])

        ymax = max(h1data_0.GetMaximum(), h1mc_0.GetMaximum(), h1data_1.GetMaximum(), h1mc_1.GetMaximum())* 1.2;
        ymin = min(h1data_0.GetMinimum(), h1mc_0.GetMinimum(),h1data_1.GetMinimum(), h1mc_1.GetMinimum())* -0.3;
        print("YMAX, ymax", ymax)

    #canvas plotting
        c1 = TCanvas("c0","c0",0,0,800,800);
        c1.Divide(1,2,1e-3,1e-3);
        p1 = c1.cd(1);
        p1.SetPad(0,0.3,1,1);
        p1.SetMargin(0.15,0.02,0.,0.15);
        p1.SetTicks(1,1);
        p1.SetLogy()
        #p1.SetLogx()

        frame1 = p1.DrawFrame(5*1e-2, 1e-9, 10., 5*1e-2);
        frame1.GetXaxis().SetTitle("transverse impulse #it{p}_{T} (GeV/c)");
        frame1.GetYaxis().SetTitle("#frac{1}{<#it{N}_{ch}^{PV}>} #frac{1}{#it{N}_{ev}} #frac{d#it{N}_{#gamma}}{d#it{p_{T}}} (GeV/c)^{#minus1}");
        FrameSettings(frame1)

        color = [kRed+1, kOrange+1, kYellow+1, kGreen+2, kBlue+1, kCyan+2, kMagenta+2]
        color = [kMagenta+2, kOrange+1, kYellow+1, kRed+1, kBlue+1, kGreen+2, kCyan+1,]
        for i in range(len(data_list) - 1, -1, -1):
            
            for j in range(2):
                make_common_style(data_list[i][j], 20, 0.6, color[i], 1, 0);
                make_common_style(mc_list[i][j]  , 25, 0.6, color[i], 1, 0);
                data_list[i][j].Draw("E0same");
                mc_list[i][j].Draw("E0same");

        make_common_style(h1wire_data, 20, 0.6, kCyan+1, 1, 0)
        make_common_style(h1wire_mc, 25, 0.6, kCyan+1, 1, 0)
        h1wire_data.Draw("E0same")
        h1wire_mc.Draw("E0same")

        txt = TPaveText(0.,0.9,1.,0.95,"NDC");
        txt.SetFillColor(kWhite);
        txt.SetFillStyle(0);
        txt.SetBorderSize(0);
        txt.SetTextAlign(22);#middle,left
        txt.SetTextFont(42);#helvetica
        txt.SetTextSize(0.045);
        txt.AddText("#it{p}_{T} spectra for different radii");
        txt.Draw();
        ROOT.SetOwnership(txt,False);

        txt = TPaveText(0.90,0.72,0.95,0.82,"NDC");
        txt.SetFillColor(kWhite);
        txt.SetFillStyle(0);
        txt.SetBorderSize(0);
        txt.SetTextAlign(32);#middle,left
        txt.SetTextFont(42);#helvetica
        txt.SetTextSize(0.03);
        # txt.AddText("ALICE");
        txt.AddText("this thesis")
        txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
        txt.Draw();
        ROOT.SetOwnership(txt,False);

        leg = TLegend(0.70,0.50,0.95,0.70);
        leg.SetBorderSize(0);
        leg.SetFillColor(kWhite);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.03);
        leg.SetTextAlign(32);
        leg.SetTextFont(42);#helvetica
        for rid in range(len(r_bins)-1):
            leg.AddEntry(data_list[rid][0], "{0} cm < #it{{r}}_{{xy}} < {1} cm".format(r_bins[rid], r_bins[rid+1]),"LP");
        leg.AddEntry(h1wire_data, "wire data")
        leg.AddEntry(h1wire_mc, "wire mc")
        leg.Draw("");
        ROOT.SetOwnership(leg,False);

        leg = TLegend(0.25,0.72,0.40,0.82);
        leg.SetBorderSize(0);
        leg.SetFillColor(kWhite);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.03);
        leg.SetTextAlign(12);
        leg.SetTextFont(42);#helvetica
        leg.AddEntry(data_list[3][0], "Data #gamma candidates (LHC22f)","LP");
        leg.AddEntry(mc_list[3][0] , "M.C. rec. #gamma (LHC23d1k)","LP");
        leg.Draw("");
        ROOT.SetOwnership(leg,False);

        p2 = c1.cd(2);
        p2.SetPad(0,0,1,0.3);
        p2.SetMargin(0.15,0.02,0.22,0.0);
        p2.SetTicks(1,1);

        frame2 = p2.DrawFrame(5*1e-2,0.2,10.,2.4);
        frame2.GetXaxis().SetTitle("#it{p}_{T} (GeV/c)");
        frame2.GetYaxis().SetTitle("#frac{Data}{M.C.}");

        FrameSettingsRatio(frame2)

        color = [kRed+1, kOrange+7, kYellow+1, kGreen+2, kBlue+1,  kMagenta+2, kCyan+1, kCyan+4]
        for i in range(len(data_list)):
            for j in range(2):
                h1ratio = data_list[i][j].Clone("h1ratio");
                h1ratio.Reset();
                h1ratio.Sumw2();
                h1ratio.Divide(data_list[i][j], mc_list[i][j], 1., 1., "G");
                h1ratio.Draw("E0same");
                h1ratio.SetDirectory(0);
                ROOT.SetOwnership(h1ratio,False);
        
        h1ratio_wire = h1wire_data.Clone("h1ratio_wire");
        h1ratio_wire.Reset();
        h1ratio_wire.Sumw2();
        h1ratio_wire.Divide(h1wire_data, h1wire_mc, 1., 1., "G");
        h1ratio_wire.Draw("E0same");
        h1ratio_wire.SetDirectory(0);
        ROOT.SetOwnership(h1ratio_wire,False);

        line1 = TLine(0.,1,10.,1);
        line1.SetLineColor(kBlack);
        line1.SetLineStyle(1);
        line1.SetLineWidth(2);
        line1.Draw("");
        ROOT.SetOwnership(line1,False);

        c1.Update();
        ROOT.SetOwnership(c1,False);
        filepath = os.path.join(self.folder, "{0}_material_budget_vs_Pt_combined_{1}.pdf".format(date, self.suffix));    
        c1.SaveAs(filepath);

    def draw_material_pt(self, filename_data, filename_mc, rid, cutname, date):
        rootfile_data = TFile.Open(filename_data, "READ");
        rootfile_mc   = TFile.Open(filename_mc  , "READ");

        list_data = rootfile_data.Get(cutname);
        list_mc = rootfile_mc.Get(cutname);
        r_bins = [0, 14, 30, 42, 58, 69, 90, 180]
        
        h1data_complete = list_data.FindObject("h1pt0_r{0:d}".format(rid));
        h1mc_complete = list_mc.FindObject("h1pt0_r{0:d}".format(rid));
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

        ymax = max(h1data_complete.GetMaximum() , h1mc_complete.GetMaximum()) * 1.2;
        ymin = max(h1data_complete.GetMinimum() , h1mc_complete.GetMinimum()) * 0.1;

    #canvas plotting
        c1 = TCanvas("c0","c0",0,0,800,800);
        c1.Divide(1,2,1e-3,1e-3);
        p1 = c1.cd(1);
        p1.SetPad(0,0.3,1,1);
        p1.SetMargin(0.15,0.02,0.,0.15);
        p1.SetTicks(1,1);
        # p1.SetLogy();

        frame1 = p1.DrawFrame(0., ymin, 10., ymax); #(0., 1e-20, 10., 1e-1);#
        frame1.GetXaxis().SetTitle("transverse impulse #it{p}_{T} (GeV/c)");
        frame1.GetYaxis().SetTitle("#frac{1}{<#it{N}_{ch}^{PV}>} #frac{1}{#it{N}_{ev}} #frac{d#it{N}_{#gamma}}{d#it{p_{T}}} (GeV/c)^{#minus1}");
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
        txt.AddText("Reconstructed photons in {0:2.1f} < #it{{r}}_{{xy}} < {1:2.1f} cm".format(r_bins[rid], r_bins[rid+1]));
        txt.Draw();
        ROOT.SetOwnership(txt,False);

        ALICEtext("thesis")

        leg = TLegend(0.17,0.72,0.35,0.82);
        leg.SetBorderSize(0);
        leg.SetFillColor(kWhite);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.045);
        leg.AddEntry(h1data_complete, "Data #gamma candidates (LHC22f)","LP");
        leg.AddEntry(h1mc_complete  , "M.C. rec. #gamma (LHC23d1k)","LP");
        leg.Draw("");
        ROOT.SetOwnership(leg,False);

        p2 = c1.cd(2);
        p2.SetPad(0,0,1,0.3);
        p2.SetMargin(0.15,0.02,0.22,0.0);
        p2.SetTicks(1,1);
        cut_in_ratio = 0.5
        if rid ==5:
            cut_in_ratio = 0.0

        frame2 = p2.DrawFrame(0.,cut_in_ratio,10.,2.);
        frame2.GetXaxis().SetTitle("#it{p}_{T} (GeV/c)");
        frame2.GetYaxis().SetTitle("#frac{Data}{M.C.}");

        FrameSettingsRatio(frame2)

        h1ratio = h1data_complete.Clone("h1ratio");
        h1ratio.Reset();
        h1ratio.Sumw2();
        h1ratio.Divide(h1data_complete, h1mc_complete, 1., 1., "G");
        h1ratio.Draw("E0same");
        h1ratio.SetDirectory(0);
        ROOT.SetOwnership(h1ratio,False);

        line1 = TLine(0.,1,10.,1);
        line1.SetLineColor(kBlack);
        line1.SetLineStyle(1);
        line1.SetLineWidth(2);
        line1.Draw("");
        ROOT.SetOwnership(line1,False);

        leg = TLegend(0.17,0.82,0.35,0.97);
        leg.SetBorderSize(0);
        leg.SetFillColor(kWhite);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.05);
        leg.AddEntry(h1ratio   ,"Data / M.C. rec.","LP");
        leg.Draw("");
        ROOT.SetOwnership(leg,False);

        c1.Modified();
        c1.Update();
        ROOT.SetOwnership(c1,False);
        filepath = os.path.join(self.folder, "{0}_material_budget_vs_Pt_r{1}_{2}.pdf".format(date, rid, self.suffix));    
        c1.SaveAs(filepath);