# This code is based on the code provided by Daiki Sekihata in his repository
# https://github.com/dsekihat/photon_sw_run3 for photon analysis in ALICE Run 3


import os, sys, shutil
import math
import argparse
import numpy as np
import ctypes
import datetime
import yaml
import ROOT
ROOT.gROOT.SetBatch(True);
from ROOT import TFile, THashList, TF1
from ROOT import TFile, TDirectory, THashList, TH1F, TH2F, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TLine, TPython
from ROOT import gStyle, gROOT, gSystem, gPad
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan
gStyle.SetOptStat(0);
gStyle.SetOptTitle(0);

#________________________________________________
def make_common_style(g1,marker,size,color,width=1,fill=0):
    g1.SetMarkerStyle(marker);
    g1.SetMarkerColor(color);
    g1.SetMarkerSize(size);
    g1.SetLineColor(color);
    g1.SetLineWidth(width);
    g1.SetFillColor(color);
    g1.SetFillStyle(fill);

class make_2D_Eta_vs_Phi:
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

    def analyze(self, list_ev, list_v0, arr_rxy):
        h1mult = list_ev.FindObject("hMultNTracksPV").Clone("h1mult");
        nev = h1mult.GetEntries();
        nch = h1mult.GetMean();
        hs = list_v0.FindObject("hs_conv_point").Clone("hs"); #Clone is called to avoid acess deleted object.
        outlist = THashList();
        outlist.SetOwner(True);
        outlist.SetName("outlist");
        outlist.Add(h1mult);
        outlist.Add(hs);
        arr_phi = np.array(self.get_bin_edges(hs.GetAxis(1)), dtype=float);
        arr_eta = np.array(self.get_bin_edges(hs.GetAxis(2)), dtype=float);
        hs.Sumw2();
        hs.Scale(1./nev);
        hs.Scale(1./nch);

        #first, eta loop
        for ir in range(0, len(arr_rxy)-1):
            r1 = arr_rxy[ir];
            r2 = arr_rxy[ir+1];
            dr = r2 - r1;
            hs.GetAxis(1).SetRange(r1, r2); #don't call SetRangeUser
            h2etaphi = hs.Projection(3, 2, "");
            h2etaphi.SetName("h2etaphi_r{0:d}".format(ir));
            h2etaphi.SetTitle("conversion point #it{{#eta}} vs. #it{{#varphi}} in {0:3.2f} < #it{{r}}_{{xy}} < {1:3.2f}".format(r1, r2));
            h2etaphi.SetXTitle("#varphi (rad.)");
            h2etaphi.SetYTitle("#eta");
            h2etaphi.SetZTitle("#frac{1}{<N_{ch}^{PV}>} #frac{1}{N_{ev}} #frac{dN_{#gamma}}{dr_{xy}}");
            h2etaphi.Scale(1/dr);
            outlist.Add(h2etaphi);

        hs.GetAxis(2).SetRange(0,0);
        return outlist;
    #________________________________________________
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
        outname = os.path.join(self.folder, "{0}_material_budget_eta_vs_phi_{1}_{2}_{3}TeV_{4}{5}.root".format(date, type, self.config["common"]["system"], self.config["common"]["energy"], self.config["common"]["period"], self.suffix));
        outfile = TFile(outname, "RECREATE");

        arr_rxy = self.config["common"]["rxy_bin"];

        cutnames = self.config[type]["subsystems"][0]['cutnames']
        nc = len(cutnames);
        for ic in range(0,nc):
            cutname = cutnames[ic];
            list_v0_cut = list_v0.FindObject(cutname);
            outlist = self.analyze(list_ev, list_v0_cut, arr_rxy);
            outlist.SetName(cutname);
            outlist.SetOwner(True);
            outfile.WriteTObject(outlist);
            outlist.Clear();
        outfile.Close();

        rootfile.Close();

    def find_phi_for_eta(self, hist, eta_value, output_file):
        eta_bin = hist.GetYaxis().FindBin(eta_value)
        x_bins = hist.GetNbinsX()

        with open(output_file, 'w') as file:
            for x_bin in range(1, x_bins + 1):
                phi = hist.GetXaxis().GetBinCenter(x_bin)
                nev = hist.GetBinContent(x_bin, eta_bin)

                if nev > 4e-7:
                    output_line = f"phi: {phi}, nev: {nev}\n"
                    file.write(output_line)

    def draw_material_z_vs_phi(self, filename_data, filename_mc, rid, cutname, log):
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
        #p1.SetPad(0,0.35,1,1);
        p1.SetMargin(0.1,0.13,0.13,0.2);
        p1.SetTicks(1,1);


        frame1 = p1.DrawFrame(0, -2, TMath.TwoPi(), 2);
        frame1.GetXaxis().SetTitle("#it{#varphi} (rad.)");
        frame1.GetYaxis().SetTitle("pseudorapidity #it{#eta}");
        frame1.GetXaxis().SetTitleSize(0.05);
        frame1.GetYaxis().SetTitleSize(0.05);
        frame1.GetXaxis().SetTitleOffset(1.0);
        frame1.GetYaxis().SetTitleOffset(1.0);
        frame1.GetXaxis().SetLabelSize(0.05);
        frame1.GetYaxis().SetLabelSize(0.05);
        frame1.GetYaxis().SetMaxDigits(3);
        frame1.GetXaxis().SetLabelOffset(0.01);
        frame1.GetYaxis().SetLabelOffset(0.01);
        ROOT.SetOwnership(frame1,False);

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
        txt.AddText("#it{{#eta}} vs. #it{{#varphi}} for pp at #sqrt{{#it{{s}}}}  = 13.6 TeV in {0:3.2f} cm < #it{{r}}_{{xy}} < {1:3.2f} cm".format(r_bins[rid], r_bins[rid+1]))
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

        p2 = c1.cd(2);
        p2.SetMargin(0.1,0.13,0.22,0.05);
        p2.SetTicks(1,1);

        frame2 = p2.DrawFrame(0, -2, TMath.TwoPi(), 2);
        frame2.GetXaxis().SetTitle("#it{#varphi} (rad.)");
        frame2.GetYaxis().SetTitle("pseudorapidity #it{#eta}");
        frame2.GetXaxis().SetTitleSize(0.05);
        frame2.GetYaxis().SetTitleSize(0.05);
        frame2.GetXaxis().SetTitleOffset(1.0);
        frame2.GetYaxis().SetTitleOffset(1.0);
        frame2.GetXaxis().SetLabelSize(0.05);
        frame2.GetYaxis().SetLabelSize(0.05);
        frame2.GetYaxis().SetMaxDigits(3);
        frame2.GetXaxis().SetLabelOffset(0.01);
        frame2.GetYaxis().SetLabelOffset(0.01);
        ROOT.SetOwnership(frame2,False);

        if log == "log":
            gPad.SetLogz()
        h2mc_complete.SetMaximum(max_z_value_data)
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

        date = datetime.date.today().strftime("%Y%m%d");
        c1.Modified();
        c1.Update();
        ROOT.SetOwnership(c1,False);
        if log == "log":
            outname = os.path.join(self.folder,"{0}_material_budget_eta_vs_phi_log_r{1}.png".format(date, rid))
            c1.SaveAs(outname);

        elif log == "":
            outname = os.path.join(self.folder, "{0}_material_budget_eta_vs_phi_r{1}.png".format(date, rid))
            c1.SaveAs(outname);

        if rid == 3:
            self.find_phi_for_eta(h2data_complete,0, os.path.join(self.folder, "{0}_2D_Eta_vs_Phi_phi_values_for_r3.txt".format(date)));