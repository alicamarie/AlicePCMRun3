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
from ROOT import TFile, TDirectory, THashList, TH1F, TH2F, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TLine, TPython
from ROOT import gStyle, gROOT, gSystem
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kTRUE
gStyle.SetOptStat(0);
gStyle.SetOptTitle(0);
from FomattingMaterialBudget import FrameSettings, ALICEtext, FrameSettingsRatio, RatioLegendSettings
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
    x_bins = hist.GetNbinsX()
    mean =  hist.GetMean();
    with open(output_file, 'w') as file:
        for x_bin in range(1, x_bins + 1):
            phi = hist.GetXaxis().GetBinCenter(x_bin)
            nev = hist.GetBinContent(x_bin)

            if nev > 1.1*mean:
                output_line = f"phi: {phi}, nev: {nev}\n"
                file.write(output_line)

class draw_1D_sliced_rxy:
    def __init__(self):
        print("default constructor is called");
    def __init__(self,  config, suffix, folder, period_data, period_mc):
        print("period_data = {0} , period_mc = {1} , config = {2}, suffix = {3}".format(period_data,period_mc, config, suffix));
        self.period_data = period_data;
        self.period_mc = period_mc
        self.suffix = suffix;
        self.config = config
        self.folder = folder;
        
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
    def analyze(self, list_ev, list_v0, arr_rxy, arr_eta):
        h1mult = list_ev.FindObject("hMultNTracksPV").Clone("h1mult");
        nev = h1mult.GetEntries();
        nch = h1mult.GetMean();
        hs = list_v0.FindObject("hs_conv_point").Clone("hs"); #Clone is called to avoid acess deleted object.
        outlist = THashList();
        outlist.SetOwner(True);
        outlist.SetName("outlist");
        outlist.Add(h1mult);

        hs_orig = list_v0.FindObject("hs_conv_point").Clone("hs_orig")
        outlist.Add(hs);
        hs_orig.Sumw2();
        outlist.Add(hs_orig)

        hRxy = hs_orig.Projection(1,"")
        hRxy.Scale(1,"width")
        hRxy.Scale(1./nev);
        hRxy.Scale(1./nch);
        hRxy.SetName("hRxy")
        outlist.Add(hRxy)

        hPhi = hs_orig.Projection(2,"")
        hPhi.Scale(1,"width")
        hPhi.Scale(1./nev);
        hPhi.Scale(1./nch);
        hPhi.SetName("hPhi")
        outlist.Add(hPhi)

        hEta = hs_orig.Projection(3,"")
        hEta.Scale(1,"width")
        hEta.Scale(1./nev);
        hEta.Scale(1./nch);
        hEta.SetName("hEta")
        outlist.Add(hEta)      
          
        arr_phi = np.array(self.get_bin_edges(hs.GetAxis(2)), dtype=float);
        hs.Sumw2();
        hs.Scale(1./nev);
        hs.Scale(1./nch);
        
        eta_bin = [[11, 15],[15, 20], [20, 25], [25, 29], [11, 20], [20, 29]]
        eta_values = [-2.,  -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1., -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,  0.,  0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9,  1.0,  1.1,  1.2,  1.3,  1.4,  1.5,  1.6,  1.7,  1.8,  1.9,  2. ]
        r_bins = [0, 14, 30, 42, 58, 69, 90]

        #first, eta loop
        for ieta in range(0, len(arr_eta)):
            eta1, eta2 = arr_eta[ieta];
            deta = eta2 - eta1;
            hs.GetAxis(3).SetRange(eta1, eta2+1); #don't call SetRangeUser
            h2rphi = hs.Projection(1, 2, "");
            h2rphi.SetName("h2rphi_eta{0:d}".format(ieta));
            h2rphi.SetTitle("conversion point #it{{r}}_{{xy}} vs. #it{{#varphi}} in {0:3.2f} < #it{{#eta}} < {1:3.2f}".format(eta_values[eta1], eta_values[eta2]));
            h2rphi.SetXTitle("#varphi (rad.)");
            h2rphi.SetYTitle("r_{xy} (cm)");
            h2rphi.SetZTitle("#frac{1}{<N_{ch}^{PV}>} #frac{1}{N_{ev}} #frac{dN_{#gamma}}{d#eta}");
            h2rphi.Scale(1/deta);
            outlist.Add(h2rphi);

            #second, rxy loop
            for ir in range(0, len(arr_rxy)-1):
                r1 = arr_rxy[ir];
                r2 = arr_rxy[ir+1];
                dr = r2 - r1;
                bin_r1 = h2rphi.GetYaxis().FindBin(r1 + 1e-6);
                bin_r2 = h2rphi.GetYaxis().FindBin(r2 - 1e-6);

                h1phi = h2rphi.ProjectionX("h1phi_eta{0:d}_r{1:d}".format(ieta, ir), bin_r1, bin_r2, "");
                h1phi.SetTitle("conversion point #it{{#varphi}} in {0:3.2f} < #it{{#eta}} < {1:3.2f} , {2:3.2f} < #it{{r}}_{{xy}} < {3:3.2f} cm".format(eta1, eta2, r1, r2));
                h1phi.SetXTitle("#varphi (rad.)");
                h1phi.SetYTitle("#frac{1}{<N_{ch}^{PV}>} #frac{1}{N_{ev}} #frac{d^{3}N_{#gamma}}{dr_{xy} d#eta d#varphi} (cm #upoint rad.)^{-1}");
                h1phi.Scale(dr);
                h1phi.Scale(1, "width");
                outlist.Add(h1phi);
            hs.GetAxis(3).SetRange(0,0);
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
        outname = os.path.join(self.folder, "{0}_material_budget_dR_{1}_{2}_{3}TeV_{4}{5}.root".format(date, type, self.config["common"]["system"], self.config["common"]["energy"], self.config["common"]["period"], self.suffix));
        outfile = TFile(outname, "RECREATE");

        arr_rxy = self.config["common"]["rxy_bin"];
        arr_eta = self.config["common"]["eta_bin"];

        cutnames = self.config[type]["subsystems"][0]['cutnames']
        nc = len(cutnames);
        for ic in range(0,nc):
            cutname = cutnames[ic];
            list_v0_cut = list_v0.FindObject(cutname);
            outlist = self.analyze(list_ev, list_v0_cut, arr_rxy, arr_eta);
            outlist.SetName(cutname);
            outlist.SetOwner(True);
            outfile.WriteTObject(outlist);
            outlist.Clear();
        
        outfile.Close();
        rootfile.Close();

    def draw_material_phi(self, filename_data, filename_mc, cutname, etaid, rid, date):
        rootfile_data = TFile.Open(filename_data, "READ");
        rootfile_mc   = TFile.Open(filename_mc  , "READ");

        list_data = rootfile_data.Get(cutname);
        list_mc = rootfile_mc.Get(cutname);

        eta_bin = [[11, 15],[15, 20], [20, 25], [25, 29], [11, 20], [20, 29]]
        eta_values = [-2.,  -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1., -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,  0.,  0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9,  1.0,  1.1,  1.2,  1.3,  1.4,  1.5,  1.6,  1.7,  1.8,  1.9,  2. ]
        r_bins = [0, 14, 30, 42, 58, 69, 90, 180]
        eta_min, eta_max = eta_bin[etaid];
            
        h1data_complete = list_data.FindObject("h1phi_eta{0:d}_r{1:d}".format(etaid, rid));
        h1mc_complete = list_mc.FindObject("h1phi_eta{0:d}_r{1:d}".format(etaid, rid));
        h1data_complete.SetDirectory(0);
        h1mc_complete.SetDirectory(0);
    
        #normalization
        h1data_complete.Sumw2()
        h1mc_complete.Sumw2()    
        h1data_complete.Scale(1./(eta_max-eta_min))
        h1mc_complete.Scale(1./(eta_max-eta_min))

        if rid == 4:
            x_bins = h1data_complete.GetNbinsX()
            mean = [15*1e-6, 10*1e-6, 11*1e-6, 18*1e-6, 5*1e-6, 5.5*1e-6]
            output_file = os.path.join(self.folder,"{0}_1D_sliced_in_R_phi_values_eta{1}_r{2}.txt".format(date, etaid, rid));
            with open(output_file, 'w') as file:
                for x_bin in range(1, x_bins + 1):
                    phi = h1data_complete.GetXaxis().GetBinCenter(x_bin)
                    nev_ = h1data_complete.GetBinContent(x_bin)
                    if nev_ > mean[etaid]:
                        output_line = f"phi: {phi}, nev: {nev_}\n"
                        file.write(output_line)

        if rid == 3:
            phi_array = []
            x_bins = h1data_complete.GetNbinsX()
            mean = [0.205*1e-3, 0.159*1e-3, 0.16*1e-3, 0.23*1e-3, 0.07*1e-3, 0.07*1e-3]
            output_file = os.path.join(self.folder,"{0}_1D_sliced_in_R_phi_values_eta{1}_r{2}.txt".format(date, etaid, rid));
            with open(output_file, 'w') as file:
                file.write("bin, phi, nev \n") 
                for x_bin in list(range(12, 26)) + list(range(48, 61)):
                    phi = h1data_complete.GetXaxis().GetBinCenter(x_bin)
                    nev_ = h1data_complete.GetBinContent(x_bin)
                    if nev_ > 1.1*mean[etaid]:
                        output_line = f"{x_bin}, {phi}, {nev_}\n"
                        file.write(output_line)
                        phi_array.append(phi)

        make_common_style(h1data_complete, 20, 1.0, kBlue+1, 1, 0);
        make_common_style(h1mc_complete  , 20, 1.0, kRed+1, 1, 0);
        ROOT.SetOwnership(h1data_complete, False);
        ROOT.SetOwnership(h1mc_complete, False);

        ymax = max(h1data_complete.GetMaximum() , h1mc_complete.GetMaximum()) * 1.7;
        ymin = max(h1data_complete.GetMinimum() , h1mc_complete.GetMinimum()) * -0.2;

        c1 = TCanvas("c0","c0",0,0,800,800);
        c1.Divide(1,2,1e-3,1e-3);
        p1 = c1.cd(1);
        p1.SetPad(0,0.3,1,1);
        p1.SetMargin(0.15,0.02,0.,0.15);
        p1.SetTicks(1,1);

        frame1 = p1.DrawFrame(0, ymin, TMath.TwoPi(), ymax);
        frame1.GetXaxis().SetTitle("conversion point #it{#varphi} (rad.)");
        frame1.GetYaxis().SetTitle("#frac{1}{<#it{N}_{ch}^{PV}>} #frac{1}{#it{N}_{ev}} #frac{d^{3}#it{N}_{#gamma}}{d#it{r}_{xy} d#it{#eta} d#it{#varphi}} (cm #upoint rad.)^{#minus1}");
        FrameSettings(frame1)
        ROOT.SetOwnership(frame1,False);

        h1mc_complete.Draw("E0Hsame");
        h1data_complete.Draw("E0Hsame");

        txt = TPaveText(0.,0.85,1.,0.95,"NDC");
        txt.SetFillColor(kWhite);
        txt.SetFillStyle(0);
        txt.SetBorderSize(0);
        txt.SetTextAlign(22);#middle,left
        txt.SetTextFont(42);#helvetica
        txt.SetTextSize(0.045);
        txt.AddText("Reconstructed photons in {0:2.1f} < #it{{#eta}} < {1:2.1f}".format(eta_values[eta_min], eta_values[eta_max]));
        txt.AddText("and {0:2.1f} < #it{{r}}_{{xy}} < {1:2.1f} cm".format(r_bins[rid], r_bins[rid+1]));
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

        frame2 = p2.DrawFrame(0,cut_in_ratio,TMath.TwoPi(),2.);
        frame2.GetXaxis().SetTitle("#it{#varphi} (rad.)");
        frame2.GetYaxis().SetTitle("#frac{Data}{M.C.}");
        FrameSettingsRatio(frame2)

        h1ratio = h1data_complete.Clone("h1ratio");
        h1ratio.Reset();
        h1ratio.Sumw2();
        h1ratio.Divide(h1data_complete, h1mc_complete, 1., 1., "G");
        h1ratio.Draw("E0same");
        h1ratio.SetDirectory(0);
        ROOT.SetOwnership(h1ratio,False);

        line1 = TLine(0,1,TMath.TwoPi(),1);
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

        leg = RatioLegendSettings()
        leg.AddEntry(h1ratio   ,"Data / M.C. rec.","LP");
        leg.AddEntry(line2  , "ratio \pm 5%", "LP")
        leg.Draw("");
        ROOT.SetOwnership(leg,False);

        c1.Modified();
        c1.Update();
        ROOT.SetOwnership(c1,False);
        filepath = os.path.join(self.folder, "{0}_material_budget_vs_phi_eta{1}_r{2}_{3}.pdf".format(date, etaid, rid, self.suffix));    
        c1.SaveAs(filepath);