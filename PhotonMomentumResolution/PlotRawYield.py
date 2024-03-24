# This code is based on the code provided by Daiki Sekihata in his repository
# https://github.com/dsekihat/photon_sw_run3 for photon analysis in ALICE Run 3


import numpy as np
import datetime
import ROOT
from ROOT import TFile, TDirectory, THashList, TH1F, TH2F, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TGaxis, TPad, TLatex, TBox, TString, TLine
from ROOT import gStyle, gROOT, gSystem
from ROOT import TFile, TDirectory, THashList, TF1, TH1F, TCanvas, TH1
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kFullCircle
import re
import numpy as np
import datetime
import math
import ctypes
gStyle.SetOptStat(0);
gStyle.SetOptTitle(0);
gStyle.SetPadTickX(1);
gStyle.SetPadTickY(1);
gStyle.SetErrorX(0)
gStyle.SetEndErrorSize(5)


from HistoFormatting import FrameSettings, CanvasSettings, PadSettings, DrawHisto, SetTitle, SetStyleTLatex

class PlotRawYieldInvMass:
    def __init__(self):
        print("default constructor is called");
    def __init__(self, meson, filename, dirname):
        print("target meson = {0} , filename = {1} , dirname = {2}".format(meson, filename, dirname));
        self.meson = meson;
        self.rootfile = TFile.Open(filename, "READ");
        self.arr_pt = np.array([0,1,2,3,4,5], dtype=float);
        self.f1total = TF1("GaussExpLinear", 
               "(x<[1]) * ([0]*(TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))))+[4]+[5]*x)+\
               (x>=[1])*([0]*TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))+[4]+[5]*x)", 0,1);
        self.f1total.SetNpx(1000);
        self.fit_min = 0.04;
        self.fit_max = 0.24;
        self.integral_min = 0.18;
        self.integral_max = 0.25;
        self.xtitle = "#it{m}_{#gamma#gamma} (GeV/#it{c}^{2})";
        self.ytitle = "#it{p}_{T,#gamma#gamma} (GeV/#it{c}^{2})";

    def __del__(self):
        if self.rootfile.IsOpen():
            print("close input root file.");
            self.rootfile.Close();

    def set_arr_pt(self, arr_pt):
        print("pT array = ", arr_pt);
        self.arr_pt = arr_pt;

    def set_subsystem(self, ssname):
        self.ssname = ssname;
        self.list_ss   = self.rootfile.Get(ssname);

    def set_cutname(self, cutname):
        self.cutname = cutname;
        self.list_ss_cut = self.list_ss.FindObject(cutname);
    
    def set_fitname(self, fitname):
        self.fitname = fitname;
        self.list_fitname  = self.list_ss_cut.FindObject(fitname);

    def set_fit_range(self, fit_min, fit_max):
        self.fit_min = fit_min;
        self.fit_max = fit_max;

    def set_fit_list(self):
        self.list_fitrange  = self.list_fitname.FindObject("fit_{0:3.2f}_{1:3.2f}_GeVc2".format(self.fit_min, self.fit_max));
        return self.list_fitrange
    
    def make_common_style(self, g1,marker,size,color,width=3,fill=0):
        g1.SetMarkerStyle(marker);
        g1.SetMarkerColor(color);
        g1.SetMarkerSize(size);
        g1.SetLineColor(color);
        g1.SetLineWidth(width);
        g1.SetFillColor(color);
        g1.SetFillStyle(fill);

    def CanvasSettings(self, c1, leftMargin, rightMargin, topMargin, bottomMargin):
        c1.SetTickx();
        c1.SetTicky();
        c1.SetLogy(0);
        c1.SetMargin(leftMargin, rightMargin, topMargin, bottomMargin);
        c1.SetFillColor(0);

    def PadSettings(self, pad1, leftMargin, rightMargin, topMargin, bottomMargin):
        pad1.SetFillColor(0);
        pad1.GetFrame().SetFillColor(0);
        pad1.SetMargin(leftMargin, rightMargin, topMargin, bottomMargin);
        pad1.SetTickx();
        pad1.SetTicky();

    def DrawHisto(self, pad1, histo1, Title, XTitle, YTitle,  markerColor, drawsettings, markerSize =1.):

        histo1.GetYaxis().SetLabelSize(0.02);
        histo1.GetYaxis().SetTitleSize(0.025);
        histo1.GetYaxis().SetDecimals();
        histo1.GetXaxis().SetTitleSize(0.025);
        histo1.GetXaxis().SetLabelSize(0.02);
        histo1.SetMarkerStyle(kFullCircle);
        histo1.GetYaxis().SetLabelSize(0.05);
        histo1.GetXaxis().SetLabelSize(0.05);
        histo1.GetXaxis().SetNdivisions(507, True);

        histo1.SetMarkerStyle(kFullCircle);
        histo1.SetMarkerColor(markerColor);
        histo1.SetMarkerSize(markerSize);
        histo1.SetLineColor(markerColor);
        histo1.SetLineWidth(1);
        histo1.SetFillColor(markerColor);
        histo1.SetFillStyle(0);
        histo1.DrawCopy(drawsettings);

    def SetStyleTLatex( self, text, textSize, lineWidth, textColor = 1, textFont = 42, kNDC = True, align = 11):
        # if kNDC:
        #      text.SetNDC();
        text.SetTextFont(textFont);
        text.SetTextColor(textColor);
        text.SetTextSize(textSize);
        text.SetLineWidth(lineWidth);
        text.SetTextAlign(align);   
        text.SetFillStyle(0)
        text.SetBorderSize(0)


    def PlotLabelsInvMassInPtPlots(self, startTextX, startTextY, textHeight, differenceText, 
                  textAlice, dateDummy, fEnergy, fDecayChannel, fDetectionChannel, 
                  textEvents = TString(""), fNEvents = 0, textAlign = 11):

        alice           = TPaveText(startTextX, startTextY, startTextX + textHeight, startTextY+0.2, "NDC")
        alice.AddText(textAlice)
        self.SetStyleTLatex( alice, textHeight*1.3, 1, 1, 42, True, textAlign);
        alice.Draw();
        latexDate    = TPaveText(startTextX, startTextY, startTextX + textHeight, startTextY+0.2, "NDC")
        latexDate.AddText(dateDummy)
        self.SetStyleTLatex( latexDate, textHeight, 1, 1, 42, True, textAlign);
        latexDate.Draw();
        energy       = TPaveText(startTextX, startTextY, startTextX + textHeight, startTextY+0.2, "NDC")
        energy.AddText(fEnergy)
        self.SetStyleTLatex( energy, textHeight*1, 1, 1, 42, True, textAlign);
        energy.Draw();
        process     = TPaveText(startTextX, startTextY, startTextX + textHeight, startTextY+0.2, "NDC")
        process.AddText(fEnergy)
        self.SetStyleTLatex( process, textHeight*1, 1, 1, 42, True, textAlign);
        process.Draw();
        detprocess  = TPaveText(startTextX, startTextY, startTextX + textHeight, startTextY+0.2, "NDC")
        detprocess.AddText(fDetectionChannel)
        self.SetStyleTLatex( detprocess, textHeight*1, 1, 1, 42, True, textAlign);
        detprocess.Draw();
    
    def GetAndSetLegend2(self, positionX, positionY, positionXRight, positionYUp, textSize, columns = 1, header= TString(""), textFont = 43, margin = 0):

        legend = TLegend(positionX,positionY,positionXRight,positionYUp);
        legend.SetNColumns(columns);
        legend.SetLineColor(0);
        legend.SetLineWidth(0);
        legend.SetFillColor(0);
        legend.SetFillStyle(0);
        legend.SetLineStyle(0);
        legend.SetBorderSize(0);
        legend.SetTextFont(textFont);
        legend.SetTextSize(textSize);
        if margin != 0:
            legend.SetMargin(margin);
        if header.CompareTo("")!= 0:
             legend.SetHeader(header);
        return legend;

    def SetHistoRange(self, iParam):
        yMin_array = [0, 0.1, 0, 0.1*1e-3, 0.1*1e-3, 0, 0, 0]
        yMax_array = [0, 0.15,0, 40*1e-3, 30*1e-3, 0, 0, 0]
        return yMin_array[iParam], yMax_array[iParam]

    def PlotHistoYield(self, fHistoParameter,
                            namePlot, nameCanvas, namePad, Period, numberRowsPlot,
                            numberColumnsPlot,
                            fDecayChannel, fMonteCarloInfo, cutnames, decayChannel = "#gamma#gamma", fDetectionChannel = "#gamma#gamma",
                            fEnergy = "pp at #sqrt{#it{s}} = 13.6 TeV", isVsPtConv = False, BckNmb = 0, 
                            fPlottingType = TString("thesis")):   

        TGaxis.SetMaxDigits(3);
        npt = len(self.arr_pt);

        titlePt = fHistoParameter[0].GetTitle(); # ACCESS TITLE

        yMin_, yMax_ = 1e-8, 1.4*1e-3
        color = [kRed+1, kBlue+1, kGreen+2, kMagenta+2, kCyan+1]

        c1 = TCanvas("c0","c0",0,0,800,800);
        c1.Divide(1,2,1e-3,1e-3);
        p1 = c1.cd(1);
        p1.SetPad(0,0.3,1,1);
        p1.SetMargin(0.15,0.02,0.,0.15);
        p1.SetTicks(1,1);
        p1.SetLogy();


        frame1 = p1.DrawFrame(0., yMin_, 12., yMax_);
        FrameSettings(frame1, "", "#frac{1}{#it{N}_{ev}}#frac{d#it{N}^{#pi^{0}}}{d#it{p}_{T}} counts #upoint (GeV/#it{c})^{-1}",0., 12.)
        frame1.GetXaxis().SetTitleSize(0.045);
        frame1.GetYaxis().SetTitleSize(0.045);
        frame1.GetXaxis().SetTitleOffset(1.0);
        frame1.GetYaxis().SetTitleOffset(1.52);
        frame1.GetXaxis().SetLabelSize(0.045);
        frame1.GetYaxis().SetLabelSize(0.045);
        frame1.GetYaxis().SetMaxDigits(3);
        frame1.GetXaxis().SetLabelOffset(0.01);
        frame1.GetYaxis().SetLabelOffset(0.01);

        for icut in range(len(fHistoParameter)):
            DrawHisto(fHistoParameter[icut], color[icut], "E1,same");

        p2 = c1.cd(2);
        p2.SetPad(0,0,1,0.3);
        p2.SetMargin(0.15,0.02,0.22,0.0);
        p2.SetTicks(1,1);
        p2.SetLogy();
        frame2 = p2.DrawFrame(0., 5*1e-3, 12., 1e1);
        FrameSettings(frame2, "p_{T} (GeV/c)", "#frac{cut}{analysis}",0., 12.)
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

        line1 = TLine(0,1,12.,1);
        line1.SetLineColor(kBlack);
        line1.SetLineStyle(1);
        line1.SetLineWidth(1);
        line1.Draw("");
        ROOT.SetOwnership(line1,False);

        for iratio in range(len(fHistoParameter)-1):
            h1ratio1 = fHistoParameter[iratio+1].Clone("h1ratio");
            h1ratio1.Reset();
            h1ratio1.Sumw2();
            h1ratio1.Divide(fHistoParameter[iratio+1], fHistoParameter[0], 1., 1., "G");
            h1ratio1.DrawCopy("E1,same");
            del h1ratio1

        c1.cd()

        txt = TPaveText(0.,0.9,1.,0.95,"NDC");
        txt.SetFillColor(kWhite);
        txt.SetFillStyle(0);
        txt.SetBorderSize(0);
        txt.SetTextAlign(22);#middle,left
        txt.SetTextFont(42);#helvetica
        txt.SetTextSize(0.045);
        txt.AddText("Raw yield as a function of #it{p}_{T}");
        txt.Draw();
        ROOT.SetOwnership(txt,False);       

        txt = TPaveText(0.65,0.77,0.95,0.87,"NDC");
        txt.SetFillColor(kWhite);
        txt.SetFillStyle(0);
        txt.SetBorderSize(0);
        txt.SetTextAlign(12);#middle,left
        txt.SetTextFont(42);#helvetica
        txt.SetTextSize(0.02);
        txt.AddText("this thesis");
        txt.AddText(Period)
        txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
        txt.Draw();
        ROOT.SetOwnership(txt,False);

        leg = TLegend(0.65,0.62,0.95,0.77);
        leg.SetBorderSize(0);
        leg.SetFillColor(kWhite);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.02);
        leg.SetTextFont(42);#helvetica
        markersize       = fHistoParameter[0].GetMarkerSize();
        for icut in range(len(fHistoParameter)):
            fHistoParameter[icut].SetMarkerSize(markersize);
            leg.AddEntry(fHistoParameter[icut],"{}".format(cutnames[icut]),"ep");
        leg.Draw("");
        ROOT.SetOwnership(leg,False);

        c1.Modified();
        c1.Update();
        ROOT.SetOwnership(c1,False);
        c1.SaveAs(namePlot.Data());
        c1.Close();

        del c1;