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
ROOT.TGaxis.SetMaxDigits(2)

from HistoFormatting import FrameSettings, CanvasSettings, PadSettings, DrawHisto, SetTitle, SetStyleTLatex

class PlotHistoInvMass:
    def __init__(self):
        print("default constructor is called");
    def __init__(self, meson, filename, dirname):
        print("target meson = {0} , filename = {1} , dirname = {2}".format(meson, filename, dirname));
        self.meson = meson;
        self.rootfile = TFile.Open(filename, "READ");
        self.arr_pt = np.array([0,1,2,3,4,5], dtype=float);
        self.fit_min = 0.04;
        self.fit_max = 0.24;
        self.xtitle = "#it{m}_{#gamma#gamma} (GeV/#it{c}^{2})";
        self.ytitle = "#it{p}_{T,#gamma#gamma} (GeV/#it{c}^{2})";

    def __del__(self):
        if self.rootfile.IsOpen():
            print("close input root file.");
            self.rootfile.Close();
 
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
        yMin_array = [1e-8, 120, 0, 0.1*1e-3,0]
        yMax_array = [1e-3, 150, 0, 30*1e-3, 0]
        return yMin_array[iParam], yMax_array[iParam]

    def PlotHistoParameters(self, fHistoParameter,
                            namePlot,  Period, numberRowsPlot,
                            numberColumnsPlot,
                            fDecayChannel, fMonteCarloInfo, decayChannel = "#gamma#gamma", fDetectionChannel = "#gamma#gamma",
                            fEnergy = "pp at #sqrt{#it{s}} = 13.6 TeV", isVsPtConv = False, BckNmb = 0, 
                            fPlottingType = TString("thesis")):   

        if fHistoParameter[len(fHistoParameter)-1].GetEntries() == 0:
            length = len(fHistoParameter)-2
        else:
            length = len(fHistoParameter)
        length = len(fHistoParameter)
        TGaxis.SetMaxDigits(3);
        npt = len(self.arr_pt);
        canvas = TCanvas("c1", "", 3000, 2000)
        canvas.SetTicks(1,1);
        CanvasSettings(canvas, 0, 0, 0.15, 0);
        canvas.cd()
    
        pad = TPad("p1", "", -0.0, 0.0, 1.0, 0.85, 0)
        PadSettings(pad, 0, 0, 0.1, 0);
        pad.Divide(numberColumnsPlot, numberRowsPlot, 0.0, 0.0)
        pad.Draw()

        place        = 0;
        legendPlace  = [numberColumnsPlot, numberRowsPlot]; # right, top
        
        for iParam in range(5):
            place += 1;

            if place == numberColumnsPlot:
                place +=1;
            
            if fHistoParameter[iParam].GetEntries() == 0.0:
                break
            else:
                pass

            pad.cd(place);
            pad.cd(place).SetTopMargin(0.15);
            pad.cd(place).SetBottomMargin(0.15);
            pad.cd(place).SetRightMargin(0.15);
            pad.cd(place).SetLeftMargin(0.15); 

            titlePt = fHistoParameter[iParam].GetTitle();
            yMin, yMax = self.SetHistoRange(iParam);

            if yMax == 0:
                yMax = fHistoParameter[iParam].GetMaximum();
                if yMax >= 0:
                    yMax_ = 1.4*(yMax)
                else:
                    yMax_ = 0.5*(yMax)
            else:
                yMax_ = yMax
            if yMin ==0:
                yMin = fHistoParameter[iParam].GetMinimum();
                if yMin <= 0:
                    yMin_ = 1.4*(yMin)
                else:
                    yMin_ = 0.5*(yMin)
                
                if yMin_ < -1000:
                    yMax_ = +1000.
                if yMax_ > +1000:
                    yMin_ = -500
            else: yMin_ = yMin

            frame = pad.cd(place).DrawFrame(0., yMin_, 12., yMax_);
            FrameSettings(frame, "p_{T} (GeV/c)", "dN{}/dM{}".format(decayChannel, decayChannel),
                        0., 12.)
            if iParam ==0:
                pad.cd(place).SetLogy();
            if iParam == 1:
                fHistoParameter[iParam].Scale(1000.)

            yTitle_list = ["#frac{1}{#it{N}_{ev}}#frac{d#it{N}^{#pi^{0}}}{d#it{p}_{T}} counts #upoint (GeV/#it{c})^{-1}",
                           "MeV/#it{c}^{2}","#lambda (GeV/#it{c}^{2})", "#sigma (GeV/#it{c}^{2})",
                           "fwhm/2.36(GeV/#it{c}^{2})" ]

            frame.GetYaxis().SetTitle(yTitle_list[iParam]);
            DrawHisto(fHistoParameter[iParam], kBlue, "E1 ,same");

            if iParam == 1:
                line2 = TLine(0,134.97,12,134.97);
                line2.SetLineColor(kRed);
                line2.SetLineStyle(2);
                line2.SetLineWidth(1);
                line2.Draw("");
                ROOT.SetOwnership(line2,False); 

            SetTitle(titlePt)

        canvas.cd();

        txt = TPaveText(0.,0.85,1.,0.95,"NDC");
        txt.SetFillColor(kWhite);
        txt.SetFillStyle(0);
        txt.SetBorderSize(0);
        txt.SetTextAlign(22);#middle,left
        txt.SetTextFont(42);#helvetica
        txt.SetTextSize(0.045);
        txt.AddText("Parameters for asymmetric Gaussian fit on invariant mass of #it{#pi}^{0}");
        txt.Draw();
        ROOT.SetOwnership(txt,False);       

        nPixels        = 13;
        textHeight     = 0.08;
        startTextX     = 0.10;
        columnsLegend     = 1;
        widthLegend    = 1./numberColumnsPlot;
        heightLegend   = 1./numberRowsPlot;
        marginWidthLeg = 0.15;
        exampleBin        = numberColumnsPlot-1;
        if numberColumnsPlot > 7:
            startTextX          = 0.05;
            nPixels             = 12;
            widthLegend         = 2./numberColumnsPlot;
            marginWidthLeg      = 0.25;

        # plotting Legend
        padLegend                = TPad("dummyPad","",1-widthLegend,1-heightLegend-0.05,1.,0.85,0);   # gives the size of the histo areas
        PadSettings(padLegend, 0, 0, 0, 0);
        padLegend.Draw();
        padLegend.cd();

        textAlice = "";
        if fPlottingType.CompareTo("wip")==0:
            textAlice       = "ALICE work in progress";
        elif fPlottingType.CompareTo("thesis")==0:
            textAlice       = "this thesis"; # only this thesis
        elif fPlottingType.CompareTo("performance")==0:
            textAlice       = "ALICE performance";
        else:
            textAlice       = "ALICE";

        textEvents="";
        if fMonteCarloInfo:
            textEvents          = "MC";
        else:
            textEvents          = "Data";

        if  padLegend.XtoPixel(padLegend.GetX2()) < padLegend.YtoPixel(padLegend.GetY1()):
            textHeight          = nPixels/padLegend.XtoPixel(padLegend.GetX2()) ;
        else:
            textHeight          = nPixels/padLegend.YtoPixel(padLegend.GetY1());
        
        startTextY     = 0.7;
        differenceText = textHeight*1.05;

        # plot labels
        alice           = TPaveText(startTextX, startTextY, startTextX + 2*textHeight, startTextY - 2*differenceText+0.2, "NDC")
        alice.AddText(textAlice)
        SetStyleTLatex( alice, textHeight*3, 1, 1, 42, True);
        alice.Draw();
        latexPeriod    = TPaveText(startTextX, (startTextY-4*differenceText), startTextX + 2*textHeight, (startTextY-3*differenceText)+0.2, "NDC")
        latexPeriod.AddText(Period)
        SetStyleTLatex( latexPeriod, textHeight*3, 1, 1, 42, True);
        latexPeriod.Draw();
        energy       = TPaveText(startTextX, (startTextY-8*differenceText), startTextX + 2*textHeight, (startTextY-4*differenceText)+0.2, "NDC")
        energy.AddText(fEnergy)
        SetStyleTLatex( energy, textHeight*3, 1, 1, 42, True);
        energy.Draw();
        events = TPaveText(startTextX, (startTextY-12*differenceText), startTextX + 2*textHeight, (startTextY-5*differenceText)+0.2, "NDC")
        events.AddText(fMonteCarloInfo)
        SetStyleTLatex( events, textHeight*3, 1, 1, 42, True);
        events.Draw();
        
        #plot legend
        legendData     = self.GetAndSetLegend2(  startTextX, startTextY-8*differenceText, 0.8,  startTextY-20*differenceText-0.02, 3*nPixels, columnsLegend, TString(""), 43, marginWidthLeg);
        markersize       = fHistoParameter[exampleBin].GetMarkerSize();
        fHistoParameter[exampleBin].SetMarkerSize(markersize);
        legendData.AddEntry(fHistoParameter[exampleBin],"fitparameter for p_{T} ranges","ep");
        linesize         = fHistoParameter[exampleBin].GetLineWidth();
        fHistoParameter[exampleBin].SetLineWidth(linesize);
        legendData.AddEntry(line2, "literature value of #it{#pi}^{0} mass")
        legendData.Draw();
        
        canvas.SaveAs(namePlot.Data());
        ROOT.SetOwnership(canvas, False)
        del padLegend;
        del pad;
        del canvas;

if __name__ == "__main__":
    cuts = True
    generated = True
    cutname = "qc"
    period = "LHC23d1k";
    suffix = "AnyTrack";
    file = "/Users/alicamarieenderich/20231015_invariant_mass_plots/pi0_data_ptspectrum_pp_13.6TeV_LHC22qAnyTrack_copied_run.root"
    dirname = ""
    draw = PlotHistoInvMass("pi0", file, dirname)
    plottingRange = [0., 0.3]