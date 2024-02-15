# This code is based on the code provided by Daiki Sekihata in his repository
# https://github.com/dsekihat/photon_sw_run3 for photon analysis in ALICE Run 3


import numpy as np
import datetime
import ROOT
from ROOT import TFile, TDirectory, THashList, TH1F, TH2F, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TGaxis, TPad, TLatex, TBox, TString, TLine
from ROOT import gStyle, gROOT, gSystem
from ROOT import TFile, TDirectory, THashList, TF1, TH1F, TCanvas, TH1
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kFullCircle, kGray
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

class PlotSameMixed:
    def __init__(self):
        print("default constructor is called");
    def __init__(self, meson, filename, dirname):
        print("target meson = {0} , filename = {1} , dirname = {2}".format(meson, filename, dirname));
        self.meson = meson;
        self.rootfile = TFile.Open(filename, "READ");
        self.arr_pt = np.array([0.4, 0.8, 1.2, 1.6, 2.0, 3.0, 4.0, 6.0, 12.0], dtype=float);
        self.f1total = TF1("GaussExpLinear", 
               "(x<[1]) * ([0]*(TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))))+[4]+[5]*x)+\
               (x>=[1])*([0]*TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))+[4]+[5]*x)", 0,1);
        self.f1total.SetNpx(1000);
        self.fit_min = 0.06;
        self.fit_max = 0.20;
        self.integral_min = 0.18;
        self.integral_max = 0.25;
        self.yield_min = 0.035;
        self.yield_max= 0.01;
        self.xtitle = "#it{m}_{#gamma#gamma} (GeV/#it{c}^{2})";
        self.ytitle = "#it{p}_{T,#gamma#gamma} (GeV/#it{c}^{2})";

    def __del__(self):
        if self.rootfile.IsOpen():
            print("close input root file.");
            self.rootfile.Close();

    def set_arr_pt(self, arr_pt):
        print("pT array = ", arr_pt);
        self.arr_pt = arr_pt;

    def set_fit_range(self, fit_min, fit_max):
        self.fit_min = fit_min;
        self.fit_max = fit_max;
    
    def set_integral_range(self, integral_min, integral_max):
        self.integral_min = integral_min;
        self.integral_max = integral_max;
    
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

    def PlotSameMixedInPtBins(self, fHistoSamePtBinPlot, fHistoMixedPtBinPlot, 
                            namePlot, nameCanvas, namePad, fPlottingRangeMeson, Period, numberRowsPlot,
                            numberColumnsPlot, fStartBinPtRange, fNumberPtBins, fRangeBinsPt,
                            fDecayChannel, fMonteCarloInfo, decayChannel = "#gamma#gamma", fDetectionChannel = "#gamma#gamma",
                            fEnergy = "pp at #sqrt{#it{s}} = 13.6 TeV", optionBackground = "Asymmetric Gaussian fit", isVsPtConv = False, BckNmb = 0, 
                            fPlottingType = TString("thesis")):   
                 
        TGaxis.SetMaxDigits(3);
        npt = len(self.arr_pt);
        canvas = TCanvas("c1", "", 3000, 2000)
        canvas.SetTicks(1,1);
        CanvasSettings(canvas, 0, 0, 0.15, 0);
        canvas.cd()
    
        pad = TPad("p1", "", -0.0, 0.0, 1.0,0.85, 0)
        PadSettings(pad, 0, 0, 0.1, 0);
        pad.Divide(numberColumnsPlot, numberRowsPlot, 0.0, 0.0)
        pad.Draw()

        place        = 0;
        legendPlace  = [numberColumnsPlot, numberRowsPlot]; # right, top
            
        for iPt in range(fStartBinPtRange, fRangeBinsPt-1):
            startPt = self.arr_pt[iPt];
            endPt = self.arr_pt[iPt+1];
            place += 1;
            
            fit_min = self.fit_min
            fit_max = self.fit_max
            
            if place == numberColumnsPlot:
                place +=1

            pad.cd(place);
            pad.cd(place).SetTopMargin(0.15);
            pad.cd(place).SetBottomMargin(0.15);
            pad.cd(place).SetRightMargin(0.15);
            pad.cd(place).SetLeftMargin(0.15); # -> change values later to make neater

            titlePt = "{:.2f} GeV/c < #it{{p}}_{{T}} < {:.2f} GeV/c".format(startPt, endPt);
            
            yMin = 0.;
            yMax = 0.;
            if isinstance(fHistoSamePtBinPlot[iPt], TH1):
                for i in range ( fHistoSamePtBinPlot[iPt].GetXaxis().FindBin(fPlottingRangeMeson[0]), fHistoSamePtBinPlot[iPt].GetXaxis().FindBin(fPlottingRangeMeson[1])):
                    if fHistoSamePtBinPlot[iPt].GetBinContent(i) < yMin:
                        yMin = fHistoSamePtBinPlot[iPt].GetBinContent(i);
                    
                    if fHistoSamePtBinPlot[iPt].GetBinContent(i) > yMax:
                        yMax = fHistoSamePtBinPlot[iPt].GetBinContent(i);
            else:
                yMax = fHistoSamePtBinPlot[iPt].GetMaximum(fPlottingRangeMeson[0], fPlottingRangeMeson[1]);
                yMin = fHistoSamePtBinPlot[iPt].GetMinimum(fPlottingRangeMeson[0], fPlottingRangeMeson[1]);
            
            frame = pad.cd(place).DrawFrame(fPlottingRangeMeson[0], yMin, fPlottingRangeMeson[1], 1.6*yMax);

            FrameSettings(frame, "M{} (GeV/#it{{c}}^{{2}})".format(decayChannel), "dN{}/dM{}".format(decayChannel, decayChannel),
                        fPlottingRangeMeson[0], fPlottingRangeMeson[1])

            DrawHisto(fHistoSamePtBinPlot[iPt], kBlue, "Hsame", 0.9, 24);
            DrawHisto(fHistoMixedPtBinPlot[iPt], kRed, "Hsame", 0.9, 24);
        #DrawHisto(list_data0[0], kGreen+2, "Hsame", 0.9, 24); 
            # line2 = TLine(fit_min,yMin,fit_min,1.6*yMax);
            # line2.SetLineColor(kRed);
            # line2.SetLineStyle(2);
            # line2.SetLineWidth(1);
            # line2.Draw("");
            # ROOT.SetOwnership(line2,False); 

            # line2 = TLine(fit_max,yMin,fit_max,1.6*yMax);
            # line2.SetLineColor(kRed);
            # line2.SetLineStyle(2);
            # line2.SetLineWidth(1);
            # line2.Draw("");
            # ROOT.SetOwnership(line2,False); 

            line4 = TLine(self.integral_min,yMin,self.integral_min,1.6*yMax);
            line4.SetLineColor(kGray+2);
            line4.SetLineStyle(1);
            line4.SetLineWidth(1);
            line4.Draw("");
            ROOT.SetOwnership(line4,False); 

            line4 = TLine(self.integral_max,yMin,self.integral_max,1.6*yMax);
            line4.SetLineColor(kGray+2);
            line4.SetLineStyle(1);
            line4.SetLineWidth(1);
            line4.Draw("");
            ROOT.SetOwnership(line4,False); 

            SetTitle(titlePt)

        canvas.cd();

        txt = TPaveText(0.,0.85,1.,0.95,"NDC");
        txt.SetFillColor(kWhite);
        txt.SetFillStyle(0);
        txt.SetBorderSize(0);
        txt.SetTextAlign(22);#middle,left
        txt.SetTextFont(42);#helvetica
        txt.SetTextSize(0.045);
        txt.AddText("Same and (scaled) mixed event");
        txt.Draw();
        ROOT.SetOwnership(txt,False);       

        nPixels        = 13;
        textHeight     = 0.08;
        startTextX     = 0.10;
        columnsLegend     = 1;
        widthLegend    = 1./numberColumnsPlot;
        heightLegend   = 1./numberRowsPlot;
        marginWidthLeg = 0.15;
        exampleBin        = numberColumnsPlot+fStartBinPtRange-1;
        if numberColumnsPlot > 7:
            startTextX          = 0.05;
            nPixels             = 12;
            widthLegend         = 2./numberColumnsPlot;
            marginWidthLeg      = 0.25;

        # plotting Legend
        padLegend                = TPad("dummyPad","",1-widthLegend,1-heightLegend -0.05,1., 0.85 ,0);   # gives the size of the histo areas
        PadSettings(padLegend, 0, 0, 0, 0);
        padLegend.Draw();
        padLegend.cd();

        textAlice = "";
        if fPlottingType.CompareTo("wip")==0:
            textAlice       = "ALICE work in progress";
        elif fPlottingType.CompareTo("thesis")==0:
            textAlice       = "ALICE this thesis";
        elif fPlottingType.CompareTo("performance")==0:
            textAlice       = "ALICE performance";
        else:
            textAlice       = "ALICE";

        # textEvents="";
        # if fMonteCarloInfo:
        #     textEvents          = "MC";
        # else:
        #     textEvents          = "Data";

        if  padLegend.XtoPixel(padLegend.GetX2()) < padLegend.YtoPixel(padLegend.GetY1()):
            textHeight          = nPixels/padLegend.XtoPixel(padLegend.GetX2()) ;
        else:
            textHeight          = nPixels/padLegend.YtoPixel(padLegend.GetY1());
        
        startTextY     = 0.7;
        differenceText = textHeight*1.20;
        # plot labels
        alice           = TPaveText(startTextX, startTextY, startTextX + 2*textHeight, startTextY - 2*differenceText+0.2, "NDC")
        alice.AddText(textAlice)
        # TLatex(startTextX, startTextY, textAlice);
        SetStyleTLatex( alice, textHeight*3, 1, 1, 42, True);
        alice.Draw();
        # self.PlotLabelsInvMassInPtPlots ( startTextX, startTextY, textHeight, differenceText, textAlice, dateDummy, fEnergy, fDecayChannel, fDetectionChannel);
        # latexDate       = TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy);
        latexPeriod    = TPaveText(startTextX, (startTextY-4*differenceText), startTextX + 2*textHeight, (startTextY-3*differenceText)+0.2, "NDC")
        latexPeriod.AddText(Period)
        SetStyleTLatex( latexPeriod, textHeight*3, 1, 1, 42, True);
        latexPeriod.Draw();
        # energy          = TLatex(startTextX, (startTextY-2.25*differenceText), fEnergy);
        energy       = TPaveText(startTextX, (startTextY-8*differenceText), startTextX + 2*textHeight, (startTextY-4*differenceText)+0.2, "NDC")
        energy.AddText(fEnergy)
        SetStyleTLatex( energy, textHeight*3, 1, 1, 42, True);
        energy.Draw();

        events = TPaveText(startTextX, (startTextY-12*differenceText), startTextX + 2*textHeight, (startTextY-5*differenceText)+0.2, "NDC")
        events.AddText(fMonteCarloInfo)
        SetStyleTLatex( events, textHeight*3, 1, 1, 42, True);
        events.Draw();
        # # process         = TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
        # process     = TPaveText(startTextX, (startTextY-3.25*differenceText), startTextX + textHeight, (startTextY-3.25*differenceText)+0.2, "NDC")
        # process.AddText(fDecayChannel)
        # self.SetStyleTLatex( process, textHeight*1, 1, 1, 42, True);
        # process.Draw();
        # # detprocess      = TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionChannel);
        # detprocess  = TPaveText(startTextX, (startTextY-4.25*differenceText), startTextX + textHeight, (startTextY-4.25*differenceText)+0.2, "NDC")
        # detprocess.AddText(fDetectionChannel)
        # self.SetStyleTLatex( detprocess, textHeight*1, 1, 1, 42, True);
        # detprocess.Draw();
        #        
        legendData     = self.GetAndSetLegend2(  startTextX, startTextY-8*differenceText, 0.8,  startTextY-20*differenceText -0.02, 3*nPixels, columnsLegend, TString(""), 43, marginWidthLeg);
        markersize       = fHistoSamePtBinPlot[exampleBin].GetMarkerSize();
        fHistoSamePtBinPlot[exampleBin].SetMarkerSize(markersize);
        legendData.AddEntry(fHistoSamePtBinPlot[exampleBin],"same evt. #it{{M}}{}".format(decayChannel),"ep");
        markersize       = fHistoSamePtBinPlot[exampleBin].GetMarkerSize();
        fHistoMixedPtBinPlot[exampleBin].SetMarkerSize(markersize);
        legendData.AddEntry(fHistoMixedPtBinPlot[exampleBin], "mixed evt. #it{{M}}{} (scaled)".format(decayChannel), "ep");
        # legendData.AddEntry(line2, "fitting range", "l")
        legendData.AddEntry(line4, "integration rang for scaling", "l")
        legendData.Draw();

        canvas.SaveAs(namePlot.Data());
        ROOT.SetOwnership(canvas, False)
        del padLegend;
        del pad;
        del canvas;

# if __name__ == "__main__":
#     cuts = True
#     generated = True
#     cutname = "qc"
#     period = "LHC23d1k";
#     suffix = "AnyTrack";
#     file = "/Users/alicamarieenderich/20231015_invariant_mass_plots/pi0_data_ptspectrum_pp_13.6TeV_LHC22qAnyTrack_copied_run.root"
#     dirname = ""
#     draw = PlotInvMass("pi0", file, dirname)
#     plottingRange = [0., 0.3]
