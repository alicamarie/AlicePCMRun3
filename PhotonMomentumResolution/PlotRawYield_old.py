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

        # self.make_common_style(histo1, kFullCircle, markerSize, markerColor, 1)
        histo1.SetMarkerStyle(kFullCircle);
        histo1.SetMarkerColor(markerColor);
        histo1.SetMarkerSize(markerSize);
        histo1.SetLineColor(markerColor);
        histo1.SetLineWidth(1);
        histo1.SetFillColor(markerColor);
        histo1.SetFillStyle(0);
        histo1.DrawCopy(drawsettings);
        # histo1.SetLineStyle(1);
        # histo1.SetLineColor(4);
        # histo1.SetMarkerColor(4);
        # histo1.SetMarkerStyle(24);
        # histo1.SetLineWidth(1);
        #histo1.DrawCopy(drawsettings);

        # TitlePlot = TPaveText(0.1, 0.95, 0.2, 1.0, "NDC")
        # TitlePlot.AddText("{}".format(Title))
        # # TitlePlot = TLatex(0.1, 0.95, "{}".format(Title));
        # # TitlePlot.SetNDC();
        # TitlePlot.SetTextColor(1);
        # TitlePlot.SetTextSize(0.062);
        # TitlePlot.Draw();
    

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
        # TLatex(startTextX, startTextY, textAlice);
        self.SetStyleTLatex( alice, textHeight*1.3, 1, 1, 42, True, textAlign);
        alice.Draw();
        # latexDate       = TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy);
        latexDate    = TPaveText(startTextX, startTextY, startTextX + textHeight, startTextY+0.2, "NDC")
        latexDate.AddText(dateDummy)
        self.SetStyleTLatex( latexDate, textHeight, 1, 1, 42, True, textAlign);
        latexDate.Draw();
        # energy          = TLatex(startTextX, (startTextY-2.25*differenceText), fEnergy);
        energy       = TPaveText(startTextX, startTextY, startTextX + textHeight, startTextY+0.2, "NDC")
        energy.AddText(fEnergy)
        self.SetStyleTLatex( energy, textHeight*1, 1, 1, 42, True, textAlign);
        energy.Draw();
        # process         = TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
        process     = TPaveText(startTextX, startTextY, startTextX + textHeight, startTextY+0.2, "NDC")
        process.AddText(fEnergy)
        self.SetStyleTLatex( process, textHeight*1, 1, 1, 42, True, textAlign);
        process.Draw();
        # detprocess      = TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionChannel);
        detprocess  = TPaveText(startTextX, startTextY, startTextX + textHeight, startTextY+0.2, "NDC")
        detprocess.AddText(fDetectionChannel)
        self.SetStyleTLatex( detprocess, textHeight*1, 1, 1, 42, True, textAlign);
        detprocess.Draw();
        # events          = TLatex();
        # if textEvents.CompareTo("") != 0 and fNEvents > 0:
        #     events              = TLatex(startTextX, (startTextY-5.25*differenceText), "{}: {:.1e} events".format(textEvents, fNEvents));
        #     self.SetStyleTLatex( events, textHeight*1, 1, 1, 42, True, textAlign);
        #     events.Draw();
    
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
        # canvas = TCanvas("c1", "", 3000, 2000)
        # canvas.SetTicks(1,1);
        # self.CanvasSettings(canvas, 0.15, 0.15, 0.15, 0.15);
        # canvas.cd()
    
        # pad = TPad("p1", "", -0.0, 0.0, 1.0, 1.0, 0)

        # pad.SetTopMargin(0.15);
        # pad.SetBottomMargin(0.15);
        # pad.SetRightMargin(0.15);
        # pad.SetLeftMargin(0.15); # -> change values later to make neater

        titlePt = fHistoParameter[0].GetTitle(); # ACCESS TITLE

        yMin_, yMax_ = 1e-8, 1.4*1e-3
        color = [kRed+1, kBlue+1, kGreen+2, kMagenta+2, kCyan+1]
        # pad.SetLogy();
        # pad.DrawFrame(0., yMin_, 12., yMax_);
        
        # for icut in range(len(fHistoParameter)):
        #     self.DrawHisto(pad,fHistoParameter[icut], titlePt,
        #                 "M_{} (GeV/c^2)".format(decayChannel), "dN_{}/dM_{}".format(decayChannel, decayChannel), color[icut], "E1,same");

        c1 = TCanvas("c0","c0",0,0,800,800);
        c1.Divide(1,1,1e-3,1e-3);
        p1 = c1.cd(1);
        p1.SetPad(0,0,1,1);
        p1.SetMargin(0.17,0.05,0.1,0.15);
        p1.SetTicks(1,1);
        p1.SetLogy();

        frame1 = p1.DrawFrame(0, yMin_, 12, yMax_);
        frame1.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})");
        frame1.GetYaxis().SetTitle("#frac{d#it{N}}{d#it{p}_{T}} (GeV/#it{c})^{-1}");        
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

        for icut in range(len(fHistoParameter)):
            #fHistoParameter[icut].SetMarkersize(1)
            self.DrawHisto(frame1,fHistoParameter[icut], titlePt,
                        "M_{} (GeV/c^2)".format(decayChannel), "dN_{}/dM_{}".format(decayChannel, decayChannel), color[icut], "E1,same");

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
        
        # txt = TPaveText(0.90,0.77,0.95,0.82,"NDC");
        # txt.SetFillColor(kWhite);
        # txt.SetFillStyle(0);
        # txt.SetBorderSize(0);
        # txt.SetTextAlign(32);#middle,left
        # txt.SetTextFont(42);#helvetica
        # txt.SetTextSize(0.02);
        # txt.AddText("ALICE this thesis");
        # txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
        # txt.Draw();
        # ROOT.SetOwnership(txt,False);

        # leg = TLegend(0.17,0.72,0.35,0.82);
        # leg.SetBorderSize(0);
        # leg.SetFillColor(kWhite);
        # leg.SetFillStyle(0);
        # leg.SetTextSize(0.045);
        # leg.AddEntry(h1_mc_rec ,"M.C. rec. primary #gamma (LHC23d1k)","LP");
        # leg.AddEntry(h1_data   ,"Data #gamma candidates (LHC22f pass4)","LP");
        # leg.Draw("");
        # ROOT.SetOwnership(leg,False);


        # c1.cd();
        # nPixels        = 13;
        # textHeight     = 0.08;
        # startTextX     = 0.10;
        # columnsLegend     = 1;
        # widthLegend    = 1./numberColumnsPlot;
        # heightLegend   = 1./numberRowsPlot;
        # marginWidthLeg = 0.15;
        # exampleBin        = 1 #numberColumnsPlot-1;
        # if numberColumnsPlot > 7:
        #     startTextX          = 0.05;
        #     nPixels             = 12;
        #     widthLegend         = 2./numberColumnsPlot;
        #     marginWidthLeg      = 0.25;

        # # plotting Legend
        # padLegend                = TPad("dummyPad","",1-widthLegend,1-heightLegend,1.,1.,0);   # gives the size of the histo areas
        # self.PadSettings(padLegend, 0, 0, 0, 0);
        # padLegend.Draw();
        # padLegend.cd();

        # textAlice = "";
        # if fPlottingType.CompareTo("wip")==0:
        #     textAlice       = "ALICE work in progress";
        # elif fPlottingType.CompareTo("thesis")==0:
        #     textAlice       = "ALICE this thesis";
        # elif fPlottingType.CompareTo("performance")==0:
        #     textAlice       = "ALICE performance";
        # else:
        #     textAlice       = "ALICE";

        # textEvents="";
        # if fMonteCarloInfo:
        #     textEvents          = "MC";
        # else:
        #     textEvents          = "Data";
    

        # if  padLegend.XtoPixel(padLegend.GetX2()) < padLegend.YtoPixel(padLegend.GetY1()):
        #     textHeight          = nPixels/padLegend.XtoPixel(padLegend.GetX2()) ;
        # else:
        #     textHeight          = nPixels/padLegend.YtoPixel(padLegend.GetY1());
        
        # startTextY     = 0.7;
        # differenceText = textHeight*1.05;
        # # plot labels
        # alice           = TPaveText(startTextX, startTextY, startTextX + 2*textHeight, startTextY - 2*differenceText+0.2, "NDC")
        # alice.AddText(textAlice)
        # # TLatex(startTextX, startTextY, textAlice);
        # self.SetStyleTLatex( alice, textHeight*2, 1, 1, 42, True);
        # alice.Draw();
        # # self.PlotLabelsInvMassInPtPlots ( startTextX, startTextY, textHeight, differenceText, textAlice, dateDummy, fEnergy, fDecayChannel, fDetectionChannel);
        # # latexDate       = TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy);
        # latexPeriod    = TPaveText(startTextX, (startTextY-2*differenceText), startTextX + 2*textHeight, (startTextY-3*differenceText)+0.2, "NDC")
        # latexPeriod.AddText(Period)
        # self.SetStyleTLatex( latexPeriod, textHeight*2, 1, 1, 42, True);
        # latexPeriod.Draw();
        # # energy          = TLatex(startTextX, (startTextY-2.25*differenceText), fEnergy);
        # energy       = TPaveText(startTextX, (startTextY-4*differenceText), startTextX + 2*textHeight, (startTextY-4*differenceText)+0.2, "NDC")
        # energy.AddText(fEnergy)
        # self.SetStyleTLatex( energy, textHeight*2, 1, 1, 42, True);
        # energy.Draw();
        # # # process         = TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
        # # process     = TPaveText(startTextX, (startTextY-3.25*differenceText), startTextX + textHeight, (startTextY-3.25*differenceText)+0.2, "NDC")
        # # process.AddText(fDecayChannel)
        # # self.SetStyleTLatex( process, textHeight*1, 1, 1, 42, True);
        # # process.Draw();
        # # # detprocess      = TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionChannel);
        # # detprocess  = TPaveText(startTextX, (startTextY-4.25*differenceText), startTextX + textHeight, (startTextY-4.25*differenceText)+0.2, "NDC")
        # # detprocess.AddText(fDetectionChannel)
        # # self.SetStyleTLatex( detprocess, textHeight*1, 1, 1, 42, True);
        # # detprocess.Draw();

        # color = [kRed, kBlue, kGreen, kMagenta, kCyan]   
        # legendData     = self.GetAndSetLegend2(  startTextX, startTextY-4*differenceText, 0.85,  startTextY-10*differenceText, 2*nPixels, columnsLegend, TString(""), 43, marginWidthLeg);
        # markersize       = fHistoParameter[0].GetMarkerSize();
        # for icut in range(len(fHistoParameter)):
        #     fHistoParameter[icut].SetMarkerSize(markersize);
        #     legendData.AddEntry(fHistoParameter[icut],"fitparameters for {}".format(cutnames[icut]),"ep");
        
        # linesize         = fHistoParameter[0].GetLineWidth();
        # fHistoParameter[0].SetLineWidth(linesize);
        # legendData.Draw();


        txt = TPaveText(0.87,0.76,0.92,0.82,"NDC");
        txt.SetFillColor(kWhite);
        txt.SetFillStyle(0);
        txt.SetBorderSize(0);
        txt.SetTextAlign(32);#middle,left
        txt.SetTextFont(42);#helvetica
        txt.SetTextSize(0.02);
        txt.AddText("ALICE this thesis");
        txt.AddText(Period)
        txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
        txt.Draw();
        ROOT.SetOwnership(txt,False);

        leg = TLegend(0.2,0.72,0.35,0.82);
        leg.SetBorderSize(0);
        leg.SetFillColor(kWhite);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.02);
        markersize       = fHistoParameter[0].GetMarkerSize();
        for icut in range(len(fHistoParameter)):
            fHistoParameter[icut].SetMarkerSize(markersize);
            leg.AddEntry(fHistoParameter[icut],"fitparameters for {}".format(cutnames[icut]),"ep");
        leg.Draw("");
        ROOT.SetOwnership(leg,False);

        # date = datetime.date.today().strftime("%Y%m%d");
        c1.Modified();
        c1.Update();
        ROOT.SetOwnership(c1,False);
        c1.SaveAs(namePlot.Data());
        c1.Close();

        # del padLegend;
        # del pad;
        del c1;

# if __name__ == "__main__":
#     cuts = True
#     generated = True
#     cutname = "qc"
#     period = "LHC23d1k";
#     suffix = "AnyTrack";
#     file = "/Users/alicamarieenderich/20231015_invariant_mass_plots/pi0_data_ptspectrum_pp_13.6TeV_LHC22qAnyTrack_copied_run.root"
#     dirname = ""
#     draw = PlotHistoInvMass("pi0", file, dirname)
#     plottingRange = [0., 0.3]

   # draw.PlotInvMass(histogram, normalizedHistogram, "InvMass.pdf", plottingRange, 5, 5)
# access the directory to get to the actual histograms