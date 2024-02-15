import numpy as np
import datetime
import ROOT
from  old_code.file_manager import FileManager
from ROOT import TFile, TDirectory, THashList, TH1F, TH2F, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TGaxis, TPad, TLatex, TBox
from ROOT import gStyle, gROOT, gSystem
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan
from histo_manager import slice_histogram, rebin_histogram
#from painter import make_common_style, CanvasSettings, PadSettings, DrawHisto
import re
import numpy as np
import datetime
import math
import ctypes
gStyle.SetOptStat(0);
gStyle.SetOptTitle(0);

class PlotInvMass:
    def __init__(self):
        print("default constructor is called");
    def __init__(self, meson, filename, dirname):
        print("target meson = {0} , filename = {1} , dirname = {2}".format(meson, filename, dirname));
        self.meson = meson;
        self.rootfile = TFile.Open(filename, "READ");
        # self.list_pcm = self.rootfile.Get("PCMPCM")
        # self.list_qc_qc = self.list_pcm.FindObject("qc_qc")
        # self.list_gauss = self.list_qc_qc.FindObject("gausexplinear");
        # self.list_fit = self.list_gauss.FindObject("fit_0.04_0.24_GeVc2");
        # self.rootdir = self.rootfile.Get(dirname);
        # self.list_ev = self.rootdir.Get("Event");
        # self.list_pair = self.rootdir.Get("Pair");
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
        #self.list_pair_ss = self.list_pair.FindObject(ssname);

    def set_cutname(self, cutname):
        self.cutname = cutname;
        # if self.list_ev_ss is None or self.list_pair_ss is None:
        #     print("Please define subsystem name first!");
        #     return None;
        #self.list_ev_ss_cut   = self.list_ev_ss.FindObject(cutname);
        self.list_ss_cut = self.list_ss.FindObject(cutname);
    
    def set_fitname(self, fitname):
        self.fitname = fitname;
        self.list_fitname  = self.list_ss_cut.FindObject(fitname);

    def set_fit_range(self, fit_min, fit_max):
        self.fit_min = fit_min;
        self.fit_max = fit_max;
    
    def set_fitrange(self):
        #self.fitname = fitname;
        self.list_fitrange  = self.list_fitname.FindObject("fit_{0:3.2f}_{1:3.2f}_GeVc2".format(self.fit_min, self.fit_max));

    def make_common_style(g1,marker,size,color,width=1,fill=0):
        g1.SetMarkerStyle(marker);
        g1.SetMarkerColor(color);
        g1.SetMarkerSize(size);
        g1.SetLineColor(color);
        g1.SetLineWidth(width);
        g1.SetFillColor(color);
        g1.SetFillStyle(fill);

    def CanvasSettings(c1, leftMargin, rightMargin, topMargin, bottomMargin):
            c1.SetTicks(1,1);
            c1.SetGrid(0,0);
            c1.SetLogy(0);
            c1.SetMargin(leftMargin, rightMargin, topMargin, bottomMargin);
            c1.SetFillColor(0);

    def PadSettings(pad1, leftMargin, rightMargin, topMargin, bottomMargin):
        pad1.SetFillColor(0);
        pad1.GetFrame().SetFillColor(0);
        pad1.SetMargin(leftMargin, rightMargin, topMargin, bottomMargin);
        pad1.SetTicks(1,1);

    def DrawHisto(self, histo1, Title, XTitle, YTitle, xMin, xMax,  markerSize = 0.2):
        histo1.GetXaxis().SetRangeUser(xMin, xMax);
        yMin = 0;
        yMax = 0;
        for binx in range ( histo1.GetXaxis().FindBin(xMin), histo1.GetXaxis().FindBin(xMax)):
            for biny in range(1, histo1.GetNbinsY() + 1):
                bin_content = histo1.GetBinContent(binx, biny);
                if bin_content < yMin:
                    yMin = bin_content;

        histo1.GetYaxis().SetRangeUser(yMin, 1.6*yMax);
        histo1.SetXTitle(XTitle);
        histo1.SetYTitle(YTitle);
        histo1.GetYaxis().SetLabelSize(0.05);
        histo1.GetYaxis().SetTitleSize(0.025);
        histo1.GetYaxis().SetDecimals();
        histo1.GetYaxis().SetTitleOffset(0.5);
        histo1.GetXaxis().SetTitleSize(0.025);
        histo1.GetXaxis().SetLabelSize(0.05);
        histo1.SetMarkerStyle(20);
        histo1.SetMarkerColor(1);
        histo1.SetLineColor(1);
        histo1.SetLineWidth(1);
        histo1.SetMarkerSize(markerSize);
        histo1.DrawCopy("hist,same");

        HistoTitle = TLatex(0.1, 0.95, "{}".format(Title));
        HistoTitle.SetNDC();
        HistoTitle.SetTextColor(1);
        HistoTitle.SetTextSize(0.062);
        HistoTitle.Draw();

    def AddLabels(self, text, textSize, textFont = 42, align = 11, textX = 0.1, textY = 0.1):
        txt = TLatex(textX, textY, text);
        txt.SetFillColor(kWhite);
        txt.SetFillStyle(0);
        txt.SetBorderSize(0);
        txt.SetTextAlign(align);#middle,left
        txt.SetTextFont(textFont);#helvetica
        txt.SetTextSize(textSize);
        txt. AddText("pp at #sqrt{#it{s}} = 13.6 TeV,|#it{#eta}_{#gamma}| < 0.9 ");
        
    def PlotInvMass(self, fHistoMappingGGInvMassPtBinPlot,fHistoMappingBackNormInvMassPtBinPlot,
                            namePlot, nameCanvas,namePad, fPlottingRangeMeson,
                            numberRowsPlot, numberColumnsPlot, fStartBinPtRange):

            TGaxis.SetMaxDigits(3);
            npt = len(self.arr_pt);
            c1 = TCanvas(nameCanvas, "", 1500, 1000)
            c1.SetTicks(1,1);
            self.CanvasSettings(c1, 0, 0, 0, 0);
            c1.cd()
        
            p1 = TPad(namePad, "", -0.0, 0.0, 1.0, 1.0, 0)
            self.PadSettings(p1, 0, 0, 0, 0);
            p1.Divide(numberColumnsPlot, numberRowsPlot, 0.0, 0.0)
            p1.Draw()

            place        = 0;
            legendPlace  = [numberColumnsPlot, numberRowsPlot]; # right, bottom
                
            for iPt in range(fStartBinPtRange, npt):
                startPt = self.arr_pt[iPt];
                endPt = self.arr_pt[iPt + 1];

                place += 1;

                if place > legendPlace[0] - 1 and place < legendPlace[1] + 1:
                    iPt -= 1;
                else:
                    p1.cd(place);
                    p1.cd(place).SetTopMargin(0.15);
                    p1.cd(place).SetBottomMargin(0.15);
                    p1.cd(place).SetRightMargin(0.15);
                    p1.cd(place).SetLeftMargin(0.15); # -> change values later to make neater

                    titlePt = "{:.2f} GeV/c < pT < {:.2f} GeV/c".format(startPt, endPt);

                    self.DrawHisto(fHistoMappingGGInvMassPtBinPlot[iPt], titlePt,
                                "M_{\gamma \gamma} (GeV/c^2)", "dN_{\gamma \gamma}/dM_{\gamma \gamma}",
                                fPlottingRangeMeson[0], fPlottingRangeMeson[1]);

                    self.DrawHisto(fHistoMappingBackNormInvMassPtBinPlot[iPt], titlePt,
                                "M_{\gamma \gamma} (GeV/c^2)", "dN_{\gamma \gamma}/dM_{\gamma \gamma}",
                                fPlottingRangeMeson[0], fPlottingRangeMeson[1]);

                    box = TBox(self.fit_min, fHistoMappingGGInvMassPtBinPlot[iPt].GetMaximum() * 0.93,
                                    self.fit_max, fHistoMappingGGInvMassPtBinPlot[iPt].GetMaximum() * 0.91);
                    box.SetFillStyle(1001);
                    box.SetFillColor(ROOT.kAzure + 9);
                    box.Draw("same");
            
            c1.cd();
            # columnsLegend     = 1;
            widthLegend    = 1./numberColumnsPlot;
            heightLegend   = 1./numberRowsPlot;
            # marginWidthLeg = 0.15;
            exampleBin        = numberColumnsPlot+fStartBinPtRange-1;
            # if numberColumnsPlot > 7:
            #     widthLegend         = 2./numberColumnsPlot;
            #     marginWidthLeg      = 0.25;

            # plotting Legend
            padLegend                = TPad("dummyPad","",1-widthLegend,1-heightLegend,1.,1.,0);   # gives the size of the histo areas
            self.PadSettings( padLegend, 0, 0, 0, 0);
            padLegend.Draw();
            padLegend.cd();

            textY    = 0.9;
            textX = 0.1;
            textSize = 0.05;
            text = "pp at #sqrt{#it{s}} = 13.6 TeV";
            self.AddLabels(text, textSize,42, 11, textX, textY)

            legendData = ROOT.TLegend(textX, textY, textX+0.5, textY+0.5 , "", 43)
            #legendData.SetNColumns(columnsLegend)
            #legendData.SetMargin(marginWidthLeg)

            markersize = 1.0;
            fHistoMappingGGInvMassPtBinPlot[exampleBin].SetMarkerSize(3 * markersize)
            legendData.AddEntry(fHistoMappingGGInvMassPtBinPlot[exampleBin], "same evt. M_{\gamma \gamma} (BG+Signal)")

            linesize = fHistoMappingBackNormInvMassPtBinPlot[exampleBin].GetLineWidth()
            fHistoMappingBackNormInvMassPtBinPlot[exampleBin].SetLineWidth(5 * linesize)

            legendData.AddEntry(fHistoMappingBackNormInvMassPtBinPlot[exampleBin], "#it{{M}}_{\gamma \gamma}");
            legendData.Draw();

            c1.SaveAs(namePlot);
            del padLegend;
            del p1;
            del c1;

if __name__ == "__main__":
    cuts = True
    generated = True
    cutname = "qc"
    period = "LHC23d1k";
    suffix = "AnyTrack";
    file = "/Users/alicamarieenderich/20231015_invariant_mass_plots/pi0_data_ptspectrum_pp_13.6TeV_LHC22qAnyTrack_copied_run.root"
    dirname = ""
    draw = PlotInvMass("pi0", file, dirname)

    run = draw.PlotInvMass()
# access the directory to get to the actual histograms
