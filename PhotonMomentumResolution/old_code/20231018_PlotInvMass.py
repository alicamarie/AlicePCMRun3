import numpy as np
import datetime
import ROOT
from  old_code.file_manager import FileManager
from ROOT import TFile, TDirectory, THashList, TH1F, TH2F, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TGaxis, TPad, TLatex, TBox
from ROOT import gStyle, gROOT, gSystem
from ROOT import TFile, TDirectory, THashList, TF1, TH1F, TCanvas, TH1
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan
from histo_manager import slice_histogram, rebin_histogram
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
    
    def make_common_style(self, g1,marker,size,color,width=1,fill=0):
        g1.SetMarkerStyle(marker);
        g1.SetMarkerColor(color);
        g1.SetMarkerSize(size);
        g1.SetLineColor(color);
        g1.SetLineWidth(width);
        g1.SetFillColor(color);
        g1.SetFillStyle(fill);

    def CanvasSettings(self, c1, leftMargin, rightMargin, topMargin, bottomMargin):
            c1.SetTicks(1,1);
            c1.SetGrid(0,0);
            c1.SetLogy(0);
            c1.SetMargin(leftMargin, rightMargin, topMargin, bottomMargin);
            c1.SetFillColor(0);

    def PadSettings(self, pad1, leftMargin, rightMargin, topMargin, bottomMargin):
        pad1.SetFillColor(0);
        pad1.GetFrame().SetFillColor(0);
        pad1.SetMargin(leftMargin, rightMargin, topMargin, bottomMargin);
        pad1.SetTicks(1,1);

    def DrawHisto(self, plot1, histo1,histofit, Title, XTitle, YTitle, xMin, xMax,  markerSize = 0.2):
        print(type(histo1))
        histo1.GetXaxis().SetRangeUser(xMin, xMax);
        yMin = 0.;
        yMax = 0.;
        if isinstance(histo1, TH1):
            ny = histo1.GetNbinsY() + 1

            for binx in range ( histo1.GetXaxis().FindBin(xMin), histo1.GetXaxis().FindBin(xMax)):
                for biny in range(1, ny):
                    bin_content = histo1.GetBinContent(binx, biny);
                    if bin_content < yMin:
                        yMin = bin_content;
                    if bin_content > yMax:
                        yMax = bin_content;
        else:
            yMax = histo1.GetMaximum(0.0, 0.3);
            yMin = histo1.GetMinimum(0.0, 0.3);
        
        frame1 = plot1.DrawFrame(0., yMin, 0.3, 1.6*yMax);
        frame1.GetXaxis().SetTitle(XTitle);
        frame1.GetYaxis().SetTitle(YTitle);
        frame1.GetXaxis().SetTitleSize(0.05);
        frame1.GetYaxis().SetTitleSize(0.05);
        frame1.GetXaxis().SetTitleOffset(1.0);
        frame1.GetYaxis().SetTitleOffset(1.3);
        frame1.GetXaxis().SetLabelSize(0.05);
        frame1.GetYaxis().SetLabelSize(0.05);
        frame1.GetYaxis().SetMaxDigits(3);
        frame1.GetXaxis().SetLabelOffset(0.01);
        frame1.GetYaxis().SetLabelOffset(0.01);
        ROOT.SetOwnership(frame1,False);
        histo1.SetLineWidth(1);
        frame1.SetMarkerSize(markerSize);
        self.make_common_style(histo1, 1, 1.0, kBlue+1, 1, 0)
        self.make_common_style(histofit, 1, 1.0, kRed+1, 1, 0)        
        histo1.Draw("same")
        histofit.Draw("same")

        HistoTitle = TPaveText(0.,0., 0.5,0.5,"NDC");
        HistoTitle.AddText("haaallo");
        # HistoTitle.SetFillColor(kWhite);
        # HistoTitle.SetFillStyle(0);
        # HistoTitle.SetBorderSize(0);
        # HistoTitle.SetTextAlign(12);#middle,left
        # HistoTitle.SetTextFont(42);#helvetica
        # HistoTitle.SetTextSize(0.5);
        HistoTitle.Draw("");

    def AddLabels(self, text, textSize, textFont = 42, align = 11, textX = 0.1, textY = 0.1):
        txt = TPaveText(textX, textY, textX+0.5, textY+0.5, text);
        # txt.SetFillColor(kWhite);
        # txt.SetFillStyle(0);
        txt.SetBorderSize(0);
        txt.SetTextAlign(align);#middle,left
        txt.SetTextFont(textFont);#helvetica
        txt.SetTextSize(textSize);
        txt. AddText("pp at #sqrt{#it{s}} = 13.6 TeV,|#it{#eta}_{#gamma}| < 0.9 ");
        
    def PlotInvMass(self, fHistoMappingGGInvMassPtBinPlot,fHistoMappingBackNormInvMassPtBinPlot,
                            outputName, PlottingRange, numberRowsPlot, numberColumnsPlot):

            TGaxis.SetMaxDigits(3);
            npt = len(self.arr_pt);
            c1 = TCanvas("c1", "", 1500, 1000)
            c1.SetTicks(1,1);
            self.CanvasSettings(c1, 0, 0, 0, 0);
            c1.cd()
        
            p1 = TPad("p1", "", -0.0, 0.0, 1.0, 1.0, 0)
            self.PadSettings(p1, 0, 0, 0, 0);
            p1.Divide(numberColumnsPlot, numberRowsPlot, 0.0, 0.0)
            p1.Draw()

            place        = 0;
            legendPlace  = [numberColumnsPlot, numberRowsPlot]; # right, top
                
            for iPt in range(npt-1):
                startPt = self.arr_pt[iPt];
                endPt = self.arr_pt[iPt+1];
                place += 1;

                plot1 = p1.cd(place);
                plot1.SetTopMargin(0.15);
                plot1.SetBottomMargin(0.15);
                plot1.SetRightMargin(0.15);
                plot1.SetLeftMargin(0.15); # -> change values later to make neater

                titlePt = "{:.2f} GeV/c < pT < {:.2f} GeV/c".format(startPt, endPt);

                self.DrawHisto(plot1, fHistoMappingGGInvMassPtBinPlot[iPt],fHistoMappingBackNormInvMassPtBinPlot[iPt], titlePt,
                            "M_{\gamma \gamma} (GeV/c^2)", "dN_{\gamma \gamma}/dM_{\gamma \gamma}",
                            PlottingRange[0], PlottingRange[1]);

                # self.DrawHisto(plot1, fHistoMappingBackNormInvMassPtBinPlot[iPt], titlePt,
                #             "M_{\gamma \gamma} (GeV/c^2)", "dN_{\gamma \gamma}/dM_{\gamma \gamma}",
                #             PlottingRange[0], PlottingRange[1]);


            
            c1.cd();
            widthLegend = 1./numberColumnsPlot;
            heightLegend = 1./numberRowsPlot;
            exampleBin = numberColumnsPlot;

            # plotting Legend
            LegendPad = TPad("pad","",1-widthLegend,1-heightLegend,1.,1.,0);   # gives the size of the histo areas
            self.PadSettings( LegendPad, 0, 0, 0, 0);
            LegendPad.Draw();
            LegendPad.cd();

            textY = 0.1;
            textX = 0.1;
            textSize = 0.1;
            text = "pp at #sqrt{#it{s}} = 13.6 TeV";
            self.AddLabels(text, textSize,42, 11, textX, textY)

            # legendData = TLegend(textX, textY, textX+0.5, textY+0.5 )#, "", 43)
            # legendData.AddEntry(fHistoMappingGGInvMassPtBinPlot[exampleBin], "same evt. M_{\gamma \gamma} (BG+Signal)")
            # legendData.AddEntry(fHistoMappingBackNormInvMassPtBinPlot[exampleBin], "#it{M}_{\gamma \gamma}");
            # legendData.Draw();

            c1.SaveAs(outputName);
            del LegendPad;
            # del legendData
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
    plottingRange = [0., 0.3]

   # draw.PlotInvMass(histogram, normalizedHistogram, "InvMass.pdf", plottingRange, 5, 5)
# access the directory to get to the actual histograms
