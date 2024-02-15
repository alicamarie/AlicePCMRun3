import ROOT
from ROOT import TH1D, TH2D, TH3D, TGraph, TGraphErrors, TGraphAsymmErrors, TLatex
from ROOT import gROOT, gSystem, gStyle, gPad
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kOrange

def make_common_style(g1,marker,size,color,width=1,fill=0):
    g1.SetMarkerStyle(marker);
    g1.SetMarkerColor(color);
    g1.SetMarkerSize(size);
    g1.SetLineColor(color);
    g1.SetLineWidth(width);
    g1.SetFillColor(color);
    g1.SetFillStyle(fill);

def CanvasSettings(     c1,
                        leftMargin,
                        rightMargin,
                        topMargin,
                        bottomMargin):
        c1.SetTicks(1,1);
        c1.SetGrid(0,0);
        c1.SetLogy(0);
        c1.SetMargin(leftMargin, rightMargin, topMargin, bottomMargin);
        c1.SetFillColor(0);

def PadSettings(        pad1,
                        leftMargin,
                        rightMargin,
                        topMargin,
                        bottomMargin):
    pad1.SetFillColor(0);
    pad1.GetFrame().SetFillColor(0);
    pad1.SetMargin(leftMargin, rightMargin, topMargin, bottomMargin);
    pad1.SetTicks(1,1);


def DrawHisto(             histo1,
                            Title,
                            XTitle,
                            YTitle,
                            xMin,
                            xMax,
                            markerSize = 0.2
                        ):

            histo1.GetXaxis().SetRangeUser(xMin, xMax);
            yMin = 0;
            yMax = 0;
            for i in range ( histo1.GetXaxis().FindBin(xMin), histo1.GetXaxis().FindBin(xMax)):
                if histo1.GetBinContent(i) < yMin:
                    yMin = histo1.GetBinContent(i);
                
                if histo1.GetBinContent(i) > yMax:
                    yMax = histo1.GetBinContent(i);
                
            histo1.GetYaxis().SetRangeUser(yMin, 1.5*yMax);

            histo1.SetXTitle(XTitle.Data());

            histo1.SetYTitle(YTitle.Data());
            
            histo1.GetYaxis().SetLabelSize(0.02);
            histo1.GetYaxis().SetTitleSize(0.025);
            histo1.GetYaxis().SetDecimals();
            histo1.GetYaxis().SetTitleOffset(0.5);
            histo1.GetXaxis().SetTitleSize(0.025);
            histo1.GetXaxis().SetLabelSize(0.02);
            histo1.SetMarkerStyle(20);
            histo1.SetMarkerColor(1);
            histo1.SetLineColor(1);
            histo1.SetLineWidth(1);
            histo1.SetMarkerSize(markerSize);
            histo1.SetTitleOffset(1.2, "XY");
            histo1.SetTitleSize(0.05, "XY");
            histo1.GetYaxis().SetLabelSize(0.05);
            histo1.GetXaxis().SetLabelSize(0.05);
            histo1.GetXaxis().SetNdivisions(507, True);

            histo1.SetLineStyle(1);
            histo1.SetLineColor(4);
            histo1.SetMarkerColor(4);
            histo1.SetMarkerStyle(24);
            histo1.SetLineWidth(1);
            histo1.DrawCopy("hist,same");

            TitlePlot = TLatex(0.1, 0.95, "{}".format(Title.Data()));
            TitlePlot.SetNDC();
            TitlePlot.SetTextColor(1);
            TitlePlot.SetTextSize(0.062);
            TitlePlot.Draw();