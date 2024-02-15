import ROOT
from ROOT import TH1D, TH2D, TH3D, TGraph, TGraphErrors, TLatex, TPaveText, TLegend
from ROOT import gROOT, gSystem, gStyle, gPad
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kOrange, kFullCircle
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kFullCircle, kTRUE

def CanvasSettings(c1, leftMargin, rightMargin, topMargin, bottomMargin):
        c1.SetTickx();
        c1.SetTicky();
        c1.SetLogy(0);
        c1.SetMargin(leftMargin, rightMargin, topMargin, bottomMargin);
        c1.SetFillColor(0);

def PadSettings( pad1, leftMargin, rightMargin, topMargin, bottomMargin):
    pad1.SetFillColor(0);
    pad1.GetFrame().SetFillColor(0);
    pad1.SetMargin(leftMargin, rightMargin, topMargin, bottomMargin);
    pad1.SetTickx();
    pad1.SetTicky();

def FrameSettings(pad1, XTitle, YTitle, xMin, xMax):

        pad1.GetYaxis().SetLabelSize(0.06);
        pad1.GetYaxis().SetTitleSize(0.06);
        pad1.GetYaxis().SetTitleOffset(1.0);
        pad1.GetYaxis().SetTitle(YTitle);
        pad1.GetYaxis().SetDecimals();

        pad1.GetXaxis().SetRangeUser(xMin, xMax);
        pad1.GetXaxis().SetLabelSize(0.06);
        pad1.GetXaxis().SetTitleSize(0.06);
        pad1.GetXaxis().SetTitleOffset(1.0);
        pad1.GetXaxis().SetTitle(XTitle);
        pad1.GetXaxis().SetNdivisions(507, True);
        pad1.GetYaxis().SetNdivisions(5, 5, 0, kTRUE)

def FrameSettingsCombined(pad1, XTitle, YTitle, xMin, xMax):

        pad1.GetYaxis().SetLabelSize(0.06);
        pad1.GetYaxis().SetTitleSize(0.06);
        pad1.GetYaxis().SetTitleOffset(1.0);
        pad1.GetYaxis().SetTitle(YTitle);
        pad1.GetYaxis().SetDecimals();

        pad1.GetXaxis().SetRangeUser(xMin, xMax);
        pad1.GetXaxis().SetLabelSize(0.06);
        pad1.GetXaxis().SetTitleSize(0.06);
        pad1.GetXaxis().SetTitleOffset(1.0);
        pad1.GetXaxis().SetTitle(XTitle);
        pad1.GetXaxis().SetNdivisions(507, True);
        pad1.GetYaxis().SetNdivisions(5, 5, 0, kTRUE)
    
def DrawHisto(histo1, markerColor, drawsettings, markerSize = 1.8, markerstyle = kFullCircle):
        histo1.SetMarkerStyle(markerstyle);
        histo1.SetMarkerColor(markerColor);
        histo1.SetMarkerSize(markerSize);
        histo1.SetLineColor(markerColor);
        histo1.SetLineWidth(1);
        histo1.SetFillColor(markerColor);
        histo1.SetFillStyle(0);
        histo1.DrawCopy(drawsettings);
            
def DrawHistoCombined(histo1, markerColor, drawsettings, markerSize = 2.2):
        histo1.SetMarkerStyle(kFullCircle);
        histo1.SetMarkerColor(markerColor);
        histo1.SetMarkerSize(markerSize);
        histo1.SetLineColorAlpha(markerColor, 0.6);
        histo1.SetLineWidth(1);
        histo1.SetFillColor(markerColor);
        histo1.SetFillStyle(0);
        histo1.DrawCopy(drawsettings);
            
def SetTitle(titlePt):       
    TitlePlot = TPaveText(0.3, 0.92, 0.6, 0.98, "NDC")
    TitlePlot.AddText("{}".format(titlePt))
    TitlePlot.SetTextColor(1);
    TitlePlot.SetTextSize(0.08);
    TitlePlot.SetFillStyle(0)
    TitlePlot.SetBorderSize(0)
    TitlePlot.Draw();
    ROOT.SetOwnership(TitlePlot, False);

def SetStyleTLatex(text, textSize, lineWidth, textColor = 1, textFont = 42, kNDC = True, align = 11):
        text.SetTextFont(textFont);
        text.SetTextColor(textColor);
        text.SetTextSize(textSize);
        text.SetLineWidth(lineWidth);
        text.SetTextAlign(align);   
        text.SetFillStyle(0)
        text.SetBorderSize(0)