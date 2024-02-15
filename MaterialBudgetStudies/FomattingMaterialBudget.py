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

def FrameSettings(frame1):

    frame1.GetXaxis().SetTitleSize(0.045);
    frame1.GetYaxis().SetTitleSize(0.045);
    frame1.GetXaxis().SetTitleOffset(1.0);
    frame1.GetYaxis().SetTitleOffset(1.52);
    frame1.GetXaxis().SetLabelSize(0.045);
    frame1.GetYaxis().SetLabelSize(0.045);
    frame1.GetYaxis().SetMaxDigits(3);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetNdivisions(507, True);
    frame1.GetYaxis().SetNdivisions(5, 5, 0, kTRUE)

def ALICEtext(thesis):
    txt = TPaveText(0.90,0.67,0.95,0.82,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(32);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.045);
    if thesis == "thesis":
        txt.AddText("ALICE");
        txt.AddText("this thesis")
    if thesis == "simulation":
        txt.AddText("ALICE simulation")
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.Draw();
    ROOT.SetOwnership(txt,False);

def RatioLegendSettings():
        leg = TLegend(0.17,0.82,0.35,0.97);
        leg.SetBorderSize(0);
        leg.SetFillColor(kWhite);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.08);
        return leg      

def FrameSettingsRatio(frame2):
        frame2.GetXaxis().SetTitleSize(0.10);
        frame2.GetYaxis().SetTitleSize(0.10);
        frame2.GetXaxis().SetTitleOffset(1.0);
        frame2.GetYaxis().SetTitleOffset(0.7);
        frame2.GetXaxis().SetLabelSize(0.10);
        frame2.GetYaxis().SetLabelSize(0.10);
        frame2.GetYaxis().CenterTitle(True);
        frame2.GetXaxis().SetLabelOffset(0.01);
        frame2.GetYaxis().SetLabelOffset(0.01);
        frame2.GetYaxis().CenterTitle(True);
        frame2.GetYaxis().SetNdivisions(5, 5, 0, kTRUE)
        frame2.GetXaxis().SetNdivisions(507, True);
        ROOT.SetOwnership(frame2,False);    

def FrameSettings2Dtwoplots(frame1):
        frame1.GetXaxis().SetTitleSize(0.05);
        frame1.GetYaxis().SetTitleSize(0.05);
        frame1.GetXaxis().SetTitleOffset(1.0);
        frame1.GetYaxis().SetTitleOffset(1.0);
        frame1.GetXaxis().SetLabelSize(0.05);
        frame1.GetYaxis().SetLabelSize(0.05);
        frame1.GetYaxis().SetMaxDigits(3);
        frame1.GetXaxis().SetLabelOffset(0.01);
        frame1.GetYaxis().SetLabelOffset(0.01);
        frame1.GetXaxis().SetNdivisions(507, True);
        frame1.GetYaxis().SetNdivisions(5, 5, 0, kTRUE)
        ROOT.SetOwnership(frame1,False);

def FrameSettings2D(frame1):
        frame1.GetXaxis().SetTitleSize(0.04);
        frame1.GetYaxis().SetTitleSize(0.04);
        frame1.GetXaxis().SetTitleOffset(1.7);
        frame1.GetYaxis().SetTitleOffset(1.7);
        frame1.GetXaxis().SetLabelSize(0.04);
        frame1.GetYaxis().SetLabelSize(0.04);
        frame1.GetYaxis().SetMaxDigits(3);
        frame1.GetXaxis().SetLabelOffset(0.01);
        frame1.GetYaxis().SetLabelOffset(0.01);
        frame1.GetXaxis().SetNdivisions(507, True);
        frame1.GetYaxis().SetNdivisions(5, 5, 0, kTRUE)
        ROOT.SetOwnership(frame1,False);

def ALICEtext2Dtwoplots(thesis,oneortwo):
        if oneortwo ==1:
                height = 0.66
        else:
              height = 0.81
        txt = TPaveText(0.8,height,0.85,height+0.12,"NDC"); #needed to change values because two plots in canvas
        txt.SetFillColor(kWhite);
        txt.SetFillStyle(0);
        txt.SetBorderSize(0);
        txt.SetTextAlign(32);#middle,left
        txt.SetTextFont(42);#helvetica
        txt.SetTextSize(0.04);
        if thesis == "thesis":
                txt.AddText("ALICE");
                txt.AddText("this thesis")
        if thesis == "simulation":
                txt.AddText("ALICE simulation")
        txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
        txt.Draw();
        ROOT.SetOwnership(txt,False);

def ALICEtext2D(thesis):
    txt = TPaveText(0.72,0.77,0.78,0.85,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(32);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.03);
    if thesis == "thesis":
        txt.AddText("ALICE");
        txt.AddText("this thesis")
    if thesis == "simulation":
        txt.AddText("ALICE simulation")
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.Draw();
    ROOT.SetOwnership(txt,False);



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