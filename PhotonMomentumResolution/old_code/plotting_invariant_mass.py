import numpy as np
import datetime
import ROOT
from  old_code.file_manager import FileManager
from ROOT import TFile, TDirectory, THashList, TH1F, TH2F, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TGaxis, TPad, TLatex, TBox
from ROOT import gStyle, gROOT, gSystem
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan
from histo_manager import slice_histogram, rebin_histogram
from painter import make_common_style
gStyle.SetOptStat(0);
gStyle.SetOptTitle(0);  
    
#__________________________________________________________________________________________________________
def DrawGammaPadSettings(   pad1,
                        leftMargin,
                        rightMargin,
                        topMargin,
                        bottomMargin):
    pad1.SetFillColor(0);
    pad1.GetFrame().SetFillColor(0);
    pad1.SetBorderMode(0);
    pad1.SetLeftMargin(leftMargin);
    pad1.SetBottomMargin(bottomMargin);
    pad1.SetRightMargin(rightMargin);
    pad1.SetTopMargin(topMargin);
    pad1.SetTickx();
    pad1.SetTicky();

def DrawGammaCanvasSettings(    c1,
                                leftMargin,
                                rightMargin,
                                topMargin,
                                bottomMargin):
        c1.SetTickx();
        c1.SetTicky();
        c1.SetGridx(0);
        c1.SetGridy(0);
        c1.SetLogy(0);
        c1.SetLeftMargin(leftMargin);
        c1.SetRightMargin(rightMargin);
        c1.SetTopMargin(topMargin);
        c1.SetBottomMargin(bottomMargin);
        c1.SetFillColor(0);

def DrawGammaHisto(     histo1,
                        Title,
                        XTitle,
                        YTitle,
                        xMin,
                        xMax,
                        bck,
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
            
        
        if (xMin > 0.2):  
            histo1.GetYaxis().SetRangeUser(yMin, 1.5*yMax);
        else: 
            histo1.GetYaxis().SetRangeUser(yMin, 1.2*yMax);
        

        if XTitle.Length() > 0:
            histo1.SetXTitle(XTitle.Data());
        
        if YTitle.Length() > 0:
            histo1.SetYTitle(YTitle.Data());
        
        histo1.GetYaxis().SetLabelSize(0.02);
        histo1.GetYaxis().SetTitleSize(0.025);
        histo1.GetYaxis().SetDecimals();
        histo1.GetYaxis().SetTitleOffset(0.5);
        histo1.GetXaxis().SetTitleSize(0.025);
        histo1.GetXaxis().SetLabelSize(0.02);
        histo1.SetMarkerStyle(20)
        histo1.SetMarkerColor(1)
        histo1.SetLineColor(1)
        histo1.SetLineWidth(1)
        histo1.SetMarkerSize(markerSize)
        histo1.SetTitleOffset(1.2, "XY")
        histo1.SetTitleSize(0.05, "XY")
        histo1.GetYaxis().SetLabelSize(0.05)
        histo1.GetXaxis().SetLabelSize(0.05)
        histo1.GetXaxis().SetNdivisions(507, True)
        if bck == 1 :
            histo1.SetLineStyle(1);
            histo1.SetLineColor(4);
            histo1.SetMarkerColor(4);
            histo1.SetMarkerStyle(24);
            histo1.SetLineWidth(1);
            histo1.DrawCopy("hist,same");
        else:
            if bck == 2 :
                histo1.DrawCopy("same");
            elif bck == 3:
              histo1.SetLineStyle(1);
              histo1.SetLineColor(kRed+2);
              histo1.SetMarkerColor(kRed+2);
              histo1.SetMarkerStyle(24);
              histo1.SetLineWidth(1);
              histo1.DrawCopy("hist,same");
            elif bck == 4:
              histo1.SetLineStyle(1);
              histo1.SetLineColor(kGreen+3);
              histo1.SetMarkerColor(kGreen+3);
              histo1.SetMarkerStyle(24);
              histo1.SetLineWidth(1);
              histo1.DrawCopy("hist,same");
            elif bck == 5:
              histo1.SetLineStyle(1);
              histo1.SetLineColor(kCyan+3);
              histo1.SetMarkerColor(kCyan+3);
              histo1.SetMarkerStyle(24);
              histo1.SetLineWidth(1);
              histo1.DrawCopy("hist,same");
            elif bck == 6:
              histo1.SetLineStyle(1);
              histo1.SetLineColor(kGreen+1);
              histo1.SetMarkerColor(kGreen+1);
              histo1.SetMarkerStyle(24);
              histo1.SetLineWidth(1);
              histo1.DrawCopy("hist,same");
            elif bck == -1:
              histo1.SetLineStyle(1);
              histo1.SetLineColor(11);
              histo1.SetMarkerColor(11);
              histo1.SetMarkerStyle(24);
              histo1.SetLineWidth(1);
              histo1.DrawCopy("hist,same");
            else:
                if Title.Length() > 0:
                    histo1.SetTitle("");
                            
            histo1.DrawCopy("e1,p")

            if Title.Length() > 0:
                alice = TLatex(0.1, 0.95, "{}".format(Title.Data()))
                alice.SetNDC()
                alice.SetTextColor(1)
                alice.SetTextSize(0.062)
                alice.Draw()

def SetStyleTLatex( text,
                    textSize,
                    lineWidth,
                    textColor = 1,
                    textFont = 42,
                    kNDC = True,
                    align = 11
                    ):
        if kNDC:
             text.SetNDC();
        text.SetTextFont(textFont);
        text.SetTextColor(textColor);
        text.SetTextSize(textSize);
        text.SetLineWidth(lineWidth);
        text.SetTextAlign(align);



def PlotLabelsInvMassInPtPlots          (   startTextX,
                                            startTextY,
                                            textHeight,
                                            differenceText,
                                            textAlice,
                                            dateDummy,
                                            fEnergy,
                                            fDecayChannel,
                                            fDetectionChannel,
                                            textEvents          = "",
                                            fNEvents           = 0,
                                            textAlign             = 11
    ):

        alice           = TLatex(startTextX, startTextY, textAlice);
        SetStyleTLatex( alice, textHeight*1.3, 1, 1, 42, True, textAlign);
        alice.Draw();
        latexDate       = TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy.Data());
        SetStyleTLatex( latexDate, textHeight, 1, 1, 42, True, textAlign);
        latexDate.Draw();
        energy          = TLatex(startTextX, (startTextY-2.25*differenceText), fEnergy);
        SetStyleTLatex( energy, textHeight*1, 1, 1, 42, True, textAlign);
        energy.Draw();
        process         = TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
        SetStyleTLatex( process, textHeight*1, 1, 1, 42, True, textAlign);
        process.Draw();
        detprocess      = TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionChannel);
        SetStyleTLatex( detprocess, textHeight*1, 1, 1, 42, True, textAlign);
        detprocess.Draw();
        events          = TLatex();
        if textEvents.CompareTo("") != 0 and fNEvents > 0:
            events              = TLatex(startTextX, (startTextY-5.25*differenceText), "{}: {:.1e} events".format(textEvents, fNEvents));
            SetStyleTLatex( events, textHeight*1, 1, 1, 42, True, textAlign);
            events.Draw();
        
    

    
    
    
    
#__________________________________________ Plotting all Invariant Mass bins _______________________________________________
def PlotInvMassInPtBins(    fHistoMappingGGInvMassPtBinPlot,
                            fHistoMappingBackNormInvMassPtBinPlot,
                            namePlot,
                            nameCanvas,
                            namePad,
                            fPlottingRangeMeson,
                            dateDummy,
                            fRowPlot,
                            fColumnPlot,
                            fStartBinPtRange,
                            fNumberPtBins,
                            fRangeBinsPt,
                            fDecayChannel,
                            fMonteCarloInfo,
                            decayChannel,
                            fDetectionChannel,
                            fEnergy,
                            optionBackground        = "mixed evt.",
                            isVsPtConv               = False,
                            BckNmb                    = 0,
                            fPlottingType           = "performance"
                            ):

        TGaxis.SetMaxDigits(3);

        canvasDataSpectra = TCanvas(nameCanvas, "", 1400, 900)
        DrawGammaCanvasSettings( canvasDataSpectra, 0, 0, 0, 0);
        canvasDataSpectra.cd()
    
        padDataSpectra = TPad(namePad, "", -0.0, 0.0, 1.0, 1.0, 0)
        DrawGammaPadSettings( padDataSpectra, 0, 0, 0, 0);
        padDataSpectra.Divide(fColumnPlot, fRowPlot, 0.0, 0.0)
        padDataSpectra.Draw()


        place                     = 0;
        legendPlace            = [fColumnPlot, fColumnPlot];
        if (fColumnPlot > 7):
            legendPlace[0]              = fColumnPlot-1;
        print(fColumnPlot);


        for iPt in range(fStartBinPtRange, fNumberPtBins):
            startPt = fRangeBinsPt[iPt];
            endPt = fRangeBinsPt[iPt + 1];

            place += 1;

            if place > legendPlace[0] - 1 and place < legendPlace[1] + 1:
                iPt -= 1;
            else:
                padDataSpectra.cd(place);
                padDataSpectra.cd(place).SetTopMargin(0.12);
                padDataSpectra.cd(place).SetBottomMargin(0.15);
                padDataSpectra.cd(place).SetRightMargin(0.02);
                remaining = int((place - 1) % fColumnPlot);

                if remaining > 0:
                    padDataSpectra.cd(place).SetLeftMargin(0.15);
                else:
                    padDataSpectra.cd(place).SetLeftMargin(0.25);

                titlePt = "{:.2f} GeV/c < pT < {:.2f} GeV/c".format(startPt, endPt);

                if isVsPtConv:
                    titlePt = "{:.2f} GeV/c < pT,conv < {:.2f} GeV/c".format(startPt, endPt);

                if "SubPiZero" in namePlot:
                    DrawGammaHisto(fHistoMappingGGInvMassPtBinPlot[iPt], titlePt,
                                "M_{} (GeV/c^2)".format(decayChannel), "dN_{}/dM_{}".format(decayChannel, decayChannel),
                                fPlottingRangeMeson[0], fPlottingRangeMeson[1], 0)

                    DrawGammaHisto(fHistoMappingBackNormInvMassPtBinPlot[iPt], titlePt,
                                "M_{} (GeV/c^2)".format(decayChannel), "dN_{}/dM_{}".format(decayChannel, decayChannel),
                                fPlottingRangeMeson[0], fPlottingRangeMeson[1], 1)
                else:
                    DrawGammaHisto(fHistoMappingGGInvMassPtBinPlot[iPt], titlePt,
                                "M_{} (GeV/c^2)".format(decayChannel), "dN_{}/dM_{}".format(decayChannel, decayChannel),
                                fPlottingRangeMeson[0], fPlottingRangeMeson[1], 0)

                    DrawGammaHisto(fHistoMappingBackNormInvMassPtBinPlot[iPt], titlePt,
                                "M_{} (GeV/c^2)".format(decayChannel), "dN_{}/dM_{}".format(decayChannel, decayChannel),
                                fPlottingRangeMeson[0], fPlottingRangeMeson[1], 1)

                fBGFitRangeLow = fBGFitRange[0]
                fBGFitRangeHigh = fBGFitRange[1]

                if "Left" in namePlot:
                    fBGFitRangeLow = fBGFitRangeLeft[0]
                    fBGFitRangeHigh = fBGFitRangeLeft[1]
                elif "SubPiZero" in namePlot:
                    fBGFitRangeLow = fBGFitRange_SubPiZero[0]
                    fBGFitRangeHigh = fBGFitRange_SubPiZero[1]
                elif "FixedPzPiZero" in namePlot:
                    fBGFitRangeLow = fBGFitRange_FixedPzPiZero[0]
                    fBGFitRangeHigh = fBGFitRange_FixedPzPiZero[1]

                box = TBox(fBGFitRangeLow, fHistoMappingGGInvMassPtBinPlot[iPt].GetMaximum() * 0.93,
                                fBGFitRangeHigh, fHistoMappingGGInvMassPtBinPlot[iPt].GetMaximum() * 0.91)
                box.SetFillStyle(1001)
                box.SetFillColor(ROOT.kAzure + 9)
                box.Draw("same")

        canvasDataSpectra.cd();
        nPixels        = 13;
        textHeight     = 0.08;
        startTextX     = 0.10;
        columnsLegend     = 1;
        widthLegend    = 1./fColumnPlot;
        heightLegend   = 1./fRowPlot;
        marginWidthLeg = 0.15;
        exampleBin        = fColumnPlot+fStartBinPtRange-1;
        if fColumnPlot > 7:
            startTextX          = 0.05;
            nPixels             = 12;
            widthLegend         = 2./fColumnPlot;
            marginWidthLeg      = 0.25;

        # plotting Legend
        padLegend                = TPad("dummyPad","",1-widthLegend,1-heightLegend,1.,1.,0);   # gives the size of the histo areas
        DrawGammaPadSettings( padLegend, 0, 0, 0, 0);
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

        textEvents="";
        if fMonteCarloInfo:
            textEvents          = "MC";
        else:
            textEvents          = "Data";
    

        if  padLegend.XtoPixel(padLegend.GetX2()) < padLegend.YtoPixel(padLegend.GetY1()):
            textHeight          = nPixels/padLegend.XtoPixel(padLegend.GetX2()) ;
        else:
            textHeight          = nPixels/padLegend.YtoPixel(padLegend.GetY1());
        
        startTextY     = 0.9;
        differenceText = textHeight*1.05;
        # plot labels
        PlotLabelsInvMassInPtPlots ( startTextX, startTextY, textHeight, differenceText, textAlice, dateDummy, fEnergy, fDecayChannel, fDetectionChannel, textEvents, fNEvents);


        legendData = ROOT.TLegend(startTextX, startTextY - 5.75 * differenceText, 0.85, startTextY - (5.75 + 2 / columnsLegend) * differenceText, "", 43)
        legendData.SetTextSizePixels(nPixels)
        legendData.SetNColumns(columnsLegend)
        legendData.SetMargin(marginWidthLeg)

        markersize = fHistoMappingGGInvMassPtBinPlot[exampleBin].GetMarkerSize()
        fHistoMappingGGInvMassPtBinPlot[exampleBin].SetMarkerSize(3 * markersize)
        legendData.AddEntry(fHistoMappingGGInvMassPtBinPlot[exampleBin], "same evt. M_{} (BG+Signal)".format(decayChannel), "ep")

        linesize = fHistoMappingBackNormInvMassPtBinPlot[exampleBin].GetLineWidth()
        fHistoMappingBackNormInvMassPtBinPlot[exampleBin].SetLineWidth(5 * linesize)

        if BckNmb == 0:
            if namePlot.Contains("FixedPzPiZero") == True:
                legendData.AddEntry(fHistoMappingBackNormInvMassPtBinPlot[exampleBin], "{} #it{{M}}_{} (p_{{z}} of #pi^{{0}} fixed)".format(optionBackground, decayChannel), "l")
            else:
                legendData.AddEntry(fHistoMappingBackNormInvMassPtBinPlot[exampleBin], "{} #it{{M}}_{}".format(optionBackground, decayChannel), "l")
        else:
            if namePlot.Contains("FixedPzPiZero") == True:
                legendData.AddEntry(fHistoMappingBackNormInvMassPtBinPlot[exampleBin], "{} #it{{M}}_{} group {} (p_{{z}} of #pi^{{0}} fixed)".format(optionBackground, decayChannel, BckNmb), "l")
            else:
                legendData.AddEntry(fHistoMappingBackNormInvMassPtBinPlot[exampleBin], "{} #it{{M}}_{} group {}".format(optionBackground, decayChannel, BckNmb), "l")

        legendData.Draw();

        canvasDataSpectra.SaveAs(namePlot.Data());
        del padLegend;
        del padDataSpectra;
        del canvasDataSpectra;