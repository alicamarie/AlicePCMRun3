    //__________________________________________ Plotting all Invariant Mass bins _______________________________________________
    void PlotInvMassInPtBins(   TH1D** fHistoMappingGGInvMassPtBinPlot,
                                TH1D** fHistoMappingBackNormInvMassPtBinPlot,
                                TString namePlot,
                                TString nameCanvas,
                                TString namePad,
                                Double_t* fPlottingRangeMeson,
                                TString dateDummy,
                                Int_t fRowPlot,
                                Int_t fColumnPlot,
                                Int_t fStartBinPtRange,
                                Int_t fNumberPtBins,
                                Double_t* fRangeBinsPt,
                                TString fDecayChannel,
                                Bool_t fMonteCarloInfo,
                                TString decayChannel,
                                TString fDetectionChannel,
                                TString fEnergy,
                                TString optionBackground        = "mixed evt.",
                                Bool_t isVsPtConv               = kFALSE,
                                Int_t BckNmb                    = 0,
                                TString fPlottingType           = "performance"
                            ){
        TGaxis::SetMaxDigits(3);

        TCanvas *canvasDataSpectra          = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
        DrawGammaCanvasSettings( canvasDataSpectra, 0, 0, 0, 0);
        canvasDataSpectra->cd();

        TPad * padDataSpectra               = new TPad(namePad.Data(),"",-0.0,0.0,1.,1.,0);   // gives the size of the histo areas
        DrawGammaPadSettings( padDataSpectra, 0, 0, 0, 0);
        padDataSpectra->Divide(fColumnPlot,fRowPlot,0.0,0.0);
        padDataSpectra->Draw();

        Int_t place                     = 0;
        Int_t legendPlace[2]            = {fColumnPlot, fColumnPlot};
        if (fColumnPlot > 7)
            legendPlace[0]              = fColumnPlot-1;
        cout << fColumnPlot << endl;
        for(Int_t iPt=fStartBinPtRange;iPt<fNumberPtBins;iPt++){
            //         cout<<"Pt: "<<iPt<<" of "<<fNumberPtBins<<endl;
            Double_t startPt            = fRangeBinsPt[iPt];
            Double_t endPt              = fRangeBinsPt[iPt+1];

            place                       = place + 1; //give the right place in the page
            if ( place> legendPlace[0]-1 && place < legendPlace[1]+1 ){
                iPt--;
            } else {
                padDataSpectra->cd(place);
                padDataSpectra->cd(place)->SetTopMargin(0.12);
                padDataSpectra->cd(place)->SetBottomMargin(0.15);
                padDataSpectra->cd(place)->SetRightMargin(0.02);
                int remaining           = (int)((place-1)%fColumnPlot);
                if (remaining > 0) padDataSpectra->cd(place)->SetLeftMargin(0.15);
                else padDataSpectra->cd(place)->SetLeftMargin(0.25);
    //             cout << "here" << endl;
                TString titlePt = Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt);
                if (isVsPtConv)
                    titlePt     = Form("%3.2f GeV/#it{c} < #it{p}_{T,#gamma_{conv}} < %3.2f GeV/#it{c}",startPt,endPt);
                if(namePlot.Contains("SubPiZero") == kTRUE){
                  DrawGammaHisto( fHistoMappingGGInvMassPtBinPlot[iPt],
                                  titlePt,
                                  Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                                  fPlottingRangeMeson[0],fPlottingRangeMeson[1],0);
                  //             cout << "here" << endl;
                  DrawGammaHisto( fHistoMappingBackNormInvMassPtBinPlot[iPt],
                                  titlePt,
                                  Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                                  fPlottingRangeMeson[0],fPlottingRangeMeson[1],1);
                  //             cout << "here" << endl;
                } else{
                  DrawGammaHisto( fHistoMappingGGInvMassPtBinPlot[iPt],
                                  titlePt,
                                  Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                                  fPlottingRangeMeson[0],fPlottingRangeMeson[1],0);
                  //             cout << "here" << endl;
                  DrawGammaHisto( fHistoMappingBackNormInvMassPtBinPlot[iPt],
                                  titlePt,
                                  Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                                  fPlottingRangeMeson[0],fPlottingRangeMeson[1],1);
                  //             cout << "here" << endl;
                }
                Double_t fBGFitRangeLow     = fBGFitRange[0];
                Double_t fBGFitRangeHigh    = fBGFitRange[1];
                if (namePlot.Contains("Left")){
                  fBGFitRangeLow          = fBGFitRangeLeft[0];
                  fBGFitRangeHigh         = fBGFitRangeLeft[1];

                } else if (namePlot.Contains("SubPiZero")){
                  fBGFitRangeLow          = fBGFitRange_SubPiZero[0];
                  fBGFitRangeHigh         = fBGFitRange_SubPiZero[1];
                } else if (namePlot.Contains("FixedPzPiZero")){
                  fBGFitRangeLow          = fBGFitRange_FixedPzPiZero[0];
                  fBGFitRangeHigh         = fBGFitRange_FixedPzPiZero[1];
                }
                TBox *box               = new TBox(fBGFitRangeLow,fHistoMappingGGInvMassPtBinPlot[iPt]->GetMaximum()*0.93,fBGFitRangeHigh,fHistoMappingGGInvMassPtBinPlot[iPt]->GetMaximum()*0.91);
                box->SetFillStyle(1001);
                box->SetFillColor(kAzure+9);
                box->Draw("same");
            }
        }

        canvasDataSpectra->cd();
        Double_t nPixels        = 13;
        Double_t textHeight     = 0.08;
        Double_t startTextX     = 0.10;
        Int_t columnsLegend     = 1;
        Double_t widthLegend    = 1./fColumnPlot;
        Double_t heightLegend   = 1./fRowPlot;
        Double_t marginWidthLeg = 0.15;
        Int_t exampleBin        = fColumnPlot+fStartBinPtRange-1;
        if (fColumnPlot > 7){
            startTextX          = 0.05;
            nPixels             = 12;
            widthLegend         = 2./fColumnPlot;
            marginWidthLeg      = 0.25;
        }

        // plotting Legend
        TPad * padLegend                = new TPad("dummyPad","",1-widthLegend,1-heightLegend,1.,1.,0);   // gives the size of the histo areas
        DrawGammaPadSettings( padLegend, 0, 0, 0, 0);
        padLegend->Draw();
        padLegend->cd();

        TString textAlice;
        if(fPlottingType.CompareTo("wip")==0){
            textAlice       = "ALICE work in progress";
        } else if (fPlottingType.CompareTo("thesis")==0){
            textAlice       = "ALICE this thesis";
        } else if (fPlottingType.CompareTo("performance")==0){
            textAlice       = "ALICE performance";
        } else{
            textAlice       = "ALICE";
        }
        TString textEvents;
        if(fMonteCarloInfo){
            textEvents          = "MC";
        } else {
            textEvents          = "Data";
        }

        if (padLegend->XtoPixel(padLegend->GetX2()) < padLegend->YtoPixel(padLegend->GetY1())){
            textHeight          = (Double_t)nPixels/padLegend->XtoPixel(padLegend->GetX2()) ;
        } else {
            textHeight          = (Double_t)nPixels/padLegend->YtoPixel(padLegend->GetY1());
        }
        Double_t startTextY     = 0.9;
        Double_t differenceText = textHeight*1.05;
        // plot labels
        PlotLabelsInvMassInPtPlots ( startTextX, startTextY, textHeight, differenceText, textAlice, dateDummy, fEnergy, fDecayChannel, fDetectionChannel, textEvents, fNEvents);

        TLegend* legendData     = GetAndSetLegend2(  startTextX, startTextY-5.75*differenceText, 0.85,  startTextY-(5.75+2/columnsLegend)*differenceText, nPixels, columnsLegend, "", 43, marginWidthLeg);
        Size_t markersize       = fHistoMappingGGInvMassPtBinPlot[exampleBin]->GetMarkerSize();
        fHistoMappingGGInvMassPtBinPlot[exampleBin]->SetMarkerSize(3*markersize);
        legendData->AddEntry(fHistoMappingGGInvMassPtBinPlot[exampleBin],Form("same evt. #it{M}_{%s} (BG+Signal)",decayChannel.Data()),"ep");
        Size_t linesize         = fHistoMappingBackNormInvMassPtBinPlot[exampleBin]->GetLineWidth();
        fHistoMappingBackNormInvMassPtBinPlot[exampleBin]->SetLineWidth(5*linesize);
        if(BckNmb==0){
            if(namePlot.Contains("FixedPzPiZero") == kTRUE){
                legendData->AddEntry(fHistoMappingBackNormInvMassPtBinPlot[exampleBin],Form("%s #it{M}_{%s} (p_{z} of #pi^{0} fixed)",optionBackground.Data(), decayChannel.Data()),"l");

            } else{
                legendData->AddEntry(fHistoMappingBackNormInvMassPtBinPlot[exampleBin],Form("%s #it{M}_{%s}",optionBackground.Data(), decayChannel.Data()),"l");
            }
        } else{
            if(namePlot.Contains("FixedPzPiZero") == kTRUE){
                legendData->AddEntry(fHistoMappingBackNormInvMassPtBinPlot[exampleBin],Form("%s #it{M}_{%s} group %d (p_{z} of #pi^{0} fixed)",optionBackground.Data(), decayChannel.Data(),BckNmb),"l");
            } else{
                legendData->AddEntry(fHistoMappingBackNormInvMassPtBinPlot[exampleBin],Form("%s #it{M}_{%s} group %d",optionBackground.Data(), decayChannel.Data(),BckNmb),"l");
            }
        }
        legendData->Draw();

        canvasDataSpectra->SaveAs(namePlot.Data());
        delete padLegend;
        delete padDataSpectra;
        delete canvasDataSpectra;
    }

    TLegend *GetAndSetLegend2(  Double_t positionX,
                                Double_t positionY,
                                Double_t positionXRight,
                                Double_t positionYUp,
                                Size_t textSize,
                                Int_t columns               = 1,
                                TString header              = "",
                                Font_t textFont             = 43,
                                Double_t margin             = 0
    ){

        TLegend *legend = new TLegend(positionX,positionY,positionXRight,positionYUp);
        legend->SetNColumns(columns);
        legend->SetLineColor(0);
        legend->SetLineWidth(0);
        legend->SetFillColor(0);
        legend->SetFillStyle(0);
        legend->SetLineStyle(0);
        legend->SetBorderSize(0);
        legend->SetTextFont(textFont);
        legend->SetTextSize(textSize);
        if (margin != 0) legend->SetMargin(margin);
        if (header.CompareTo("")!= 0) legend->SetHeader(header);
        return legend;
    }

    //__________________________________________ Plotting all Invariant Mass bins _______________________________________________
    void PlotInvMassInPtBins(   TH1D** fHistoMappingBackWithRemainNormInvMassPtBinPlot,
                                TH1D** fHistoMappingBackNormInvMassPtBinPlot,
                                TH1D** fHistoMappingTrueAllBackNormInvMassPtBinPlot,
                                TH1D** fHistoMappingTrueGGBackNormInvMassPtBinPlot,
                                TH1D** fHistoMappingTrueContBackNormInvMassPtBinPlot,
                                TH1D** fHistoMappingTrueMesonContainedInvMassPtBins,
                                TH1D** fHistoMappingTrueAsymEClusInvMassPtBins,
                                TString namePlot,
                                TString nameCanvas,
                                TString namePad,
                                Double_t* fPlottingRangeMeson,
                                TString dateDummy,
                                TString fMesonType,
                                Int_t fRowPlot,
                                Int_t fColumnPlot,
                                Int_t fStartBinPtRange,
                                Int_t fNumberPtBins,
                                Double_t* fRangeBinsPt,
                                TString fDecayChannel,
                                Bool_t fMonteCarloInfo,
                                TString decayChannel,
                                TString fDetectionChannel,
                                TString fEnergy,
                                TString optionBackground        = "mixed evt.",
                                Bool_t isVsPtConv               = kFALSE,
                                TString fPlottingType ="performance",
                                Bool_t doPtDependentIntDeltaRange = kFALSE,
                                Double_t*** dArrPtDependentIntDeltaRange = NULL
                            ){
        TGaxis::SetMaxDigits(3);

        TCanvas *canvasDataSpectra          = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
        DrawGammaCanvasSettings( canvasDataSpectra, 0, 0, 0, 0);
        canvasDataSpectra->cd();

        TPad * padDataSpectra               = new TPad(namePad.Data(),"",-0.0,0.0,1.,1.,0);   // gives the size of the histo areas
        DrawGammaPadSettings( padDataSpectra, 0, 0, 0, 0);
        padDataSpectra->Divide(fColumnPlot,fRowPlot,0.0,0.0);
        padDataSpectra->Draw();

        Int_t place                     = 0;
        Int_t legendPlace[2]            = {fColumnPlot, fColumnPlot};
        if (fColumnPlot > 7)
            legendPlace[0]              = fColumnPlot-1;
        for(Int_t iPt=fStartBinPtRange;iPt<fNumberPtBins;iPt++){
            //         cout<<"Pt: "<<iPt<<" of "<<fNumberPtBins<<endl;
            Double_t startPt            = fRangeBinsPt[iPt];
            Double_t endPt              = fRangeBinsPt[iPt+1];

            place                       = place + 1; //give the right place in the page
            if (place == legendPlace[0]){
                iPt--;
                padDataSpectra->cd(place);

    //             cout << "entered ALICE plotting" << endl;
                TString textAlice;
                if(fPlottingType.CompareTo("wip")==0){
                    textAlice       = "ALICE work in progress";
                } else if (fPlottingType.CompareTo("thesis")==0){
                    textAlice       = "ALICE this thesis";
                } else if (fPlottingType.CompareTo("performance")==0){
                    textAlice       = "ALICE performance";
                } else{
                    textAlice       = "ALICE";
                }
                TString textEvents;
                if(fMonteCarloInfo){
                    textEvents          = "MC";
                } else {
                    textEvents          = "Data";
                }
                Double_t nPixels        = 13;
                Double_t textHeight     = 0.08;
                if (padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) < padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1())){
                    textHeight          = (Double_t)nPixels/padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) ;
                } else {
                    textHeight          = (Double_t)nPixels/padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1());
                }

                Double_t startTextX     = 0.10;
                Double_t startTextY     = 0.9;
                Double_t differenceText = textHeight*1.25;

                PlotLabelsInvMassInPtPlots ( startTextX, startTextY, textHeight, differenceText, textAlice, dateDummy, fEnergy, fDecayChannel, fDetectionChannel, textEvents, fNEvents);

                TLegend* legendData     = new TLegend(startTextX,startTextY-5.75*differenceText,1,startTextY-(5.75+5.)*differenceText);
                legendData->SetTextSize(textHeight);
                legendData->SetTextFont(62);
                legendData->SetFillColor(0);
                legendData->SetFillStyle(0);
                legendData->SetLineWidth(0);
                legendData->SetLineColor(0);
                legendData->SetMargin(0.15);
                Size_t markersize       = fHistoMappingBackWithRemainNormInvMassPtBinPlot[iPt]->GetMarkerSize();
                fHistoMappingBackWithRemainNormInvMassPtBinPlot[iPt]->SetMarkerSize(3*markersize);
                legendData->AddEntry(fHistoMappingBackWithRemainNormInvMassPtBinPlot[iPt],Form("%s + rem. bck.", optionBackground.Data()),"ep");

                fHistoMappingBackNormInvMassPtBinPlot[iPt]->SetMarkerSize(3*markersize);
                legendData->AddEntry(fHistoMappingBackNormInvMassPtBinPlot[iPt],optionBackground,"ep");

                Size_t linesize         = fHistoMappingTrueAllBackNormInvMassPtBinPlot[iPt]->GetLineWidth();
                fHistoMappingTrueAllBackNormInvMassPtBinPlot[iPt]->SetLineWidth(5*linesize);
                legendData->AddEntry(fHistoMappingTrueAllBackNormInvMassPtBinPlot[iPt],Form("true %s ALL", optionBackground.Data()),"l");
                legendData->Draw();

                fHistoMappingTrueGGBackNormInvMassPtBinPlot[iPt]->SetLineWidth(5*linesize);
                legendData->AddEntry(fHistoMappingTrueGGBackNormInvMassPtBinPlot[iPt],Form("true %s GG", optionBackground.Data()),"l");
                legendData->Draw();

                fHistoMappingTrueContBackNormInvMassPtBinPlot[iPt]->SetLineWidth(5*linesize);
                legendData->AddEntry(fHistoMappingTrueContBackNormInvMassPtBinPlot[iPt],Form("true %s cont", optionBackground.Data()),"l");
                legendData->Draw();

                if(fHistoMappingTrueMesonContainedInvMassPtBins[iPt]){
                  fHistoMappingTrueMesonContainedInvMassPtBins[iPt]->SetLineWidth(5*linesize);
                  legendData->AddEntry(fHistoMappingTrueMesonContainedInvMassPtBins[iPt],"true full meson contained","l");
                  legendData->Draw();
                }
                if(fHistoMappingTrueAsymEClusInvMassPtBins[iPt]){
                  fHistoMappingTrueAsymEClusInvMassPtBins[iPt]->SetLineWidth(5*linesize);
                  legendData->AddEntry(fHistoMappingTrueAsymEClusInvMassPtBins[iPt],"true asym E_clus cont","l");
                  legendData->Draw();
                }
            } else if (place < legendPlace[1]+1 && place > legendPlace[0] ){
                iPt--;
            } else {

                padDataSpectra->cd(place);
                padDataSpectra->cd(place)->SetTopMargin(0.12);
                padDataSpectra->cd(place)->SetBottomMargin(0.15);
                padDataSpectra->cd(place)->SetRightMargin(0.02);
    //             cout << "place: " << place << endl;
                int remaining           = (place-1)%fColumnPlot;
                if (remaining > 0) padDataSpectra->cd(place)->SetLeftMargin(0.15);
                else padDataSpectra->cd(place)->SetLeftMargin(0.25);
    //             cout << "here" << endl;
                TString titlePt = Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt);
                if (isVsPtConv)
                    titlePt     = Form("%3.2f GeV/#it{c} < #it{p}_{T,#gamma_{conv}} < %3.2f GeV/#it{c}",startPt,endPt);

                DrawGammaHisto( fHistoMappingBackWithRemainNormInvMassPtBinPlot[iPt],
                                titlePt,
                                Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                                fPlottingRangeMeson[0],fPlottingRangeMeson[1],0);
                if (fHistoMappingBackWithRemainNormInvMassPtBinPlot[iPt]!=0x00){
                    TString nameOfPlot = fHistoMappingBackWithRemainNormInvMassPtBinPlot[iPt]->GetName();
                    Double_t mass = fMesonMass[iPt];
                    if (nameOfPlot.Contains("Left"))
                        mass                        = fMesonMassLeft[iPt];
                    if (nameOfPlot.Contains("True"))
                        mass                        = fMesonTrueMass[iPt];
                    Double_t intRangeLow            = mass + fMesonIntDeltaRange[0];
                    Double_t intRangeWideLow        = mass + fMesonIntDeltaRangeWide[0];
                    Double_t intRangeNarrowLow      = mass + fMesonIntDeltaRangeNarrow[0];
                    Double_t intRangeHigh           = mass + fMesonIntDeltaRange[1];
                    Double_t intRangeWideHigh       = mass + fMesonIntDeltaRangeWide[1];
                    Double_t intRangeNarrowHigh     = mass + fMesonIntDeltaRangeNarrow[1];

                    if (doPtDependentIntDeltaRange){
                        intRangeLow            = mass + dArrPtDependentIntDeltaRange[0][0][iPt];
                        intRangeWideLow        = mass + dArrPtDependentIntDeltaRange[1][0][iPt];
                        intRangeNarrowLow      = mass + dArrPtDependentIntDeltaRange[2][0][iPt];
                        intRangeHigh           = mass + dArrPtDependentIntDeltaRange[0][1][iPt];
                        intRangeWideHigh       = mass + dArrPtDependentIntDeltaRange[1][1][iPt];
                        intRangeNarrowHigh     = mass + dArrPtDependentIntDeltaRange[2][1][iPt];
                    }

                    Double_t normalLow              = intRangeLow-(intRangeLow-fHistoMappingBackWithRemainNormInvMassPtBinPlot[iPt]->GetXaxis()->GetBinLowEdge(fHistoMappingBackWithRemainNormInvMassPtBinPlot[iPt]->GetXaxis()->FindBin(intRangeLow)));
                    Double_t normalUp               = intRangeHigh+(fHistoMappingBackWithRemainNormInvMassPtBinPlot[iPt]->GetXaxis()->GetBinUpEdge(fHistoMappingBackWithRemainNormInvMassPtBinPlot[iPt]->GetXaxis()->FindBin(intRangeHigh))-intRangeHigh);
                    Double_t wideLow                = intRangeWideLow-(intRangeWideLow-fHistoMappingBackWithRemainNormInvMassPtBinPlot[iPt]->GetXaxis()->GetBinLowEdge(fHistoMappingBackWithRemainNormInvMassPtBinPlot[iPt]->GetXaxis()->FindBin(intRangeWideLow)));
                    Double_t wideUp                 = intRangeWideHigh+(fHistoMappingBackWithRemainNormInvMassPtBinPlot[iPt]->GetXaxis()->GetBinUpEdge(fHistoMappingBackWithRemainNormInvMassPtBinPlot[iPt]->GetXaxis()->FindBin(intRangeWideHigh))-intRangeWideHigh);
                    Double_t narrowLow              = intRangeNarrowLow-(intRangeNarrowLow-fHistoMappingBackWithRemainNormInvMassPtBinPlot[iPt]->GetXaxis()->GetBinLowEdge(fHistoMappingBackWithRemainNormInvMassPtBinPlot[iPt]->GetXaxis()->FindBin(intRangeNarrowLow)));
                    Double_t narrowUp               = intRangeNarrowHigh+(fHistoMappingBackWithRemainNormInvMassPtBinPlot[iPt]->GetXaxis()->GetBinUpEdge(fHistoMappingBackWithRemainNormInvMassPtBinPlot[iPt]->GetXaxis()->FindBin(intRangeNarrowHigh))-intRangeNarrowHigh);

                    DrawGammaLines(mass, mass, fHistoMappingBackWithRemainNormInvMassPtBinPlot[iPt]->GetMinimum(), fHistoMappingBackWithRemainNormInvMassPtBinPlot[iPt]->GetMaximum()*0.4, 1, kRed+2, 1);
                    DrawGammaLines(normalLow, normalLow, fHistoMappingBackWithRemainNormInvMassPtBinPlot[iPt]->GetMinimum(), fHistoMappingBackWithRemainNormInvMassPtBinPlot[iPt]->GetMaximum(), 1, kGray+1, 1);
                    DrawGammaLines(normalUp, normalUp, fHistoMappingBackWithRemainNormInvMassPtBinPlot[iPt]->GetMinimum(), fHistoMappingBackWithRemainNormInvMassPtBinPlot[iPt]->GetMaximum(), 1, kGray+1, 1);
                    DrawGammaLines(wideLow, wideLow, fHistoMappingBackWithRemainNormInvMassPtBinPlot[iPt]->GetMinimum(), fHistoMappingBackWithRemainNormInvMassPtBinPlot[iPt]->GetMaximum(), 1, kGray+1, 2);
                    DrawGammaLines(wideUp, wideUp, fHistoMappingBackWithRemainNormInvMassPtBinPlot[iPt]->GetMinimum(), fHistoMappingBackWithRemainNormInvMassPtBinPlot[iPt]->GetMaximum(), 1, kGray+1, 2);
                    DrawGammaLines(narrowLow, narrowLow, fHistoMappingBackWithRemainNormInvMassPtBinPlot[iPt]->GetMinimum(), fHistoMappingBackWithRemainNormInvMassPtBinPlot[iPt]->GetMaximum(), 1, kGray+1, 3);
                    DrawGammaLines(narrowUp, narrowUp, fHistoMappingBackWithRemainNormInvMassPtBinPlot[iPt]->GetMinimum(), fHistoMappingBackWithRemainNormInvMassPtBinPlot[iPt]->GetMaximum(), 1, kGray+1, 3);

                }
                TString xlabel;
                TString ylabel;

                if(namePlot.Contains("SubPiZero") == kTRUE){
                  xlabel = Form("#it{M}_{%s} - #it{M}_{#pi^{0}} (GeV/#it{c}^{2})",decayChannel.Data());
                  ylabel = Form("dN_{%s}/d#it{M}_{#pi^{+} #pi^{-}}", decayChannel.Data());
                } else{
                  xlabel = Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data());
                  ylabel = Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data());
                }
                DrawGammaHisto( fHistoMappingBackNormInvMassPtBinPlot[iPt],
                                titlePt,
                                xlabel.Data(),ylabel.Data(),
                                fPlottingRangeMeson[0],fPlottingRangeMeson[1],-1);
    //             cout << "here" << endl;
                DrawGammaHisto( fHistoMappingTrueAllBackNormInvMassPtBinPlot[iPt],
                                titlePt,
                                xlabel.Data(),ylabel.Data(),
                                fPlottingRangeMeson[0],fPlottingRangeMeson[1],1);
    //             cout << "here" << endl;
                DrawGammaHisto( fHistoMappingTrueGGBackNormInvMassPtBinPlot[iPt],
                                titlePt,
                                xlabel.Data(),ylabel.Data(),
                                fPlottingRangeMeson[0],fPlottingRangeMeson[1],3);
    //             cout << "here" << endl;
                DrawGammaHisto( fHistoMappingTrueContBackNormInvMassPtBinPlot[iPt],
                                titlePt,
                                xlabel.Data(),ylabel.Data(),
                                fPlottingRangeMeson[0],fPlottingRangeMeson[1],4);
                if(fHistoMappingTrueMesonContainedInvMassPtBins[iPt]){
                  DrawGammaHisto( fHistoMappingTrueMesonContainedInvMassPtBins[iPt],
                                  titlePt,
                                  xlabel,ylabel.Data(),
                                  fPlottingRangeMeson[0],fPlottingRangeMeson[1],5);
                }
                if(fHistoMappingTrueAsymEClusInvMassPtBins[iPt]){
                  DrawGammaHisto( fHistoMappingTrueAsymEClusInvMassPtBins[iPt],
                                  titlePt,
                                  xlabel.Data(),ylabel.Data(),
                                  fPlottingRangeMeson[0],fPlottingRangeMeson[1],6);
                }
    //             cout << "here" << endl;
            }
        }
        canvasDataSpectra->Print(namePlot.Data());
        delete padDataSpectra;
        delete canvasDataSpectra;
    }


void DrawGammaHisto(    TH1* histo1,
                        TString Title,
                        TString XTitle,
                        TString YTitle,
                        Float_t xMin,
                        Float_t xMax,
                        Int_t bck,
                        Size_t markerSize = 0.2
                ) {

    histo1->GetXaxis()->SetRangeUser(xMin, xMax);
    Double_t yMin = 0;
    Double_t yMax = 0;
    for (Int_t i = histo1->GetXaxis()->FindBin(xMin); i < histo1->GetXaxis()->FindBin(xMax); i++){
        if (histo1->GetBinContent(i) < yMin){
            yMin = histo1->GetBinContent(i);
        }
        if (histo1->GetBinContent(i) > yMax){
            yMax = histo1->GetBinContent(i);
        }
    }
    if (xMin > 0.2)  histo1->GetYaxis()->SetRangeUser(yMin, 1.5*yMax);
        else histo1->GetYaxis()->SetRangeUser(yMin, 1.2*yMax);
    

    if(XTitle.Length() > 0){
        histo1->SetXTitle(XTitle.Data());
    }
    if(YTitle.Length() > 0){
        histo1->SetYTitle(YTitle.Data());
    }
    histo1->GetYaxis()->SetLabelSize(0.02);
    histo1->GetYaxis()->SetTitleSize(0.025);
    histo1->GetYaxis()->SetDecimals();
    histo1->GetYaxis()->SetTitleOffset(0.5);
    histo1->GetXaxis()->SetTitleSize(0.025);
    histo1->GetXaxis()->SetLabelSize(0.02);
    histo1->SetMarkerStyle(20);
    histo1->SetMarkerColor(1);
    histo1->SetLineColor(1);
    histo1->SetLineWidth(1);
    histo1->SetMarkerSize(markerSize);
    histo1->SetTitleOffset(1.2,"xy");
    histo1->SetTitleSize(0.05,"xy");
    histo1->GetYaxis()->SetLabelSize(0.05);
    histo1->GetXaxis()->SetLabelSize(0.05);
    histo1->GetXaxis()->SetNdivisions(507,kTRUE);
    if( bck == 1 ){
        histo1->SetLineStyle(1);
        histo1->SetLineColor(4);
        histo1->SetMarkerColor(4);
        histo1->SetMarkerStyle(24);
        histo1->SetLineWidth(1);
        histo1->DrawCopy("hist,same");
    } else {
        if( bck == 2 ){
            histo1->DrawCopy("same");
        }else if (bck == 3){
            histo1->SetLineStyle(1);
            histo1->SetLineColor(kRed+2);
            histo1->SetMarkerColor(kRed+2);
            histo1->SetMarkerStyle(24);
            histo1->SetLineWidth(1);
            histo1->DrawCopy("hist,same");
        }else if (bck == 4){
            histo1->SetLineStyle(1);
            histo1->SetLineColor(kGreen+3);
            histo1->SetMarkerColor(kGreen+3);
            histo1->SetMarkerStyle(24);
            histo1->SetLineWidth(1);
            histo1->DrawCopy("hist,same");
        }else if (bck == 5){
            histo1->SetLineStyle(1);
            histo1->SetLineColor(kCyan+3);
            histo1->SetMarkerColor(kCyan+3);
            histo1->SetMarkerStyle(24);
            histo1->SetLineWidth(1);
            histo1->DrawCopy("hist,same");
        }else if (bck == 6){
            histo1->SetLineStyle(1);
            histo1->SetLineColor(kGreen+1);
            histo1->SetMarkerColor(kGreen+1);
            histo1->SetMarkerStyle(24);
            histo1->SetLineWidth(1);
            histo1->DrawCopy("hist,same");
        }else if (bck == -1){
            histo1->SetLineStyle(1);
            histo1->SetLineColor(11);
            histo1->SetMarkerColor(11);
            histo1->SetMarkerStyle(24);
            histo1->SetLineWidth(1);
            histo1->DrawCopy("hist,same");
        } else {
            if(Title.Length() > 0){
                histo1->SetTitle("");
            }
            histo1->DrawCopy("e1,p");
            if(Title.Length() > 0){
                TLatex *alice = new TLatex(0.1,0.95,Form("%s",Title.Data())); // Bo: this was
                alice->SetNDC();
                alice->SetTextColor(1);
                alice->SetTextSize(0.062);
                alice->Draw();
            }
        }
    }
}