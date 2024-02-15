# This code is based on the code provided by Daiki Sekihata in his repository
# https://github.com/dsekihat/photon_sw_run3 for photon analysis in ALICE Run 3


import os, sys, shutil
import numpy as np
import math
import ctypes
import ROOT
from ROOT import TFile, TDirectory, THashList, TF1, TH1F, TCanvas, TMath
from histo_manager import slice_histogram, rebin_histogram, get_bkg_subtracted, get_ratio
from ctypes import *
#import file_manager

class PairAnalyzerMC:
    def __init__(self):
        print("default constructor is called");
    def __init__(self, meson, filename, dirname):
        print("target meson = {0} , filename = {1} , dirname = {2}".format(meson, filename, dirname));
        self.meson = meson;
        self.rootfile = TFile.Open(filename, "READ");
        self.rootdir = self.rootfile.Get(dirname);
        self.list_ev = self.rootdir.Get("Event");
        self.list_pair = self.rootdir.Get("Pair");
        self.arr_pt = np.array([0,1,2,3,4,5], dtype=float);
        # self.f1total = TF1("GaussExpLinear", 
        #        "(x<[1]) * ([0]*(TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))))+[4]+[5]*x)+\
        #        (x>=[1])*([0]*TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))+[4]+[5]*x)", 0,1);

        # old remaining things, will be deleted later
        #self.f1bkg = TF1("f1bkg","pol1(0)",0,1);
        self.f1total = TF1("fGaussExp",
                 "(x<[1]) * ([0]*(TMath::Exp(-0.5*TMath::Power((x-[1])/[2],2)) + TMath::Exp((x-[1])/[3])*(1. - TMath::Exp(-0.5*TMath::Power((x-[1])/[2],2))))) + \
                 (x>=[1]) * ([0]*TMath::Exp(-0.5*TMath::Power((x-[1])/[2],2)))", 0,1);

        self.f1total.SetNpx(1000);
        self.fit_min = 0.04;
        self.fit_max = 0.24;
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
        self.arr_pt = arr_pt;

    def set_subsystem(self, ssname):
        self.ssname = ssname;
        self.list_ev_ss   = self.list_ev.FindObject(ssname);
        self.list_pair_ss = self.list_pair.FindObject(ssname);

    def set_cutname(self, cutname):
        self.cutname = cutname;
        if self.list_ev_ss is None or self.list_pair_ss is None:
            print("Please define subsystem name first!");
            return None;
        #self.list_ev_ss_cut   = self.list_ev_ss.FindObject(cutname);
        self.list_pair_ss_cut = self.list_pair_ss.FindObject(cutname);

    def set_fit_range(self, fit_min, fit_max):
        self.fit_min = fit_min;
        self.fit_max = fit_max;
    
    def set_integral_range(self, integral_min, integral_max):
        self.integral_min = integral_min;
        self.integral_max = integral_max;

    def set_fit_function(self, func):
        #if func == "gausexplinear":
            # old remaining things, will be deleted later
        self.f1total = TF1("fGaussExp",
                 "(x<[1]) * ([0]*(TMath::Exp(-0.5*TMath::Power((x-[1])/[2],2)) + TMath::Exp((x-[1])/[3])*(1. - TMath::Exp(-0.5*TMath::Power((x-[1])/[2],2))))) + \
                 (x>=[1]) * ([0]*TMath::Exp(-0.5*TMath::Power((x-[1])/[2],2)))", 0,1);
            #self.f1total = TF1("GaussExpLinear", 
            #    "(x<[1]) * ([0]*(TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))))+[4]+[5]*x)+\
            #    (x>=[1])*([0]*TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))+[4]+[5]*x)", 0,1);
        self.func = "GaussExp"
        # elif func == "gausexpquadratic":
        #     self.f1total = TF1("GaussExpQuadratic", 
        #        "(x<[1]) * ([0]*(TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))))+[4]+[5]*x +[6]*(x-[1])^2)+\
        #        (x>=[1])*([0]*TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))+[4]+[5]*x + [6]*x*x)", 0,1);
        #       # (x>=[1])*([0]*TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))+[4]+[5]*x + [6]*(x-[1])^2)", 0,1);
        #     self.func = "GaussExpQuadratic"
        #     #if bkg == "pol1":
        #     #    self.f1bkg   = TF1("f1bkg","pol1(0)", 0,1);
        #     #    self.f1total = TF1("GaussExpLinear", 
        #     #            "(x<[1]) * ([0]*(TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))))+[4]+[5]*x)+\
        #     #            (x>=[1])*([0]*TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))+[4]+[5]*x)", 0,1);
            
        #    # elif bkg == "pol2":
        #    #     self.f1bkg   = TF1("f1bkg","pol2(0)", 0,1);
        #    #      self.f1total = TF1("f1total","crystalball(0) + pol2(5)", 0,1);

        self.f1total.SetNpx(1000);
        print("initially, ", self.f1total.GetExpFormula(""));

    def set_xtitle(self, title):
        self.xtitle = title;

    def set_ytitle(self, title):
        self.ytitle = title;

    def calc_FWHM(self, histo, params, paramsErr):

        # calculate FWHM value
        tf1_fwhm     = TF1("tf1_fwhm",   "(x<[1]) * ([0]*(TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2))))))+\
            (x>=[1])*([0]*TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2))))", 0,1);
        for i in range(4):
            tf1_fwhm.SetParameter(i, params[i])
            tf1_fwhm.SetParError(i, paramsErr[i])
        tf1_fwhm.SetNpx(1000);

        maximum = tf1_fwhm.GetMaximum()
        maximum_x = tf1_fwhm.GetMaximumX()
        half_maximum = maximum / 2
        left_x = tf1_fwhm.GetX(half_maximum, 0, maximum_x);
        right_x = tf1_fwhm.GetX(half_maximum, maximum_x, 1);
        FWHM = (right_x-left_x)/TMath.Sqrt(8*TMath.Log(2))  # in equivalents to sigma

        # Calculate PLUS-error with (FWHM+FWHM_err)
        tf1_fwhm_plus_err     = TF1("tf1_fwhm_plus_err",   "(x<[1]) * ([0]*(TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2))))))+\
            (x>=[1])*([0]*TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2))))", 0,1);
        for i in range(4):
            param = params[i] + paramsErr[i]
            tf1_fwhm_plus_err.SetParameter(i, param)

        maximum_plus = tf1_fwhm_plus_err.GetMaximum()
        maximum_x_plus = tf1_fwhm_plus_err.GetMaximumX()
        half_maximum_plus = maximum_plus / 2
        left_x_plus = tf1_fwhm_plus_err.GetX(half_maximum_plus, 0, maximum_x_plus);
        right_x_plus = tf1_fwhm_plus_err.GetX(half_maximum_plus, maximum_x_plus, 1);
        FWHM_plus = right_x_plus-left_x_plus

        # Calculate MINUS-error with (FWHM-FWHM_err)
        tf1_fwhm_minus_err     = TF1("tf1_fwhm_plus_err",   "(x<[1]) * ([0]*(TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2))))))+\
            (x>=[1])*([0]*TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2))))", 0,1);
        for i in range(4):
            param = params[i] - paramsErr[i]
            tf1_fwhm_minus_err.SetParameter(i, param)

        maximum_minus = tf1_fwhm_minus_err.GetMaximum()
        maximum_x_minus = tf1_fwhm_minus_err.GetMaximumX()
        half_maximum_minus = maximum_minus / 2
        left_x_minus = tf1_fwhm_minus_err.GetX(half_maximum_minus, 0, maximum_x_minus);
        right_x_minus = tf1_fwhm_minus_err.GetX(half_maximum_minus, maximum_x_minus, 1);
        FWHM_minus = right_x_minus-left_x_minus

        # Calculate TOTAL error
        FWHM_err = TMath.Sqrt((FWHM_plus-FWHM)**2 + (FWHM_minus-FWHM)**2)/TMath.Sqrt(8*TMath.Log(2))/2  # in equivalents to sigma and divided by 2 for error propagation

        return tf1_fwhm, FWHM, FWHM_err
    

    def set_yield_range(self, yield_min, yield_max):
        self.yield_min = yield_min;
        self.yield_max = yield_max;
    
    def calculate_raw_yield_MC(self, histo, params):
        integral_min = histo.FindBin(params[1] - self.yield_min)
        integral_max = histo.FindBin(params[1] + self.yield_max)
        error_integral = c_double(0.0)
        integral_histo = histo.IntegralAndError(integral_min, integral_max, error_integral)

        # # subtract linear background
        # tf1_linear     = TF1("tf1_linear","[4]+[5]*x", 0,1);
        # tf1_linear.SetParameter(5, params[5])
        # tf1_linear.SetParameter(4, params[4])
        # error_linear = c_double(0.0)
        # linear_integral = tf1_linear.Integral((params[1] - 0.035), (params[1] + 0.01))
        # # What do I do about the error of the linear fit?

        #error = TMath.Sqrt(error_integral**2)# + error_linear**2)

        raw_yield = integral_histo #- linear_integral
        print("raw yield", raw_yield, integral_histo)#, linear_integral)
        return raw_yield, error_integral

    def get_nev(self):
        h1ev   = self.list_ev_ss.FindObject("hCollisionCounter");
        nev = h1ev.GetBinContent(4);
        return nev

    #______________________________________________________________________
    def analyze_ptspectrum(self):
        print(sys._getframe().f_code.co_name);
        outlist = THashList();
        outlist.SetName("outlist");
        print(self.arr_pt);
        npt = len(self.arr_pt);

        ########################################
        #   Collision Counter, h2same, h2mix   #
        ########################################
        h1ev   = self.list_ev_ss.FindObject("hCollisionCounter");
        if self.list_pair_ss_cut:
            h2mc_help = self.list_pair_ss_cut.FindObject("nocut")
            if h2mc_help:
                h2mc = h2mc_help.FindObject("hMggPt_Pi0_Primary").Clone("h2mc");
            else:
                print("Object 'nocut' not found in list_pair_ss_cut.")
        else:
            print("list_pair_ss_cut is not valid.")


        h2mc.Sumw2();
        h2mc.SetDirectory(0);
        if self.meson == "pi0":      
            h2mc.RebinX(2);
        elif self.meson == "eta":
            h2mc.RebinX(4);
        h2mc.Sumw2();
        nev = h1ev.GetBinContent(4);
        print("NEV = ", nev)

        outlist.Add(h1ev);
        outlist.Add(h2mc);

        ########################################
        #           h1parameter plots          #
        ########################################
    
        h1amplitude     = TH1F("h1amplitude_param",   "amplitude",                npt-1, self.arr_pt);
        h1mean          = TH1F("h1mean_param" ,       "mean",                     npt-1, self.arr_pt);
        h1sigma         = TH1F("h1sigma_param",       "sigma",                    npt-1, self.arr_pt);
        h1exponential   = TH1F("h1exponential_param", "lambda",                   npt-1, self.arr_pt);
        h1FWHM          = TH1F("h1fwhm_param",        "fwhm/#sqrt{8ln(2)}",       npt-1, self.arr_pt);      
        h1yield         = TH1F("h1yield",             "raw yield",                npt-1, self.arr_pt);       

        h1amplitude.SetXTitle("#it{p}_{T} (GeV/#it{c})");
        h1amplitude.SetYTitle("amplitude of fit");
        h1mean.SetXTitle("#it{p}_{T} (GeV/#it{c})");
        h1mean.SetYTitle("peak mean (GeV/#it{c}^{2})");
        h1sigma.SetXTitle("#it{p}_{T} (GeV/#it{c})");
        h1sigma.SetYTitle("peak sigma (GeV/#it{c}^{2})");
        h1exponential.SetXTitle("#it{p}_{T} (GeV/#it{c})");
        h1exponential.SetYTitle("#exponential coeff. of fit");
        h1FWHM.SetXTitle("#it{p}_{T} (GeV/#it{c})");
        h1FWHM.SetYTitle("fwhm/#sqrt{8ln(2)}");
        h1yield.SetXTitle("#it{p}_{T} (GeV/#it{c})");
        h1yield.SetYTitle("raw yield");        

        ########################################
        #         loop over pt slices          #
        ########################################

        for i in range(0, npt-1):
            pt1 = self.arr_pt[i];
            pt2 = self.arr_pt[i+1];
            print("pt1 = ", pt1, " pt2 = ", pt2, "npt = ", i)

            h1mc = slice_histogram(h2mc, pt1, pt2, "x", False); # implemented "e" (error calc) in ProjectionX/ProjectionY in the used function
            #h1mc.Sumw2();
            h1mc.SetName("h1mgg_mc_pt{0}".format(i));
            h1mc.SetTitle("m_{{#gamma#gamma}}^{{mc}}, {0:2.1f} < #it{{p}}_{{T,#gamma#gamma}} < {1:2.1f} GeV/#it{{c}}".format(pt1, pt2));
            h1mc.RebinX(2);
            h1mc.Scale(1./nev);
            h1mc.Sumw2();
            h1mc.SetDirectory(0);
    
            bin_min = h1mc.FindBin(self.integral_min);
            bin_max = h1mc.FindBin(self.integral_max);
            print("Integral boundaries:", " integral_min:", self.integral_min, ", integral_max:", self.integral_max);
            integral_mc = h1mc.Integral(bin_min, bin_max);
            integral_complete = h1mc.Integral()
            print("mc integral:" ," bin_min:", bin_min, ", bin_max:", bin_max, ", integral in range:", integral_mc,", integral complete:", integral_complete);

            h1subtracted = h1mc
            h1subtracted.SetTitle("m_{{#gamma#gamma}}^{{sub.}}, {0:2.1f} < #it{{p}}_{{T,#gamma#gamma}} < {1:2.1f} GeV/#it{{c}}".format(pt1, pt2));
            h1subtracted.GetXaxis().SetRangeUser(0, 0.3);
            # h1subtracted.Scale(1./nev);
            # # h1subtracted.Sumw2();
            h1subtracted .SetDirectory(0);
            # h1subtracted.RebinX(2);
            # h1subtracted.Sumw2();

            ###########################
            #         FITTING         #
            ###########################
            height = 1.0;
            mean_init = 0.135;
            sigma_init = 0.020;
            if "pi0" in self.meson:
                mean_init = 0.135;
                sigma_init = 0.005 #0.020;
                lambda_init = 0.01
            elif "eta" in self.meson:
                mean_init = 0.548;
                sigma_init = 0.012;
            bin_mean = h1subtracted.FindBin(mean_init);
            height = h1subtracted.GetBinContent(bin_mean)# - 1.0;

            fit_min = self.fit_min
            fit_max = self.fit_max
            if pt1 <= 1:
                fit_min = 0.06

            print("initial height = ", height);

            f1total = self.f1total.Clone("f1total_pt{0}".format(i));
            f1total.SetParameter(0,height); # amplitude
            f1total.SetParameter(1,mean_init); # mean
            f1total.SetParameter(2,sigma_init); # sigma
            f1total.SetParameter(3,sigma_init); # exponential 
            f1total.SetParameter(4,50); # offset
            f1total.SetParameter(5,-100); # linear
            f1total.SetParameter(6,-100); # quadratic
            f1total.SetParLimits(0, 0.8*height, 1.8*height); #0.8*height, 1.8*height);
            f1total.SetParLimits(1, 0.125,0.138);#0.14); #mean_init - 2 * sigma_init, mean_init + 2 * sigma_init); # mean +/- 3 sigma
            f1total.SetParLimits(2, 0.5 * sigma_init, 1.5*sigma_init); # sigma in range 0.5 to 2 sigma
            f1total.SetParLimits(3,0.5 * lambda_init, 3*lambda_init);#default values #change to lambda
            # f1total.SetParLimits(4,-50, +100);
            # f1total.SetParLimits(5,-50, +50);
            # f1total.SetParLimits(6,-2000, -0.001);
            h1fit = h1subtracted.Clone("h1mgg_fitted_pt{0}".format(i));
            h1fit.Fit(f1total,"QRME","",fit_min, fit_max);

            amplitude       = f1total.GetParameter(0);
            amplitude_err   = f1total.GetParError(0)
            mean            = f1total.GetParameter(1);
            mean_err        = f1total.GetParError(1);
            sigma           = f1total.GetParameter(2);
            sigma_err       = f1total.GetParError(2);
            exponential     = f1total.GetParameter(3); 
            exponential_err = f1total.GetParError(3);
            offset          = f1total.GetParameter(4); 
            offset_err      = f1total.GetParError(4);
            linear          = f1total.GetParameter(5); 
            linear_err      = f1total.GetParError(5);
            if self.func == "GaussExpQuadratic":
                quadratic = f1total.GetParameter(6);
                quadratic_err = f1total.GetParError(6);

            params = [amplitude, mean, sigma, exponential, offset, linear]
            params_err = [amplitude_err, mean_err, sigma_err, exponential_err, offset_err, linear_err]
            tf1_fwhm, FWHM, FWHM_err  = self.calc_FWHM(f1total,params, params_err )
            #tf1_fwhm.SetName("f1fwhm_pt{0}".format(i));
            #tf1_fwhm.SetTitle("fwhm/ #sqrt{8ln(2)}")
            #outlist.Add(tf1_fwhm)

            raw_yield, error_raw_yield = self.calculate_raw_yield_MC(h1subtracted, params)
            print("HERE", raw_yield, error_raw_yield)

            # Only add data points if amplitude bigger than 20
            if amplitude >= 20/1e8:
                h1amplitude.SetBinContent(i+1, amplitude)
                h1amplitude.SetBinError(i+1, amplitude_err)
                h1mean.SetBinContent(i+1, mean);
                h1mean.SetBinError(i+1, mean_err);
                h1sigma.SetBinContent(i+1, sigma);
                h1sigma.SetBinError(i+1, sigma_err);
                h1exponential.SetBinContent(i+1, exponential);
                h1exponential.SetBinError(i+1, exponential_err);
                h1FWHM.SetBinContent(i+1, FWHM)
                h1FWHM.SetBinError(i+1, FWHM_err);
            h1yield.SetBinContent(i+1, raw_yield)
            h1yield.SetBinError(i+1, error_raw_yield)
            outlist.Add(h1mc);
            outlist.Add(h1subtracted);
            outlist.Add(h1fit);
            outlist.Add(f1total);
            print("\n \n \n")

        outlist.Add(h1amplitude);
        outlist.Add(h1mean);
        outlist.Add(h1exponential);
        outlist.Add(h1FWHM);
        outlist.Add(h1yield)
        outlist.Add(h1sigma);
        del h1amplitude
        del h1mean
        del tf1_fwhm
        del h1ev
        del h1exponential
        del h1FWHM
        del h1sigma
        del h1subtracted
        del h1fit
        del f1total
        del h1yield

        return outlist;
#___________________________________________________________________
if __name__ == "__main__":
    arr_pt = np.array([0.2,0.3,0.4,0.5], dtype=float);
    #ana = PairAnalyzer("pi0", "AnalysisResults_HL_75289.root", "pi0eta-to-gammagamma", "PCMPCM", "qc_qc", arr_pt);
    ana = PairAnalyzerMC("pi0", "/Users/alicamarieenderich/20230727_Analysis_Results_from_Daiki/AnalysisResults_HL_106682.root", "pi0eta-to-gammagamma");
    del ana;

    f1total = TF1("GaussExpLinear", 
            "(x<[1]) * ([0]*(TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))))+[4]+[5]*x)+\
            (x>=[1])*([0]*TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))+[4]+[5]*x)", 0,1);
    f1bkg = TF1("f1bkg","pol1(0)",0,1);
    f1sig = TF1("fGaussExp",
                "(x<[1]) * ([0]*(TMath::Exp(-0.5*TMath::Power((x-[1])/[2],2)) + TMath::Exp((x-[1])/[3])*(1. - TMath::Exp(-0.5*TMath::Power((x-[1])/[2],2))))) + \
                (x>=[1]) * ([0]*TMath::Exp(-0.5*TMath::Power((x-[1])/[2],2)))", 0,1);
    print("initially, ", f1total.GetExpFormula("p"));

    f1total.SetParameters(1, 0.135, 0.005, 0.6, 1, 1 , 1);
    f1total.Draw();
    print("later, ", f1total.GetExpFormula("p"));