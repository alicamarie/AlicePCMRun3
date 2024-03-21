# This code is based on the code provided by Daiki Sekihata in his repository
# https://github.com/dsekihat/photon_sw_run3 for photon analysis in ALICE Run 3

import os, sys, shutil
import numpy as np
import math
import ctypes
import ROOT
from ROOT import TFile, TDirectory, THashList, TF1, TH1F, TCanvas, TMath, TH1D
from ctypes import *

class PairAnalyzer:
    def __init__(self):
        print("default constructor is called");
    def __init__(self, meson, filename, dirname, type):
        print("target meson = {0} , filename = {1} , dirname = {2}".format(meson, filename, dirname));
        self.meson = meson;
        self.type = type
        self.rootfile = TFile.Open(filename, "READ");
        self.rootdir = self.rootfile.Get(dirname);
        self.list_ev = self.rootdir.Get("Event");
        self.list_pair = self.rootdir.Get("Pair");
        self.arr_pt = np.array([0,1,2,3,4,5], dtype=float);
        self.f1total = TF1("GaussExpLinear", 
               "(x<[1]) * ([0]*(TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))))+[4]+[5]*x)+\
               (x>=[1])*([0]*TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))+[4]+[5]*x)", 0,1);

        self.f1total.SetNpx(1000);
        self.fit_min = 0.04;
        self.fit_max = 0.24;
        self.integral_min = 0.18;
        self.integral_max = 0.25;

        self.yield_min = 0.035;
        self.yield_max= 0.01;
        self.xtitle = "#it{m}_{#gamma#gamma} (GeV/#it{c}^{2})";
        self.ytitle = "#it{p}_{T,#gamma#gamma} (GeV/#it{c}^{2})";
        self.initial_params = np.array([0,1,2,3,4,5], dtype=float)
        self.initial_lower_lim = np.array([0,1,2,3,4,5], dtype=float)
        self.initial_upper_lim = np.array([0,1,2,3,4,5], dtype=float)
    def __del__(self):
        if self.rootfile.IsOpen():
            print("close input root file.");
            self.rootfile.Close();

    def slice_histogram(self, h2,x0,x1,axis,isdiff):
        h1 = 0;
        delta = 1e-6;
        if "x" in axis.lower():
            bin0 = h2.GetYaxis().FindBin(x0 + delta);
            bin1 = h2.GetYaxis().FindBin(x1 - delta);
            h1 = h2.ProjectionX("h1prjx_{0}".format(h2.GetName()),bin0,bin1,"e");
        elif "y" in axis.lower():
            bin0 = h2.GetXaxis().FindBin(x0 + delta);
            bin1 = h2.GetXaxis().FindBin(x1 - delta);
            h1 = h2.ProjectionY("h1prjy_{0}".format(h2.GetName()),bin0,bin1,"e"); # "e" for error calculation

        if isdiff and h1.Class() == TH1D.Class(): #do you want differential histogram? e.g. dN/dpT, dN/dm.
            h1.Scale(1.,"width");
        return h1;

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
        self.list_pair_ss_cut = self.list_pair_ss.FindObject(cutname);

    def set_fit_range(self, fit_min, fit_max):
        self.fit_min = fit_min;
        self.fit_max = fit_max;
    
    def set_integral_range(self, integral_min, integral_max):
        self.integral_min = integral_min;
        self.integral_max = integral_max;

    def set_fit_function(self, type, func):
        if type == "data":
            if func == "gausexplinear":
                self.f1total = TF1("GaussExpLinear", 
                "(x<[1]) * ([0]*(TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))))+[4]+[5]*x)+\
                (x>=[1])*([0]*TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))+[4]+[5]*x)", 0.,1.);
                self.func = "GaussExpLinear"
            elif func == "gausexpquadratic":
                self.f1total = TF1("GaussExpQuadratic", 
                "(x<[1]) * ([0]*(TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))))+[4]+[5]*x +[6]*(x-[1])^2)+\
                (x>=[1])*([0]*TMath::Exp(-0.5*(TMath::Power((x-[1])/[2],2)))+[4]+[5]*x + [6]*x*x)", 0,1);
                self.func = "GaussExpQuadratic"

        elif type == "mc":
            self.f1total = TF1("fGaussExp",
                 "(x<[1]) * ([0]*(TMath::Exp(-0.5*TMath::Power((x-[1])/[2],2)) + TMath::Exp((x-[1])/[3])*(1. - TMath::Exp(-0.5*TMath::Power((x-[1])/[2],2))))) + \
                 (x>=[1]) * ([0]*TMath::Exp(-0.5*TMath::Power((x-[1])/[2],2)))", 0,1);
            self.func = "GaussExp"

        self.f1total.SetNpx(1000);
    
    def set_xtitle(self, title):
        self.xtitle = title;

    def set_ytitle(self, title):
        self.ytitle = title;

    def calc_FWHM(self, histo, params, paramsErr):
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
        FWHM_plus = (right_x_plus-left_x_plus)/TMath.Sqrt(8*TMath.Log(2))

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
        FWHM_minus = (right_x_minus-left_x_minus)/TMath.Sqrt(8*TMath.Log(2))

        # Calculate TOTAL error
        FWHM_err = max(abs(FWHM_plus-FWHM), abs(FWHM_minus-FWHM))
        FWHM_err = FWHM_err/2
        print("FWHM calc: ", FWHM*1e3, "+/-", FWHM_err*1e3, ";", FWHM_plus*1e3, FWHM_minus*1e3)
        return tf1_fwhm, FWHM, FWHM_err
    
    def calc_FWHM_MC(self, histo, params, paramsErr):
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
        FWHM_plus = (right_x_plus-left_x_plus)/TMath.Sqrt(8*TMath.Log(2))

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
        FWHM_minus = (right_x_minus-left_x_minus)/TMath.Sqrt(8*TMath.Log(2))

        # Calculate TOTAL error
        FWHM_err = max(abs(FWHM_plus-FWHM), abs(FWHM_minus-FWHM))
        FWHM_err = FWHM_err/2  # in equivalents to sigma and divided by 2 for error propagation
        return tf1_fwhm, FWHM, FWHM_err

    def set_yield_range(self, yield_min, yield_max):
        self.yield_min = yield_min;
        self.yield_max = yield_max;
    
    def calculate_raw_yield(self, histo, params, params_err, binwidth, nev):
        integral_min = histo.FindBin(params[1] - self.yield_max)
        integral_max = histo.FindBin(params[1] + self.yield_max)
        error_integral = c_double(0.0)
        integral_histo = histo.IntegralAndError(integral_min, integral_max, error_integral) #without "width"
        # subtract linear background
        tf1_linear     = TF1("tf1_linear","[4]+[5]*x", 0,1);
        tf1_linear.SetParameter(5, params[5])
        tf1_linear.SetParameter(4, params[4])
        lin_integral_min = params[1] - self.yield_max
        lin_integral_max = params[1] + self.yield_max
        linear_integral = params[4]*(lin_integral_max-lin_integral_min)+0.5*params[5]*(lin_integral_max*lin_integral_max-lin_integral_min*lin_integral_min);
        error_linear = np.sqrt( ((lin_integral_max-lin_integral_min)*params_err[4])**2 + 
                               (0.5*(lin_integral_max*lin_integral_max-lin_integral_min*lin_integral_min)*params_err[5])**2)
        error = TMath.Sqrt(error_integral.value**2 + error_linear**2)
        error = error / nev
        raw_yield = (integral_histo - linear_integral)/binwidth/nev
        return raw_yield, error
    
    def calculate_raw_yield_MC(self, histo, params, params_err, binwidth, nev):
        integral_min = histo.FindBin(params[1] - self.yield_min) #
        integral_max = histo.FindBin(params[1] + self.yield_max)
        error_integral = c_double(0.0)
        integral_histo = histo.IntegralAndError(integral_min, integral_max, error_integral) #without "width"
        error = error_integral.value / nev
        raw_yield = (integral_histo)/binwidth/nev
        print("raw yield", raw_yield, integral_histo)
        return raw_yield, error
    
    def get_nev(self):
        h1ev   = self.list_ev_ss.FindObject("hCollisionCounter");
        nev = h1ev.GetBinContent(4);
        return nev
    
    def set_fit_params(self, params, params_lower_lim, params_upper_lim):
        self.initial_params = params
        self.initial_lower_lim = params_lower_lim
        self.initial_upper_lim = params_upper_lim
    
    def set_2D_histograms(self, outlist):

        ########################################
        #   Collision Counter, h2same, h2mix   #
        ########################################
        h1ev   = self.list_ev_ss.FindObject("hCollisionCounter");
        if self.list_pair_ss_cut:
            h2same_help = self.list_pair_ss_cut.FindObject("nocut")
            if h2same_help:
                h2same = h2same_help.FindObject("hMggPt_Same").Clone("h2same");
            else:
                print("Object 'nocut' not found in list_pair_ss_cut.")
        else:
            print("list_pair_ss_cut is not valid.")
        if self.list_pair_ss_cut:
            h2mix_help = self.list_pair_ss_cut.FindObject("nocut")
            if h2mix_help:
                h2mix  = h2mix_help.FindObject("hMggPt_Mixed").Clone("h2mix");
            else:
                print("Object 'nocut' not found in list_pair_ss_cut.")
        else:
            print("list_pair_ss_cut is not valid.")

        h2same.Sumw2();
        h2mix .Sumw2();
        h2same.SetDirectory(0);
        h2mix .SetDirectory(0);
        if self.meson == "pi0":      
            h2same.RebinX(2);
            h2mix .RebinX(2);
        elif self.meson == "eta":
            h2same.RebinX(4);
            h2mix .RebinX(4);
        h2same.Sumw2();
        h2mix .Sumw2();
        nev = h1ev.GetBinContent(4);

        outlist.Add(h1ev);
        outlist.Add(h2same);
        outlist.Add(h2mix);
        return h2same, h2mix, h1ev, nev

    def set_2D_histograms_MC(self, outlist):
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
        return h2mc, h1ev, nev
    
    def slice_histo_and_subtract_background(self, i, h2same, h2mix, pt1, pt2, nev):
        h1same = self.slice_histogram(h2same, pt1, pt2, "x", False); # implemented "e" (error calc) in ProjectionX/ProjectionY in the used function
        #h1same.Sumw2();
        h1same.SetName("h1mgg_same_pt{0}".format(i));
        h1same.SetTitle("m_{{#gamma#gamma}}^{{same}}, {0:2.1f} < #it{{p}}_{{T,#gamma#gamma}} < {1:2.1f} GeV/#it{{c}}".format(pt1, pt2));
        h1mix  = self.slice_histogram(h2mix , pt1, pt2, "x", False);
        #h1mix.Sumw2();
        h1mix.SetTitle("m_{{#gamma#gamma}}^{{mix}}, {0:2.1f} < #it{{p}}_{{T,#gamma#gamma}} < {1:2.1f} GeV/#it{{c}}".format(pt1, pt2));
        h1mix.SetName("h1mgg_mix_pt{0}".format(i));
        h1same.RebinX(2);
        # h1same.Scale(1./nev);
        h1same.Sumw2();
        h1mix.RebinX(2);
        # h1mix.Scale(1./nev);
        h1mix.Sumw2();
        h1same.SetDirectory(0);
        h1mix .SetDirectory(0);

        bin_min_same = h1same.FindBin(self.integral_min);
        bin_max_same = h1same.FindBin(self.integral_max);
        integral_same_event = h1same.Integral(bin_min_same, bin_max_same);
        integral_complete_same = h1same.Integral()
        bin_min_mix = h1mix.FindBin(self.integral_min);
        bin_max_mix = h1mix.FindBin(self.integral_max);
        integral_mixed_event = h1mix.Integral(bin_min_mix, bin_max_mix);
        integral_complete_mix = h1mix.Integral()
        h1mix_scaled = h1mix.Clone("h1mgg_mix_scaled_pt{0}".format(i));
        if integral_mixed_event != 0:
            h1mix_scaled.Scale(integral_same_event / integral_mixed_event)
        else:
            print("Error: Division by zero. integral_mixed_event is zero.")
        #h1mix_scaled.Sumw2();
        h1mix_scaled.SetDirectory(0);
        h1subtracted = h1same.Clone("h1mgg_pt{0}".format(i));
        h1subtracted.Add(h1mix_scaled, -1);
        #h1subtracted.Sumw2();
        h1subtracted.SetTitle("m_{{#gamma#gamma}}^{{sub.}}, {0:2.1f} < #it{{p}}_{{T,#gamma#gamma}} < {1:2.1f} GeV/#it{{c}}".format(pt1, pt2));
        h1subtracted.GetXaxis().SetRangeUser(0, 0.3);
        h1subtracted .SetDirectory(0);

        return h1subtracted, h1mix_scaled, h1same, h1mix
    
    def slice_histo(self,i, h2mc, pt1, pt2, nev):
        h1mc = self.slice_histogram(h2mc, pt1, pt2, "x", False); # implemented "e" (error calc) in ProjectionX/ProjectionY in the used function
        #h1mc.Sumw2();
        h1mc.SetName("h1mgg_pt{0}".format(i));
        h1mc.SetTitle("m_{{#gamma#gamma}}^{{mc}}, {0:2.1f} < #it{{p}}_{{T,#gamma#gamma}} < {1:2.1f} GeV/#it{{c}}".format(pt1, pt2));
        h1mc.RebinX(2);
        # h1mc.Scale(1./nev);
        h1mc.Sumw2();
        h1mc.SetDirectory(0);

        bin_min = h1mc.FindBin(self.integral_min);
        bin_max = h1mc.FindBin(self.integral_max);
        integral_mc = h1mc.Integral(bin_min, bin_max);
        integral_complete = h1mc.Integral()
        return h1mc

    #______________________________________________________________________
    def analyze_ptspectrum(self):
        outlist = THashList();
        outlist.SetName("outlist");
        npt = len(self.arr_pt);

        if self.type == "data":
            h2same, h2mix, h1ev, nev = self.set_2D_histograms(outlist)
        elif self.type == "mc":
            h2mc, h1ev, nev = self.set_2D_histograms_MC(outlist)

        ########################################
        #           h1parameter plots          #
        ########################################
    
        h1amplitude     = TH1F("h1amplitude_param",   "amplitude",                npt-1, self.arr_pt);
        h1mean          = TH1F("h1mean_param" ,       "mean",                     npt-1, self.arr_pt);
        h1sigma         = TH1F("h1sigma_param",       "sigma",                    npt-1, self.arr_pt);
        h1exponential   = TH1F("h1exponential_param", "lambda",                   npt-1, self.arr_pt);
        h1offset        = TH1F("h1offset_param",      "offset",                   npt-1, self.arr_pt);
        h1linear        = TH1F("h1linear_param",      "linear coeff.",            npt-1, self.arr_pt);
        h1quadratic     = TH1F("h1quadratic_param",   "quadratic coeff.",         npt-1, self.arr_pt);
        h1FWHM          = TH1F("h1fwhm_param",        "fwhm/2.36",                npt-1, self.arr_pt);      
        h1yield         = TH1F("h1yield_param",       "raw yield",                npt-1, self.arr_pt);       

        h1amplitude.SetXTitle("#it{p}_{T} (GeV/#it{c})");
        h1amplitude.SetYTitle("amplitude of fit");
        h1mean.SetXTitle("#it{p}_{T} (GeV/#it{c})");
        h1mean.SetYTitle("peak mean (GeV/#it{c}^{2})");
        h1sigma.SetXTitle("#it{p}_{T} (GeV/#it{c})");
        h1sigma.SetYTitle("peak sigma (GeV/#it{c}^{2})");
        h1exponential.SetXTitle("#it{p}_{T} (GeV/#it{c})");
        h1exponential.SetYTitle("#exponential coeff. of fit");
        h1offset.SetXTitle("#it{p}_{T} (GeV/#it{c})");
        h1offset.SetYTitle("offset of fit");
        h1linear.SetXTitle("#it{p}_{T} (GeV/#it{c})");
        h1linear.SetYTitle("linear coeff. of fit");
        h1quadratic.SetXTitle("#it{p}_{T} (GeV/#it{c})");
        h1quadratic.SetYTitle("quadratic coeff. of fit");
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

            if self.type == "data":
                h1prepared, h1mix_scaled, h1same, h1mix = self.slice_histo_and_subtract_background(i, h2same, h2mix, pt1, pt2, nev)

            if self.type == "mc":
                h1prepared = self.slice_histo(i, h2mc, pt1, pt2, nev)

            ###########################
            #         FITTING         #
            ###########################
            height = 1.0;
            mean_init = 0.135;
            sigma_init = 0.020;
            if "pi0" in self.meson:
                mean_init = 0.135;
                sigma_init = 0.006 #0.020;
                lambda_init = 0.016
            elif "eta" in self.meson:
                mean_init = 0.548;
                sigma_init = 0.012;
            for j in range ( h1prepared.GetXaxis().FindBin(0.125), h1prepared.GetXaxis().FindBin(0.138)):
                    if h1prepared.GetBinContent(j) > height:
                        height = h1prepared.GetBinContent(j);
            fit_min = self.fit_min
            fit_max = self.fit_max

            f1total = self.f1total.Clone("f1total_pt{0}".format(i));
            param_init = [height, mean_init, sigma_init, lambda_init]
            for i_lim, (lower_limit, upper_limit) in enumerate(zip(self.initial_lower_lim, self.initial_upper_lim)):
                print(i_lim, type(lower_limit))
                f1total.SetParLimits(i_lim, lower_limit*param_init[i_lim], upper_limit*param_init[i_lim])
                
            if self.initial_params[0] == 0: # amplitude 0 -> did not yet assign values in config

                f1total.SetParameter(0,height); # amplitude
                f1total.SetParameter(1,mean_init); # mean
                f1total.SetParameter(2,sigma_init); # sigma
                f1total.SetParameter(3,lambda_init); # exponential 
                if self.type =="data":
                    f1total.SetParameter(4,50); # offset
                    f1total.SetParameter(5,-100); # linear
                    f1total.SetParameter(6,-100); # quadratic
            
                f1total.SetParLimits(0, 0.9*height, 1.2*height);
                f1total.SetParLimits(1, 0.125,0.138);
                f1total.SetParLimits(2, 0.9 * sigma_init, 1.5*sigma_init);
                f1total.SetParLimits(3,0.9 * lambda_init, 1.6*lambda_init);

            h1fit = h1prepared.Clone("h1mgg_fitted_pt{0}".format(i));
            h1fit.Fit(f1total,"RME","",fit_min, fit_max);

            amplitude       = f1total.GetParameter(0);
            amplitude_err   = f1total.GetParError(0)
            mean            = f1total.GetParameter(1);
            mean_err        = f1total.GetParError(1);
            sigma           = f1total.GetParameter(2);
            sigma_err       = f1total.GetParError(2);
            exponential     = f1total.GetParameter(3); 
            exponential_err = f1total.GetParError(3);
            if self.type == "data":
                offset          = f1total.GetParameter(4); 
                offset_err      = f1total.GetParError(4);
                linear          = f1total.GetParameter(5); 
                linear_err      = f1total.GetParError(5);
                if self.func == "GaussExpQuadratic":
                    quadratic = f1total.GetParameter(6);
                    quadratic_err = f1total.GetParError(6);

            if self.type == "data":
                params = [amplitude, mean, sigma, exponential, offset, linear]
                params_err = [amplitude_err, mean_err, sigma_err, exponential_err, offset_err, linear_err]
                tf1_fwhm, FWHM, FWHM_err  = self.calc_FWHM(f1total,params, params_err )
            elif self.type == "mc":
                params = [amplitude, mean, sigma, exponential]
                params_err = [amplitude_err, mean_err, sigma_err, exponential_err]
                tf1_fwhm, FWHM, FWHM_err  = self.calc_FWHM_MC(f1total,params, params_err )
                print("FWHM", FWHM, type(FWHM))

            tf1_fwhm.SetName("f1fwhm_pt{0}".format(i));

            if self.type == "data":
                raw_yield, error_raw_yield = self.calculate_raw_yield(h1prepared, params, params_err, (pt2-pt1), nev)
            elif self.type == "mc":
                raw_yield, error_raw_yield = self.calculate_raw_yield_MC(h1prepared, params, params_err, (pt2-pt1), nev)
            # Only add data points if amplitude bigger than 20
            #if amplitude >= 20/1e8:
            h1amplitude.SetBinContent(i+1, amplitude)
            h1amplitude.SetBinError(i+1, amplitude_err)
            h1mean.SetBinContent(i+1, mean);
            h1mean.SetBinError(i+1, mean_err);
            h1sigma.SetBinContent(i+1, sigma);
            h1sigma.SetBinError(i+1, sigma_err);
            h1exponential.SetBinContent(i+1, exponential);
            h1exponential.SetBinError(i+1, exponential_err);
            if not math.isnan(FWHM):
                print("is not nan", FWHM)
                h1FWHM.SetBinContent(i+1, FWHM)
                h1FWHM.SetBinError(i+1, FWHM_err);
            else:
                print("is nan", FWHM)
            if self.type == "data":
                h1linear.SetBinContent(i+1, linear);
                h1linear.SetBinError(i+1, linear_err);
                h1offset.SetBinContent(i+1, offset);
                h1offset.SetBinError(i+1, offset_err);

                if self.func == "GaussExpQuadratic":
                        h1quadratic.SetBinContent(i+1, quadratic);
                        h1quadratic.SetBinError(i+1, quadratic_err);
            h1yield.SetBinContent(i+1, raw_yield)
            h1yield.SetBinError(i+1, error_raw_yield)

            if self.type == "data":
                outlist.Add(h1mix_scaled);
                outlist.Add(h1same);
                outlist.Add(h1mix);
            outlist.Add(h1prepared);
            outlist.Add(h1fit);
            outlist.Add(f1total);
            print("\n \n \n")

        outlist.Add(h1yield);
        outlist.Add(h1mean);
        outlist.Add(h1exponential);
        outlist.Add(h1sigma);
        outlist.Add(h1FWHM);
        outlist.Add(h1amplitude);
        outlist.Add(h1linear);
        outlist.Add(h1offset);
        outlist.Add(h1quadratic);

        del h1amplitude
        del h1mean
        del tf1_fwhm
        del h1ev
        del h1exponential
        del h1FWHM
        del h1sigma
        del h1linear
        del h1offset
        del h1quadratic
        if self.type == "data":
            del h1mix
            del h1same
        del h1prepared
        del h1fit
        del f1total
        del h1yield

        return outlist;