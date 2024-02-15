# This code is based on the code provided by Daiki Sekihata in his repository
# https://github.com/dsekihat/photon_sw_run3 for photon analysis in ALICE Run 3


import numpy as np
import math
import ROOT
import datetime
import os
import yaml
from ROOT import TFile
from calc_Omega_isospin import calc_Omega_isospin
from calc_omega_wire import calc_omega_wire
from calc_ratio import calc_ratio
import os, sys, shutil
import math
import argparse
import numpy as np
import ctypes
import datetime
import yaml
import ROOT
ROOT.gROOT.SetBatch(True);
from ROOT import TFile, THashList, TF1
from ROOT import TFile, TDirectory, THashList, TH1F, TH2F, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TLine, TPython
from ROOT import gStyle, gROOT, gSystem
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan
gStyle.SetOptStat(0);
gStyle.SetOptTitle(0);

cutname = "qc" 
period_mc = "LHC23d1k";
period_data = "LHC22f"
suffix = "AnyTrack";
#filename_mc = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_LHC23d1k_125889.root"
#filename_data = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_LHC22f4_new_125184.root"
# filename_data = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_147811_LHC22_pass4_lowIR.root"
# filename_mc = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_147812_LHC23d1k.root"

filename_data = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_155756_LHC22f_pass4.root"
# filename_mc = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_147812_LHC23d1k.root"
filename_mc = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_155278_LHC23d1k.root"

config_file = "config_pp_13.6TeV_LHC22f_material.yml"
with open(config_file, "r", encoding="utf-8") as config_yml:
    config = yaml.safe_load(config_yml)
date = "this_thesis" # datetime.date.today().strftime("%Y%m%d");
folder = "/Users/alicamarieenderich/{0}_calibration_weights/".format(date);  
os.makedirs(folder, exist_ok=True);

for type in ["data", "mc"]:
    if type == "data": 
        file = filename_data;
    elif type == "mc":
        file = filename_mc;
    NGamma = calc_Omega_isospin(config, filename_data, filename_mc, cutname, folder, period_data, period_mc, suffix);
    NGamma.run(file, type, date);

# NGamma = calc_Omega_isospin(config, filename_data, filename_mc, cutname, folder, period_data, period_mc, suffix);
arr_pt_cuts = config["common"]["pt_cuts"];

file_isospin = "/Users/alicamarieenderich/this_thesis_calibration_weights/this_thesis_calibration_weights_pp_13.6TeV_LHC22f_and_LHC23d1k_AnyTrack_new.root"

# for ipt in range(len(arr_pt_cuts)):
#             NGamma.draw_NGamma_isospin(file, ipt, date)
NGamma.draw_NGamma_isospin_combined(file_isospin, arr_pt_cuts, date)

# filename_data = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_147811_LHC22_pass4_lowIR.root"
# filename_mc = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_147812_LHC23d1k.root"
for type in ["data", "mc"]:
    if type == "data": 
        file = filename_data;
    elif type == "mc":
        file = filename_mc;
    NGamma_wire = calc_omega_wire(config, filename_data, filename_mc, cutname, folder, period_data, period_mc, suffix)
    NGamma_wire.run(file, type, date)

file_wire = "/Users/alicamarieenderich/this_thesis_calibration_weights/this_thesis_calibration_weights_wire_pp_13.6TeV_LHC22f_and_LHC23d1k_AnyTrack_new.root"

# for ipt in range(len(arr_pt_cuts)):
#             NGamma_wire.draw_NGamma_wire(file, ipt, date)
NGamma_wire.draw_NGamma_wire_combined(file_wire, arr_pt_cuts, date)

Ratio = calc_ratio(config, filename_data, filename_mc, cutname, folder, period_data, period_mc, suffix);
Ratio.draw_ratio_weights_combined(file_isospin, file_wire, arr_pt_cuts, date)