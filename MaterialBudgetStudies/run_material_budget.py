# This code is based on the code provided by Daiki Sekihata in his repository
# https://github.com/dsekihat/photon_sw_run3 for photon analysis in ALICE Run 3


import numpy as np
import math
import ROOT
import datetime
import os
import yaml
from ROOT import TFile
from draw_1D_integrated import draw_1D_integrated
from draw_1D_sliced_in_R import draw_1D_sliced_rxy
from draw_1D_sliced_in_eta import draw_1D_sliced_eta
from draw_1D_sliced_in_pt import draw_1D_sliced_pt
from make_2D_Eta_vs_Phi import make_2D_Eta_vs_Phi
from draw_2D_Eta_vs_Phi import draw_material_z_vs_phi
from draw_hPhotonRxy_MC_gen import draw_hPhotonRxy_mc_gen, draw_hPhotonRxy_mc_rec
from draw_hPhotonRxy_data import draw_hPhotonRxy_data
from draw_hPhotonRxy_vs_Phi import draw_hPhotonRxy_vs_Phi_mc_gen, draw_hPhotonRxy_vs_Phi_mc_rec
from draw_hPhoton_RZ import draw_material_RZ
from draw_AP_with_cuts import draw_AP_with_cuts

cutname = "qc"
period_mc = "LHC23d1k";
period_data = "LHC22f"
suffix = "AnyTrack";
filename_mc = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_155278_LHC23d1k.root"
# # filename_mc = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_LHC23d1k_125889.root"
# filename_data = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_148193_LHC22f_apass4.root"
# # filename_data = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_LHC22f4_new_125184.root"

filename_data = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_155756_LHC22f_pass4.root"
# filename_mc = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_147812_LHC23d1k.root"

config_file = "/Users/alicamarieenderich/202312_material_budget_code/config_pp_13.6TeV_LHC22f_material.yml"
with open(config_file, "r", encoding="utf-8") as config_yml:
    config = yaml.safe_load(config_yml)
date = "presentation"  #"this_thesis"# datetime.date.today().strftime("%Y%m%d");
folder = "/Users/alicamarieenderich/{0}_material_budget_plots/".format(date);  
os.makedirs(folder, exist_ok=True);

# for cuts in [True, False]:
#     for generated in [True, False]:
#         draw_integrated = draw_1D_integrated(filename_data, filename_mc, cutname, folder, period_data, suffix);
#         draw_integrated.draw_material_rxy(cuts, generated, date);
#         draw_integrated.draw_material_phi(cuts, generated, date);
#         draw_integrated.draw_material_eta(cuts, generated, date);
#         draw_integrated.draw_material_nch(cuts, date);

# for type in ["data", "mc"]:
#     if type == "data": 
#         file = filename_data;
#     elif type == "mc":
#         file = filename_mc;
#     draw_sliced = draw_1D_sliced_rxy(config, suffix, folder, period_data, period_mc);
#     draw_sliced.run(file, type, date);

# filename_mc_rxy = os.path.join(folder, "{0}_material_budget_dR_mc_{1}_{2}TeV_{3}{4}.root".format(date, config["common"]["system"], config["common"]["energy"], config["common"]["period"], suffix));
# filename_data_rxy = os.path.join(folder, "{0}_material_budget_dR_data_{1}_{2}TeV_{3}{4}.root".format(date, config["common"]["system"], config["common"]["energy"], config["common"]["period"], suffix));

# draw_sliced = draw_1D_sliced_rxy(config, suffix, folder, period_data, period_mc);
# for eta in range(6):
#     for r in range(6):
#         draw_sliced.draw_material_phi(filename_data_rxy, filename_mc_rxy, cutname, eta, r, date);

# for type in ["data", "mc"]:
#     if type == "data": 
#         file = filename_data
#     elif type == "mc":
#         file = filename_mc;
#     draw_sliced = draw_1D_sliced_eta(config, suffix, folder, period_data, period_mc, cutname);
#     draw_sliced.run(file, type, date);

# filename_mc_eta = os.path.join(folder, "{0}_material_budget_dEta_mc_{1}_{2}TeV_{3}{4}.root".format(date, config["common"]["system"], config["common"]["energy"], config["common"]["period"], suffix));
# filename_data_eta = os.path.join(folder, "{0}_material_budget_dEta_data_{1}_{2}TeV_{3}{4}.root".format(date, config["common"]["system"], config["common"]["energy"], config["common"]["period"], suffix));

# draw_sliced = draw_1D_sliced_eta(config, suffix, folder, period_data, period_mc, cutname);
# for eta in range(6): #6
#     for r in range(6): #6
#         draw_sliced.draw_material_eta(filename_data_eta, filename_mc_eta, r, cutname, date);

# for type in ["data", "mc"]:
#     if type == "data": 
#         file = filename_data
#     elif type == "mc":
#         file = filename_mc;
#     draw_sliced = draw_1D_sliced_pt(config, suffix, folder, period_data, period_mc, cutname);
#     draw_sliced.run(file, type, date);

# filename_mc_pt = os.path.join(folder, "{0}_material_budget_dPt_mc_{1}_{2}TeV_{3}{4}.root".format(date, config["common"]["system"], config["common"]["energy"], config["common"]["period"], suffix));
# filename_data_pt = os.path.join(folder, "{0}_material_budget_dPt_data_{1}_{2}TeV_{3}{4}.root".format(date, config["common"]["system"], config["common"]["energy"], config["common"]["period"], suffix));

# draw_sliced = draw_1D_sliced_pt(config, suffix, folder, period_data, period_mc, cutname);
# #for eta in range(6): #6
# # for r in range(6): #6
# draw_sliced.draw_material_pt_combined(filename_data_pt, filename_mc_pt, 1, cutname, date);
# # for eta in range(6): #6
# for r in range(6): #6
#     draw_sliced.draw_material_pt(filename_data_pt, filename_mc_pt, r, cutname, date);

# for type in ["data", "mc"]:
#     if type == "data": 
#         file = filename_data
#     elif type == "mc":
#         file = filename_mc;
#     make_eta_vs_phi = make_2D_Eta_vs_Phi(config, suffix, folder, period_data, period_mc, cutname);
#     make_eta_vs_phi.run(file, type, date);

# filename_mc_eta_vs_phi = os.path.join(folder, "{0}_material_budget_eta_vs_phi_mc_{1}_{2}TeV_{3}{4}.root".format(date, config["common"]["system"], config["common"]["energy"], config["common"]["period"], suffix));
# filename_data_eta_vs_phi = os.path.join(folder, "{0}_material_budget_eta_vs_phi_data_{1}_{2}TeV_{3}{4}.root".format(date, config["common"]["system"], config["common"]["energy"], config["common"]["period"], suffix));

# for r in range(6):
#     for log in ["log", ""]:
#         draw_material_z_vs_phi(folder, filename_data_eta_vs_phi, filename_mc_eta_vs_phi, cutname, period_data, r, suffix, log, date);

# for log in ["log", ""]:
#     for mctype in ["gen"]:
#         for zoom in [100]:
#             for circle in ["circle", ""]:
#                 draw_hPhotonRxy_mc_gen(suffix, log, mctype , zoom, circle, folder, date);
# log = "log"
# mctype = "gen"
# zoom = 100
# circle ="circle"
# draw_hPhotonRxy_mc_gen(suffix, "log", "gen" , 100, "", folder, date);
draw_hPhotonRxy_mc_gen(suffix, "log", "gen" , 100, "circle", folder, date);
# print("draw_hPhotonRxy_mc_gen done")
# draw_hPhotonRxy_mc_rec(suffix, "log", "rec", 100, "", folder, date)
# print("draw_hPhotonRxy_mc_rec done")

draw_hPhotonRxy_data(filename_data, filename_mc, cutname, suffix, folder, date);
print("draw_hPhotonRxy_data done")
outname         = os.path.join(folder,"{0}_PhotonPhivsRxy.root".format(date));
outfile         = TFile(outname, "RECREATE");
draw_hPhotonRxy_vs_Phi_mc_gen(suffix, "log", outfile, "gen", folder, date);
print("draw_hPhotonRxy_vs_Phi_mc_gen done")
draw_hPhotonRxy_vs_Phi_mc_rec(suffix, "log", outfile, "rec", folder, date); 
# draw_material_RZ(filename_data, filename_mc, suffix, 40, "wire", False, folder, date);
# draw_material_RZ(filename_data, filename_mc, suffix, 100, "complete_comparison", False,folder, date);
draw_material_RZ(filename_data, filename_mc, suffix, 100, "strcutures", True,folder, date);
# draw_AP_with_cuts(filename_data, filename_mc, cutname, suffix, folder, date)