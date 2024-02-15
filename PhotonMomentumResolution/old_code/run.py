# This code is based on the code provided by Daiki Sekihata in his repository
# https://github.com/dsekihat/photon_sw_run3 for photon analysis in ALICE Run 3


import os, sys, shutil
import math
import argparse
import numpy as np
import ctypes
import yaml
import ROOT
import datetime
ROOT.gROOT.SetBatch(True);
from ROOT import TFile, THashList, TF1, TString
#from analyze_pair import analyze_ptspectrum
from FitInvMassForPt import PairAnalyzer
from PlotInvMass import PlotInvMass
from PlotParameterHistoInvMass import PlotHistoInvMass
from PlotParametersCombined import PlotHistoParametersCombined
from PlotRawYield import PlotRawYieldInvMass
from FitInvMassForPtMC import PairAnalyzerMC

#_________________________________________________________________________________________
def run(filename, config, type, suffix, folder):
    # date = "this_thesis" #datetime.date.today().strftime("%Y - %m - %d");
    print(sys._getframe().f_code.co_name);
    arr_pt = np.array(config["common"]["pt_bin"],dtype=float);
    print("pT binning = ",arr_pt);
    print("type = ",type);
    print("input = ",filename);
    rootfile = TFile.Open(filename,"READ");
    meson = config["common"]["meson"];

    list_fit_func = config["common"]["fit_func"];

    list_fit_min = config["common"]["fit_min"];
    list_fit_max = config["common"]["fit_max"];
    if len(list_fit_min) != len(list_fit_max):
        return;

    list_integral_min = config["common"]["integral_min"];
    list_integral_max = config["common"]["integral_max"];
    if len(list_integral_min) != len(list_integral_max):
        return;

    list_yield_min = config["common"]["yield_min"];
    list_yield_max = config["common"]["yield_max"];
    if len(list_yield_min) != len(list_yield_max):
        return;

    nsys = len(config[type]['subsystems']);
    print(nsys); 

    meson = config["common"]["meson"];

    if config["common"]["do_ptspectrum"] == True:
        outname = os.path.join(folder, "{0}_{1}_{2}_{3}_ptspectrum_{4}_{5}TeV_{6}_{7}.root".format(date,period, meson, type, config["common"]["system"], config["common"]["energy"], config["common"]["period"], suffix));
        print("output file name = ",outname);
        outfile = TFile(outname,"RECREATE");

        if type == "mc":
            ana_pi0_mc = PairAnalyzerMC(meson, filename, "pi0eta-to-gammagamma-mc");
            ana_pi0_mc .set_arr_pt(arr_pt);
            #nev = ana_pi0_mc .get_nev()
            for isys in range(0,nsys):
                ssname = config[type]['subsystems'][isys]['name']; #subsystem name
                ana_pi0_mc .set_subsystem(ssname);
                ana_pi0_mc .set_xtitle("#it{m}_{#gamma#gamma} (GeV/#it{c}^{2})");
                ana_pi0_mc .set_ytitle("#it{p}_{T,#gamma#gamma} (GeV/#it{c}^{2})");
                print("analyze subsystem", ssname);
                cutnames = config[type]["subsystems"][isys]['cutnames']
                print("cutnames", cutnames); 
                nc = len(cutnames);
                outlist_ss = THashList();
                outlist_ss.SetName(ssname);
                outlist_ss.SetOwner(True);
                for ic in range(0,nc):
                    cutname = cutnames[ic];
                    ana_pi0_mc .set_cutname(cutname);
                    outlist_cut = THashList();
                    outlist_cut.SetName(cutname);
                    outlist_cut.SetOwner(True);
                    outlist_ss.Add(outlist_cut);
                    for ifunc in list_fit_func:
                        ana_pi0_mc .set_fit_function(ifunc);
                        outlist_func = THashList();
                        outlist_func.SetName(ifunc);
                        outlist_func.SetOwner(True);
                        outlist_cut.Add(outlist_func);
                        for ir in range(0, len(list_fit_min)):
                            fit_min = list_fit_min[ir];
                            fit_max = list_fit_max[ir];
                            ana_pi0_mc .set_fit_range(fit_min, fit_max);
                            integral_min = list_integral_min[0];
                            integral_max = list_integral_max[0];
                            ana_pi0_mc.set_integral_range(integral_min, integral_max);
                            yield_min = list_yield_min[0];
                            yield_max = list_yield_max[0];
                            ana_pi0_mc.set_yield_range(yield_min, yield_max);
                            outlist_fit_range = ana_pi0_mc .analyze_ptspectrum();
                            outlist_fit_range.SetName("fit_{0:3.2f}_{1:3.2f}_GeVc2".format(fit_min, fit_max));
                            outlist_func.Add(outlist_fit_range);
                outfile.WriteTObject(outlist_ss);
                outlist_ss.Clear();
            del ana_pi0_mc ;
        

            
            filename_plot = ""
            plot_pi0 = PlotInvMass(meson, outname, "pi0eta-to-gammagamma-mc");
            plot_pi0.set_arr_pt(arr_pt);
            for isys in range(0,nsys):
                ssname = config[type]['subsystems'][isys]['name']; #subsystem name
                plot_pi0.set_subsystem(ssname);
                print("plot subsystem", ssname);
                cutnames = config[type]["subsystems"][isys]['cutnames']
                print("cutnames", cutnames); 
                nc = len(cutnames);
                nfit = len(list_fit_func)
                list_parameters_linear_comparison = [];
                yield_list_lin = [];
                #list_parameters_quadratic_comparison = [];
                for ic in range(0,nc):
                    cutname = cutnames[ic];
                    print("cutname ", cutname)
                    plot_pi0.set_cutname(cutname);
                    for ifit in range(0, nfit): 
                        fitname = list_fit_func[ifit];
                        print("fitname ", fitname)
                        plot_pi0.set_fitname(fitname);
                        for ir in range(len(list_fit_min)):
                            fit_min = list_fit_min[ir];
                            fit_max = list_fit_max[ir];
                            plot_pi0.set_fit_range(fit_min, fit_max);
                            yield_min = list_yield_min[0];
                            yield_max = list_yield_max[0];
                            plot_pi0.set_yield_range(yield_min, yield_max);
                            plottingRange = [0., 0.3];
                            list_plot = plot_pi0.set_fit_list();
                            if ifit == 0:
                                yield_list_lin.append(list_plot.FindObject("h1yield"))
                            function_list = [];
                            histogram_list = [];
                            parameter_list = [];
                            same_list = [];
                            mixed_list = [];
                            it = ROOT.TIter(list_plot)
                            obj = it.Next()
                            while obj:
                                objName = obj.GetName()
                                if "param" in objName:
                                    parameter_list.append(obj)
                                # if "yield" in objName:
                                #     yield_list.append(obj)
                                obj = it.Next()

        # pdf output of all parameters for each combination of cut and fit                    
                            plot_histo = PlotHistoInvMass(meson, outname,"pi0eta-to-gammagamma")
                            outname_histo = os.path.join(folder, "{0}_{1}_InvMass_Fitparameters_{2}_{3}.pdf".format(date, period, cutname, fitname));
                            plot_histo.PlotHistoParameters(parameter_list, TString(outname_histo), "", "", period, 2,4,"#gamma#gamma", "data")

        # pdf output of all fitted histograms for each combination of cut and fit
                            for ipt in range(len(arr_pt)):
                                histo_ipt =list_plot.FindObject("h1mgg_mc_pt{0}".format(ipt));
                                histogram_list.append(histo_ipt)
                                function_ipt = list_plot.FindObject("f1total_pt{0}".format(ipt));
                                function_list.append(function_ipt)
                            output_name = os.path.join(folder, "{0}_{1}_InvMass_Overview_{2}_{3}.pdf".format(date, period, cutname, fitname))
                            plot_pi0.PlotInvMassInPtBins(histogram_list, function_list, TString(output_name), "", "", plottingRange, period, 3,3 , 0, len(arr_pt),len(arr_pt), "#gamma#gamma", "data")
        
        # # pdf output of all same and mixed scaled histograms for each cut
        #                     for ipt in range(len(arr_pt)):
        #                         histo_ipt_same =list_plot.FindObject("h1mgg_same_pt{0}".format(ipt));
        #                         same_list.append(histo_ipt_same)
        #                         histo_ipt_mixed = list_plot.FindObject("h1mgg_mix_scaled_pt{0}".format(ipt));
        #                         mixed_list.append(histo_ipt_mixed)
        #                     output_name = os.path.join(folder, "{0}_{1}_InvMass_Scaled_Mixed_{2}_{3}.pdf".format(date, period, cutname, fitname))
        #                     plot_pi0.PlotSameMixedInPtBins(same_list, mixed_list, TString(output_name), "", "", plottingRange, period, 3,3 , 0, len(arr_pt),len(arr_pt), "#gamma#gamma", "data")

        # pdf output of mass, amplitude and width for each fit and comparison of all cuts
                            if ifit ==0:
                                parameter_comparison = [];
                                for i_parameter in range(4):
                                    parameter_comparison.append(parameter_list[i_parameter]);
                                list_parameters_linear_comparison.append(parameter_comparison);
                            elif ifit ==1:
                                parameter_comparison = [];
                                for i_parameter in range(4):
                                    parameter_comparison.append(parameter_list[i_parameter]);
                                #.append(parameter_comparison);  
                          
                print(yield_list_lin, "HERE", len(yield_list_lin)) 
                plot_raw_yield = PlotRawYieldInvMass(meson, outname, "pi0eta-to-gammagamma-mc")     
                outname_raw_yield = os.path.join(folder, "{0}_{1}_InvMass_RawYield_{2}.pdf".format(date, period, fitname)); 
                plot_raw_yield.PlotHistoYield(yield_list_lin, TString(outname_raw_yield), "", "", period, 1,1, "#gamma#gamma", "data", cutnames)                   
                plot_parameters = PlotHistoParametersCombined(meson, outname,"pi0eta-to-gammagamma-mc")
                outname_histo_linear = os.path.join(folder, "{0}_{1}_InvMass_Parameters_Combined_{2}.pdf".format(date, period, fitname));
                #outname_histo_quadratic = os.path.join(folder, "InvMass_Parameters_quadratic.pdf");
                plot_parameters.PlotHistoParametersCombined(list_parameters_linear_comparison, TString(outname_histo_linear), "", "", period, 1, 5, "#gamma#gamma", "data", cutnames)
                #plot_parameters.PlotHistoParametersCombined(list_parameters_quadratic_comparison, TString(outname_histo_quadratic), "", "", date, 1, 5, "#gamma#gamma", "data", cutnames)

            del plot_pi0;
            del plot_histo;
            del parameter_list;
            del histogram_list;
            del function_list

        
            # remaining stuff from Daiki

            # return;
            # # for ic in range(0,nc):
            # #     outlist = analyze_ptspectrum_efficiency(rootfile,cutnames[ic],arr_mee,arr_ptee);
            # #     outlist.SetOwner(True);
            # #     outfile.WriteTObject(outlist);
            # #     outlist.Clear();
        else:
            ana_pi0 = PairAnalyzer(meson, filename, "pi0eta-to-gammagamma");
            ana_pi0.set_arr_pt(arr_pt);
            #nev = ana_pi0.get_nev()
            for isys in range(0,nsys):
                ssname = config[type]['subsystems'][isys]['name']; #subsystem name
                ana_pi0.set_subsystem(ssname);
                ana_pi0.set_xtitle("#it{m}_{#gamma#gamma} (GeV/#it{c}^{2})");
                ana_pi0.set_ytitle("#it{p}_{T,#gamma#gamma} (GeV/#it{c}^{2})");
                print("analyze subsystem", ssname);
                cutnames = config[type]["subsystems"][isys]['cutnames']
                print("cutnames", cutnames); 
                nc = len(cutnames);
                outlist_ss = THashList();
                outlist_ss.SetName(ssname);
                outlist_ss.SetOwner(True);
                for ic in range(0,nc):
                    cutname = cutnames[ic];
                    ana_pi0.set_cutname(cutname);
                    outlist_cut = THashList();
                    outlist_cut.SetName(cutname);
                    outlist_cut.SetOwner(True);
                    outlist_ss.Add(outlist_cut);
                    for ifunc in list_fit_func:
                        ana_pi0.set_fit_function(ifunc);
                        outlist_func = THashList();
                        outlist_func.SetName(ifunc);
                        outlist_func.SetOwner(True);
                        outlist_cut.Add(outlist_func);
                        for ir in range(0, len(list_fit_min)):
                            fit_min = list_fit_min[ir];
                            fit_max = list_fit_max[ir];
                            ana_pi0.set_fit_range(fit_min, fit_max);
                            integral_min = list_integral_min[0];
                            integral_max = list_integral_max[0];
                            ana_pi0.set_integral_range(integral_min, integral_max);
                            yield_min = list_yield_min[0];
                            yield_max = list_yield_max[0];
                            ana_pi0.set_yield_range(yield_min, yield_max);
                            outlist_fit_range = ana_pi0.analyze_ptspectrum();
                            outlist_fit_range.SetName("fit_{0:3.2f}_{1:3.2f}_GeVc2".format(fit_min, fit_max));
                            outlist_func.Add(outlist_fit_range);
                outfile.WriteTObject(outlist_ss);
                outlist_ss.Clear();
            del ana_pi0;
            
            filename_plot = ""
            plot_pi0 = PlotInvMass(meson, outname, "pi0eta-to-gammagamma");
            plot_pi0.set_arr_pt(arr_pt);
            for isys in range(0,nsys):
                ssname = config[type]['subsystems'][isys]['name']; #subsystem name
                plot_pi0.set_subsystem(ssname);
                print("plot subsystem", ssname);
                cutnames = config[type]["subsystems"][isys]['cutnames']
                print("cutnames", cutnames); 
                nc = len(cutnames);
                nfit = len(list_fit_func)
                list_parameters_linear_comparison = [];
                yield_list_lin = [];
                #list_parameters_quadratic_comparison = [];
                for ic in range(0,nc):
                    cutname = cutnames[ic];
                    print("cutname ", cutname)
                    plot_pi0.set_cutname(cutname);
                    for ifit in range(0, nfit): 
                        fitname = list_fit_func[ifit];
                        print("fitname ", fitname)
                        plot_pi0.set_fitname(fitname);
                        for ir in range(len(list_fit_min)):
                            fit_min = list_fit_min[ir];
                            fit_max = list_fit_max[ir];
                            plot_pi0.set_fit_range(fit_min, fit_max);
                            yield_min = list_yield_min[0];
                            yield_max = list_yield_max[0];
                            plot_pi0.set_yield_range(yield_min, yield_max);
                            plottingRange = [0., 0.3];
                            list_plot = plot_pi0.set_fit_list();
                            if ifit == 0:
                                yield_list_lin.append(list_plot.FindObject("h1yield_param"))
                            function_list = [];
                            histogram_list = [];
                            parameter_list = [];
                            same_list = [];
                            mixed_list = [];
                            it = ROOT.TIter(list_plot)
                            obj = it.Next()
                            while obj:
                                objName = obj.GetName()
                                if "param" in objName:
                                    parameter_list.append(obj)
                                # if "yield" in objName:
                                #     yield_list.append(obj)
                                obj = it.Next()

        # pdf output of all parameters for each combination of cut and fit                    
                            plot_histo = PlotHistoInvMass(meson, outname,"pi0eta-to-gammagamma")
                            outname_histo = os.path.join(folder, "{0}_{1}_InvMass_Fitparameters_{2}_{3}.pdf".format(date, period, cutname, fitname));
                            plot_histo.PlotHistoParameters(parameter_list, TString(outname_histo), "", "", period, 2,3,"#gamma#gamma", "data")

        # pdf output of all fitted histograms for each combination of cut and fit
                            for ipt in range(len(arr_pt)):
                                histo_ipt =list_plot.FindObject("h1mgg_subtracted_pt{0}".format(ipt));
                                histogram_list.append(histo_ipt)
                                function_ipt = list_plot.FindObject("f1total_pt{0}".format(ipt));
                                function_list.append(function_ipt)
                            output_name = os.path.join(folder, "{0}_{1}_InvMass_Overview_{2}_{3}.pdf".format(date, period, cutname, fitname))
                            plot_pi0.PlotInvMassInPtBins(histogram_list, function_list, TString(output_name), "", "", plottingRange, period, 3,3 , 0, len(arr_pt),len(arr_pt), "#gamma#gamma", "Data")
        
        # pdf output of all same and mixed scaled histograms for each cut
                            for ipt in range(len(arr_pt)):
                                histo_ipt_same =list_plot.FindObject("h1mgg_same_pt{0}".format(ipt));
                                same_list.append(histo_ipt_same)
                                histo_ipt_mixed = list_plot.FindObject("h1mgg_mix_scaled_pt{0}".format(ipt));
                                mixed_list.append(histo_ipt_mixed)
                            output_name = os.path.join(folder, "{0}_{1}_InvMass_Scaled_Mixed_{2}_{3}.pdf".format(date, period, cutname, fitname))
                            plot_pi0.PlotSameMixedInPtBins(same_list, mixed_list, TString(output_name), "", "", plottingRange, period, 3,3 , 0, len(arr_pt),len(arr_pt), "#gamma#gamma", "Data")

        # pdf output of mass, amplitude and width for each fit and comparison of all cuts
                            if ifit ==0:
                                parameter_comparison = [];
                                for i_parameter in range(4):
                                    parameter_comparison.append(parameter_list[i_parameter]);
                                list_parameters_linear_comparison.append(parameter_comparison);
                            elif ifit ==1:
                                parameter_comparison = [];
                                for i_parameter in range(4):
                                    parameter_comparison.append(parameter_list[i_parameter]);
                                #.append(parameter_comparison);  
                          
                print(yield_list_lin, "HERE", len(yield_list_lin)) 
                plot_raw_yield = PlotRawYieldInvMass(meson, outname, "pi0eta-to-gammagamma")     
                outname_raw_yield = os.path.join(folder, "{0}_{1}_InvMass_RawYield_{2}.pdf".format(date, period, fitname)); 
                plot_raw_yield.PlotHistoYield(yield_list_lin, TString(outname_raw_yield), "", "", period, 1,1, "#gamma#gamma", "data", cutnames)                   
                plot_parameters = PlotHistoParametersCombined(meson, outname,"pi0eta-to-gammagamma")
                outname_histo_linear = os.path.join(folder, "{0}_{1}_InvMass_Parameters_Combined_{2}.pdf".format(date, period, fitname));
                #outname_histo_quadratic = os.path.join(folder, "InvMass_Parameters_quadratic.pdf");
                plot_parameters.PlotHistoParametersCombined(list_parameters_linear_comparison, TString(outname_histo_linear), "", "", period, 1, 5, "#gamma#gamma", "data", cutnames)
                #plot_parameters.PlotHistoParametersCombined(list_parameters_quadratic_comparison, TString(outname_histo_quadratic), "", "", date, 1, 5, "#gamma#gamma", "data", cutnames)

            del plot_pi0;
            del plot_histo;
            del parameter_list;
            del histogram_list;
            del function_list


        outfile.Close();
    else:
        print("please check what to do in",config);

    rootfile.Close();  
#_________________________________________________________________________________________

########################################
#        Run over single Dataset       #
########################################

period = "LHC22f"
filename = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_124838_LHC22f_pass4.root"

# period = "LHC22o_small"
# filename = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_129486_LHC22o_pass4_small.root"

# period = "LHC22o_minBias"
# filename = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_135860_LHC22o_pass4_minBias_medium.root"

cutname = "qc"
suffix = "AnyTrack";
type = "data"
config_file = "/Users/alicamarieenderich/202312_invariant_mass/invariant_mass_code/configs/config_pp_13.6TeV_pi0 copy.yml"
with open(config_file, "r", encoding="utf-8") as config_yml:
    config = yaml.safe_load(config_yml)
date = "this_thesis" #datetime.date.today().strftime("%Y%m%d");
folder = "/Users/alicamarieenderich/{0}_{1}_invariant_mass_plots/".format(date, period);  
os.makedirs(folder, exist_ok=True);

run(filename,config,type,suffix, folder);

########################################
#        Run over all Datasets         #
########################################

# period_array = ["LHC22f", "LHC22o_small", "LHC22o_minBias"]
# filename_array = ["/Users/alicamarieenderich/AnalysisResults/AnalysisResults_124838_LHC22f_pass4.root", 
#             "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_129486_LHC22o_pass4_small.root", 
#             "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_135860_LHC22o_pass4_minBias_medium.root"]

# for i in range(len(period_array)):
#     period = period_array[i]
#     filename = filename_array[i]
#     cutname = "qc"
#     suffix = "AnyTrack";
#     type = "data"
#     config_file = "/Users/alicamarieenderich/invariant_mass-1/configs/config_pp_13.6TeV_pi0.yml"
#     with open(config_file, "r", encoding="utf-8") as config_yml:
#         config = yaml.safe_load(config_yml)
#     date = "this_thesis" #datetime.date.today().strftime("%Y%m%d");
#     folder = "/Users/alicamarieenderich/{0}_{1}_invariant_mass_plots/".format(date, period);  
#     os.makedirs(folder, exist_ok=True);

#     run(filename,config,type,suffix, folder);



# period = "LHC23d1kk"
# filename = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_124837_LHC23d1k.root"

# # period = "LHC22o_small"
# # filename = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_129486_LHC22o_pass4_small.root"

# # period = "LHC22o_minBias"
# # filename = "/Users/alicamarieenderich/AnalysisResults/AnalysisResults_135860_LHC22o_pass4_minBias_medium.root"

# cutname = "qc"
# suffix = "AnyTrack";
# type = "mc"
# config_file = "/Users/alicamarieenderich/invariant_mass-1/configs/config_pp_13.6TeV_pi0.yml"
# with open(config_file, "r", encoding="utf-8") as config_yml:
#     config = yaml.safe_load(config_yml)
# date = "this_thesis" #datetime.date.today().strftime("%Y%m%d");
# folder = "/Users/alicamarieenderich/{0}_{1}_invariant_mass_plots/".format(date, period);  
# os.makedirs(folder, exist_ok=True);

# run(filename,config,type,suffix, folder);