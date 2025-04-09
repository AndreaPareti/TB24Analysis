#############################################
### Analysis script for 2024 Test Beam  #####
### of the dual-readout prototype DRAGO #####
### A. Pareti, 12 Feb 2025              #####
### andrea.pareti01@ateneo.pv.it        #####
#############################################

import pandas as pd
import numpy as np
import ROOT 
import uproot
from scipy import optimize

calibfolder = "/home/storage/data_apareti/TB24/ElectronEnergyScan/"
#infolder = "/home/storage/data_apareti/TB24/PionScan_OldHVsaturated/"
infolder = "/home/storage/data_apareti/TB24/PionEnergyScan/"
#infolder = "/home/storage/data_apareti/TB24/PionScanCorrected/pion24/"
pedfolder = "/home/storage/data_apareti/TB24/PedestalRuns/"
treename = "Ftree"

#chi_value = 0.35
#chi_value = 0.25
chi_value = 0.44


# Mean Z barycenter for electromagnetic showers at 40 GeV
meanZbarS_ele = 227.718  # in mm
meanZbarC_ele = 228.643
# Mean Z barycenter for had showers at 40 GeV
meanZbarS_had = 590.164
meanZbarC_had = 575.515
#meanZbarS_had = 650.164;
#meanZbarC_had = 650.515;

# attenuation length - TB2023
# double att_length_S = 1916.; 
# double att_length_C = 3889.;
# attenuation length - 3.5m
att_length_S = 4000. 
att_length_C = 4000.

S_attenuation_correction = (ROOT.TMath.Exp(-(2500-meanZbarS_had)/att_length_S) ) / (ROOT.TMath.Exp(-(2500-meanZbarS_ele)/att_length_S) );
C_attenuation_correction = (ROOT.TMath.Exp(-(2500-meanZbarC_had)/att_length_C) ) / (ROOT.TMath.Exp(-(2500-meanZbarC_ele)/att_length_C) );




def GetPMTnoise():
    print("\n... Calculating PMT noise")

    PedestalRuns = ["0771", "0777", "0780", "0781", "0784", "0796"]

    PmtNoiseHistS = ROOT.TH1D("PmtNoiseHistS", "PmtNoiseHistS", 100, -5., 5.)
    PmtNoiseHistC = ROOT.TH1D("PmtNoiseHistC", "PmtNoiseHistC", 100, -5., 5.)

    for index, run in enumerate(PedestalRuns):
        filename = "physics_sps2024_run" + run + ".root"

        root_file = uproot.open(pedfolder+filename)
        tree = root_file[treename]
        data = tree.arrays(library="pd")    

        print(filename, data)
        for pmtS, pmtC in zip(data["totPMTSene"], data["totPMTCene"]):
            PmtNoiseHistS.Fill(pmtS)
            PmtNoiseHistC.Fill(pmtC)

    cNoise = ROOT.TCanvas("cNoise", "cNoise", 1400, 1200)
    PmtNoiseHistS.SetLineColor(ROOT.kRed); PmtNoiseHistS.SetLineWidth(2);  PmtNoiseHistS.Draw()      
    PmtNoiseHistC.SetLineColor(ROOT.kBlue); PmtNoiseHistC.SetLineWidth(2);  PmtNoiseHistC.Draw("same")

    PmtNoiseHistS.Fit("gaus"); PmtNoiseHistC.Fit("gaus")
    RmsNoiseS = PmtNoiseHistS.GetFunction("gaus").GetParameter(2)
    RmsNoiseC = PmtNoiseHistC.GetFunction("gaus").GetParameter(2)
    print("rms S: ", RmsNoiseS, "\trms C: ", RmsNoiseC)

    cNoise.Draw()
    cNoise.SaveAs("Noise.png")
    return(RmsNoiseS, RmsNoiseC)




def GetDF(run, Cut, filename, energy, containment): 
    print("Using file: ", filename, energy)
    root_file = uproot.open(infolder+filename)
    tree = root_file[treename]
    data = tree.arrays(cut=Cut, library="pd")

    # define Asymmetry variable 
    data["AsymS"] = (data["TS24"] - data["TS21"]) / (data["TS24"] + data["TS21"] )
    data["AsymC"] = (data["TC24"] - data["TC21"]) / (data["TC24"] + data["TC21"] )

    # define partial barycenter variable 
    data["BaryS"] = (data["TS00"]-28.3*data["TS11"]+28.3*data["TS15"])/(data["TS00"]+data["TS11"]+data["TS15"])
    data["BaryC"] = (data["TC00"]-28.3*data["TC11"]+28.3*data["TC15"])/(data["TC00"]+data["TC11"]+data["TC15"])

    data["pmtS_cont"] = data["totPMTSene"]/containment
    data["pmtC_cont"] = data["totPMTCene"]/containment

    data["pmtS_cont_att"] = data["totPMTSene"]/containment/S_attenuation_correction
    data["pmtC_cont_att"] = data["totPMTCene"]/containment/C_attenuation_correction

    data["totDRene"] = (data["totPMTSene"]-chi_value*data["totPMTCene"]) / (1-chi_value)
    data["totDRene_cont"] = data["totDRene"]/containment
    data["totDRene_cont_att"] = (data["totPMTSene"]/S_attenuation_correction-chi_value*data["totPMTCene"]/C_attenuation_correction) / (1-chi_value) / containment

    # input truth energy to dataframe
    data["TruthE"] = energy
    return data



#def GetAsymProfiles(data, energy):
def GetAsymProfilesMPV(data, energy, mpv_sci, mpv_cer, mpv_dr):
    # Profile normalised energy Vs Asymmetry
    DReneAsymSprof = ROOT.TProfile("DReneAsymSprof_{0}GeV".format(energy), r" DR Energy profile over Asymmetry(S) {0}GeV; TS24-TS21 / TS24 + TS21; E_DR/E".format(energy), 100, -1, 1)
    DReneAsymCprof = ROOT.TProfile("DReneAsymCprof_{0}GeV".format(energy), r"DR Energy profile over Asymmetry(C) {0}GeV; TC24-TC21 / TC24 + TC21; E_DR/E".format(energy), 100, -1, 1)

    # Profile normalised energy Vs Asymmetry
    SciAsymprof = ROOT.TProfile("sciprof_{0}GeV".format(energy), "totPMTSene profile over Asymmetry(S) {0}GeV; TS24-TS21 / TS24 + TS21; totPMTSene/E".format(energy), 100, -1, 1)
    CerAsymprof = ROOT.TProfile("cerprof_{0}GeV".format(energy), "totPMTCene profile over Asymmetry(C) {0}GeV; TC24-TC21 / TC24 + TC21; totPMTCene/E".format(energy), 100, -1, 1)

    for DRene, asymS,  asymC in zip(data["totDRene_cont"].values, data["AsymS"].values, data["AsymC"].values):
        # profile with truth energy
        #DReneAsymSprof.Fill(asymS, DRene/energy)
        #DReneAsymCprof.Fill(asymC, DRene/energy)
        # profile with channel most probable value
        DReneAsymSprof.Fill(asymS, DRene/mpv_dr)
        DReneAsymCprof.Fill(asymC, DRene/mpv_dr)


    for pmtS, pmtC, asymS, asymC, tdcS, tdcC  in zip(data["pmtS_cont"], data["pmtC_cont"], data["AsymS"], data["AsymC"], data["TDC_TS11"].values, data["TDC_TC11"].values):
        # profile with truth energy
        #SciAsymprof.Fill(asymS, pmtS/energy)
        #CerAsymprof.Fill(asymC, pmtC/energy)
        # profile with channel most probable value
        SciAsymprof.Fill(asymS, pmtS/mpv_sci)
        CerAsymprof.Fill(asymC, pmtC/mpv_cer)

    # Fit asymmetry with 5 degree polynomial
    DReneAsymSprof.Fit("pol5", "Q", "", -0.9, 0.9)
    DReneAsymCprof.Fit("pol5", "Q", "", -0.9, 0.9)
    SciAsymprof.Fit("pol5", "Q", "", -0.9, 0.9)
    CerAsymprof.Fit("pol5", "Q", "", -0.9, 0.9) 
    # Get fitted function
    fDReneAsymS = DReneAsymSprof.GetFunction("pol5")
    fDReneAsymC = DReneAsymCprof.GetFunction("pol5")
    fPMTSeneAsym = SciAsymprof.GetFunction("pol5")
    fPMTCeneAsym = CerAsymprof.GetFunction("pol5")


    ROOT.gStyle.SetOptStat(0)
    # Draw Profile with fit
    ctotPMTeneAsymProf = ROOT.TCanvas("ctotPMTeneAsymProf{0}".format(energy),"PMT signal over Asymmetry ({0} GeV)".format(energy), 1400, 1200)
    SciAsymprof.SetMarkerColor(ROOT.kRed); SciAsymprof.SetMarkerStyle(20); fPMTSeneAsym.SetLineColor(ROOT.kRed)
    CerAsymprof.SetMarkerColor(ROOT.kBlue); CerAsymprof.SetMarkerStyle(20); fPMTCeneAsym.SetLineColor(ROOT.kBlue)
    SciAsymprof.Draw()
    CerAsymprof.Draw("same")
    fPMTSeneAsym.Draw("same"); fPMTCeneAsym.Draw("same")
    leg = ROOT.TLegend(0.53, 0.75, 0.97, 0.92)
    leg.SetTextSize(0.018)   
    #leg.AddEntry(HistComb, r"(S-#chi C)/(1-#chi)")
    leg.AddEntry(SciAsymprof, "Profile over S PMTs")
    leg.AddEntry(CerAsymprof, "Profile over C PMTs")
    leg.Draw()     
    ctotPMTeneAsymProf.SaveAs("PMTsignalAsym_{0}GeV.png".format(energy))

    return data, fDReneAsymS, fDReneAsymC, DReneAsymSprof, DReneAsymCprof, fPMTSeneAsym, fPMTCeneAsym, SciAsymprof, CerAsymprof


def GetAsymProfiles(data, energy):
    # Profile normalised energy Vs Asymmetry
    DReneAsymSprof = ROOT.TProfile("DReneAsymSprof_{0}GeV".format(energy), r" DR Energy profile over Asymmetry(S) {0}GeV; TS24-TS21 / TS24 + TS21; E_DR/E".format(energy), 100, -1, 1)
    DReneAsymCprof = ROOT.TProfile("DReneAsymCprof_{0}GeV".format(energy), r"DR Energy profile over Asymmetry(C) {0}GeV; TC24-TC21 / TC24 + TC21; E_DR/E".format(energy), 100, -1, 1)

    # Profile normalised energy Vs Asymmetry
    SciAsymprof = ROOT.TProfile("sciprof_{0}GeV".format(energy), "totPMTSene profile over Asymmetry(S) {0}GeV; TS24-TS21 / TS24 + TS21; totPMTSene/E".format(energy), 100, -1, 1)
    CerAsymprof = ROOT.TProfile("cerprof_{0}GeV".format(energy), "totPMTCene profile over Asymmetry(C) {0}GeV; TC24-TC21 / TC24 + TC21; totPMTCene/E".format(energy), 100, -1, 1)

    for DRene, asymS,  asymC in zip(data["totDRene_cont"].values, data["AsymS"].values, data["AsymC"].values):
        # profile with truth energy
        DReneAsymSprof.Fill(asymS, DRene/energy)
        DReneAsymCprof.Fill(asymC, DRene/energy)
        # profile with channel most probable value
        #DReneAsymSprof.Fill(asymS, DRene/mpv_dr)
        #DReneAsymCprof.Fill(asymC, DRene/mpv_dr)


    for pmtS, pmtC, asymS, asymC, tdcS, tdcC  in zip(data["pmtS_cont"], data["pmtC_cont"], data["AsymS"], data["AsymC"], data["TDC_TS11"].values, data["TDC_TC11"].values):
        # profile with truth energy
        SciAsymprof.Fill(asymS, pmtS/energy)
        CerAsymprof.Fill(asymC, pmtC/energy)
        # profile with channel most probable value
        #SciAsymprof.Fill(asymS, pmtS/mpv_sci)
        #CerAsymprof.Fill(asymC, pmtC/mpv_cer)

    # Fit asymmetry with 5 degree polynomial
    DReneAsymSprof.Fit("pol5", "Q", "", -0.9, 0.9)
    DReneAsymCprof.Fit("pol5", "Q", "", -0.9, 0.9)
    SciAsymprof.Fit("pol5", "Q", "", -0.9, 0.9)
    CerAsymprof.Fit("pol5", "Q", "", -0.9, 0.9) 
    # Get fitted function
    fDReneAsymS = DReneAsymSprof.GetFunction("pol5")
    fDReneAsymC = DReneAsymCprof.GetFunction("pol5")
    fPMTSeneAsym = SciAsymprof.GetFunction("pol5")
    fPMTCeneAsym = CerAsymprof.GetFunction("pol5")


    ROOT.gStyle.SetOptStat(0)
    # Draw Profile with fit
    ctotPMTeneAsymProf = ROOT.TCanvas("ctotPMTeneAsymProf{0}".format(energy),"PMT signal over Asymmetry ({0} GeV)".format(energy), 1400, 1200)
    SciAsymprof.SetMarkerColor(ROOT.kRed); SciAsymprof.SetMarkerStyle(20); fPMTSeneAsym.SetLineColor(ROOT.kRed)
    CerAsymprof.SetMarkerColor(ROOT.kBlue); CerAsymprof.SetMarkerStyle(20); fPMTCeneAsym.SetLineColor(ROOT.kBlue)
    SciAsymprof.Draw()
    CerAsymprof.Draw("same")
    fPMTSeneAsym.Draw("same"); fPMTCeneAsym.Draw("same")
    leg = ROOT.TLegend(0.53, 0.75, 0.97, 0.92)
    leg.SetTextSize(0.018)   
    #leg.AddEntry(HistComb, r"(S-#chi C)/(1-#chi)")
    leg.AddEntry(SciAsymprof, "Profile over S PMTs")
    leg.AddEntry(CerAsymprof, "Profile over C PMTs")
    leg.Draw()     
    ctotPMTeneAsymProf.SaveAs("PMTsignalAsym_{0}GeV.png".format(energy))

    return data, fDReneAsymS, fDReneAsymC, DReneAsymSprof, DReneAsymCprof, fPMTSeneAsym, fPMTCeneAsym, SciAsymprof, CerAsymprof


def DrawColzPlot(outfile, ctitle, varX, varY, labelPlot, labelX, labelY, nbinX, xmin, xmax, nbinY, ymin, ymax):
    cCanva = ROOT.TCanvas(ctitle, "{0};{1};{2}".format(labelPlot, labelX, labelY), 1400, 1200)
    colHist = ROOT.TH2D(ctitle, "{0};{1};{2}".format(labelPlot, labelX, labelY), nbinX, xmin, xmax, nbinY, ymin, ymax)
    for xval, yval in zip(varX, varY):
        colHist.Fill(xval, yval)
    colHist.Draw("colz")
    cCanva.SaveAs(outfile)


def DrawProfPlot(outfile, ctitle, varX, varY, labelPlot, labelX, labelY, nbinX, xmin, xmax, markerStyle=20, markerColor=ROOT.kBlue):
    cCanva = ROOT.TCanvas(ctitle, "{0};{1};{2}".format(labelPlot, labelX, labelY), 1400, 1200)
    myProf = ROOT.TProfile(ctitle, "{0};{1};{2}".format(labelPlot, labelX, labelY), nbinX, xmin, xmax)
    for xval, yval in zip(varX, varY):
        myProf.Fill(xval, yval)
    myProf.SetMarkerStyle(markerStyle); myProf.SetMarkerColor(markerColor)
    myProf.Draw("colz")
    cCanva.SaveAs(outfile)


def GetTDCProfiles(data, energy):
    # profile energy over TS11 or TC11 TDC information
    DReneStdcProf = ROOT.TProfile("DReneStdcprof_{0}GeV".format(energy), r" DR Energy profile over TDC_TS11 {0}GeV; TDC_TS11 (adc); E_DR/E".format(energy), 100, 600, 720)
    DReneCtdcProf = ROOT.TProfile("DReneCtdcprof_{0}GeV".format(energy), r" DR Energy profile over TDC_TC11 {0}GeV; TDC_TC11 (adc); E_DR/E".format(energy), 100, 600, 720)
    ScieneStdcProf = ROOT.TProfile("ScieneStdcprof_{0}GeV".format(energy), r"S Energy profile over TDC_TS11 {0}GeV; TDC_TS11 (adc); E_S/E".format(energy), 100, 600, 720)
    CereneCtdcProf = ROOT.TProfile("CereneStdcprof_{0}GeV".format(energy), r"C Energy profile over TDC_TC11 {0}GeV; TDC_TC11 (adc); E_C/E".format(energy), 100, 600, 720)

    # Fill profiles
    for DRene, asymS,  asymC, tdcS, tdcC in zip(data["totDRene_cont"].values, data["AsymS"].values, data["AsymC"].values, data["TDC_TS11"].values, data["TDC_TC11"].values):
        DReneStdcProf.Fill(tdcS, DRene/energy)
        DReneCtdcProf.Fill(tdcC, DRene/energy)

    for pmtS, pmtC, asymS, asymC, tdcS, tdcC  in zip(data["pmtS_cont"], data["pmtC_cont"], data["AsymS"], data["AsymC"], data["TDC_TS11"].values, data["TDC_TC11"].values):
        ScieneStdcProf.Fill(tdcS, pmtS/energy)
        CereneCtdcProf.Fill(tdcC, pmtC/energy)

    UpLimit = [680, 680, 680, 680]

    # Fit TDC with a straight line
    f1 = DReneStdcProf.Fit("pol1", "SQ", "", 620, UpLimit[0])
    f2 = DReneCtdcProf.Fit("pol1", "SQ", "", 620, UpLimit[1])
    f3 = ScieneStdcProf.Fit("pol1", "SQ", "", 620, UpLimit[2])
    f4 = CereneCtdcProf.Fit("pol1", "SQ", "", 620, UpLimit[3])
    # Get fitted function
    fDReneStdc = DReneStdcProf.GetFunction("pol1")
    fDReneCtdc = DReneCtdcProf.GetFunction("pol1")
    fScieneStdc = ScieneStdcProf.GetFunction("pol1")
    fCereneCtdc = CereneCtdcProf.GetFunction("pol1")

    # Get fitted function
    fDReneStdc.SetRange(600, 750)
    fDReneCtdc.SetRange(600, 750)
    fScieneStdc.SetRange(600, 750)
    fCereneCtdc.SetRange(600, 750)

    ROOT.gStyle.SetOptStat(0)
    # Draw DR ene over S and C tdcs
    cDReneTDC = ROOT.TCanvas("cDReneTdc{0}".format(energy),"DR energy Over TDCs ({0} GeV)".format(energy), 1400, 1200)
    DReneStdcProf.SetMarkerColor(ROOT.kRed); DReneStdcProf.SetMarkerStyle(20); DReneStdcProf.SetLineColor(ROOT.kRed)
    DReneCtdcProf.SetMarkerColor(ROOT.kBlue); DReneCtdcProf.SetMarkerStyle(20); DReneCtdcProf.SetLineColor(ROOT.kBlue)
    fDReneStdc.SetLineColor(ROOT.kRed); fDReneCtdc.SetLineColor(ROOT.kBlue)
    DReneStdcProf.Draw(); DReneCtdcProf.Draw("same"); fDReneStdc.Draw("same"); fDReneCtdc.Draw("same")
    leg = ROOT.TLegend(0.53, 0.75, 0.97, 0.92)
    leg.SetTextSize(0.018)   
    #leg.AddEntry(HistComb, r"(S-#chi C)/(1-#chi)")
    leg.AddEntry(DReneStdcProf, "Profile over S PMTs")
    leg.AddEntry(DReneCtdcProf, "Profile over C PMTs")
    leg.Draw()     
    cDReneTDC.SaveAs("DRenergyOverTDCT11_{0}GeV.png".format(energy))
    ROOT.gStyle.SetOptStat(111)

    # Draw DR ene over S and C tdcs
    cPMTeneTDC = ROOT.TCanvas("cPMTeneTdc{0}".format(energy),"PMT energy Over TDCs ({0} GeV)".format(energy), 1400, 1200)
    ScieneStdcProf.SetMarkerColor(ROOT.kRed); ScieneStdcProf.SetMarkerStyle(20); ScieneStdcProf.SetLineColor(ROOT.kRed)
    CereneCtdcProf.SetMarkerColor(ROOT.kBlue); CereneCtdcProf.SetMarkerStyle(20); CereneCtdcProf.SetLineColor(ROOT.kBlue)
    fScieneStdc.SetLineColor(ROOT.kRed); fCereneCtdc.SetLineColor(ROOT.kBlue)
    ScieneStdcProf.Draw(); CereneCtdcProf.Draw("same"); fScieneStdc.Draw("same"); fCereneCtdc.Draw("same")
    leg = ROOT.TLegend(0.53, 0.75, 0.97, 0.92)
    leg.SetTextSize(0.018)   
    #leg.AddEntry(HistComb, r"(S-#chi C)/(1-#chi)")
    leg.AddEntry(ScieneStdcProf, "Profile over S PMTs")
    leg.AddEntry(CereneCtdcProf, "Profile over C PMTs")
    leg.Draw()     
    cPMTeneTDC.SaveAs("PMTenergyOverTDCT11_{0}GeV.png".format(energy))
    ROOT.gStyle.SetOptStat(111)


    # return profile and fitted function
    return DReneStdcProf, DReneCtdcProf, ScieneStdcProf, CereneCtdcProf, fDReneStdc, fDReneCtdc, fScieneStdc, fCereneCtdc


# Same as before, with considering most probable values instead of truth energy
def GetTDCProfilesMPV(data, energy, mpv_sci, mpv_cer, mpv_dr):

    # profile energy over TS11 or TC11 TDC information
    DReneStdcProf = ROOT.TProfile("DReneStdcprof_{0}GeV".format(energy), r" DR Energy profile over TDC_TS11 {0}GeV; TDC_TS11 (adc); E_DR/E".format(energy), 100, 600, 720)
    DReneCtdcProf = ROOT.TProfile("DReneCtdcprof_{0}GeV".format(energy), r" DR Energy profile over TDC_TC11 {0}GeV; TDC_TC11 (adc); E_DR/E".format(energy), 100, 600, 720)
    ScieneStdcProf = ROOT.TProfile("ScieneStdcprof_{0}GeV".format(energy), r"S Energy profile over TDC_TS11 {0}GeV; TDC_TS11 (adc); E_S/E".format(energy), 100, 600, 720)
    CereneCtdcProf = ROOT.TProfile("CereneStdcprof_{0}GeV".format(energy), r"C Energy profile over TDC_TC11 {0}GeV; TDC_TC11 (adc); E_C/E".format(energy), 100, 600, 720)

    # Fill profiles
    for DRene, asymS,  asymC, tdcS, tdcC in zip(data["totDRene_cont"].values, data["AsymS"].values, data["AsymC"].values, data["TDC_TS11"].values, data["TDC_TC11"].values):
        DReneStdcProf.Fill(tdcS, DRene/mpv_dr)
        DReneCtdcProf.Fill(tdcC, DRene/mpv_dr)

    for pmtS, pmtC, asymS, asymC, tdcS, tdcC  in zip(data["pmtS_cont"], data["pmtC_cont"], data["AsymS"], data["AsymC"], data["TDC_TS11"].values, data["TDC_TC11"].values):
        ScieneStdcProf.Fill(tdcS, pmtS/mpv_sci)
        CereneCtdcProf.Fill(tdcC, pmtC/mpv_cer)

    UpLimit = [680, 680, 680, 680]

    # Fit TDC with a straight line
    f1 = DReneStdcProf.Fit("pol1", "SQ", "", 620, UpLimit[0])
    f2 = DReneCtdcProf.Fit("pol1", "SQ", "", 620, UpLimit[1])
    f3 = ScieneStdcProf.Fit("pol1", "SQ", "", 620, UpLimit[2])
    f4 = CereneCtdcProf.Fit("pol1", "SQ", "", 620, UpLimit[3])
    # Get fitted function
    fDReneStdc = DReneStdcProf.GetFunction("pol1")
    fDReneCtdc = DReneCtdcProf.GetFunction("pol1")
    fScieneStdc = ScieneStdcProf.GetFunction("pol1")
    fCereneCtdc = CereneCtdcProf.GetFunction("pol1")

    # Get fitted function
    fDReneStdc.SetRange(600, 750)
    fDReneCtdc.SetRange(600, 750)
    fScieneStdc.SetRange(600, 750)
    fCereneCtdc.SetRange(600, 750)

    ROOT.gStyle.SetOptStat(0)
    # Draw DR ene over S and C tdcs
    cDReneTDC = ROOT.TCanvas("cDReneTdc{0}".format(energy),"DR energy Over TDCs ({0} GeV)".format(energy), 1400, 1200)
    DReneStdcProf.SetMarkerColor(ROOT.kRed); DReneStdcProf.SetMarkerStyle(20); DReneStdcProf.SetLineColor(ROOT.kRed)
    DReneCtdcProf.SetMarkerColor(ROOT.kBlue); DReneCtdcProf.SetMarkerStyle(20); DReneCtdcProf.SetLineColor(ROOT.kBlue)
    fDReneStdc.SetLineColor(ROOT.kRed); fDReneCtdc.SetLineColor(ROOT.kBlue)
    DReneStdcProf.Draw(); DReneCtdcProf.Draw("same"); fDReneStdc.Draw("same"); fDReneCtdc.Draw("same")
    leg = ROOT.TLegend(0.53, 0.75, 0.97, 0.92)
    leg.SetTextSize(0.018)   
    #leg.AddEntry(HistComb, r"(S-#chi C)/(1-#chi)")
    leg.AddEntry(DReneStdcProf, "Profile over S PMTs")
    leg.AddEntry(DReneCtdcProf, "Profile over C PMTs")
    leg.Draw()     
    cDReneTDC.SaveAs("DRenergyOverTDCT11_{0}GeV.png".format(energy))
    ROOT.gStyle.SetOptStat(111)

    # Draw DR ene over S and C tdcs
    cPMTeneTDC = ROOT.TCanvas("cPMTeneTdc{0}".format(energy),"PMT energy Over TDCs ({0} GeV)".format(energy), 1400, 1200)
    ScieneStdcProf.SetMarkerColor(ROOT.kRed); ScieneStdcProf.SetMarkerStyle(20); ScieneStdcProf.SetLineColor(ROOT.kRed)
    CereneCtdcProf.SetMarkerColor(ROOT.kBlue); CereneCtdcProf.SetMarkerStyle(20); CereneCtdcProf.SetLineColor(ROOT.kBlue)
    fScieneStdc.SetLineColor(ROOT.kRed); fCereneCtdc.SetLineColor(ROOT.kBlue)
    ScieneStdcProf.Draw(); CereneCtdcProf.Draw("same"); fScieneStdc.Draw("same"); fCereneCtdc.Draw("same")
    leg = ROOT.TLegend(0.53, 0.75, 0.97, 0.92)
    leg.SetTextSize(0.018)   
    #leg.AddEntry(HistComb, r"(S-#chi C)/(1-#chi)")
    leg.AddEntry(ScieneStdcProf, "Profile over S PMTs")
    leg.AddEntry(CereneCtdcProf, "Profile over C PMTs")
    leg.Draw()     
    cPMTeneTDC.SaveAs("PMTenergyOverTDCT11_{0}GeV.png".format(energy))
    ROOT.gStyle.SetOptStat(111)


    # return profile and fitted function
    return DReneStdcProf, DReneCtdcProf, ScieneStdcProf, CereneCtdcProf, fDReneStdc, fDReneCtdc, fScieneStdc, fCereneCtdc


def GetProfileFit(data, energy, funcname, varx, vary, nbinsX, xmin, xmax, profname, labelX, labelY, fitmin, fitmax):
    prof = ROOT.TProfile("{0}_{1}".format(profname, energy), "{0}_{1};{2};{3}".format(profname, energy, labelX, labelY), nbinsX, xmin, xmax)
    for x, y in zip(data[varx], data[vary]):
        prof.Fill(x, y/energy)
    prof.Fit(funcname, "Q", "", fitmin, fitmax)
    fit = prof.GetFunction(funcname)
    cTestProfile = ROOT.TCanvas("cProfile_{0}".format(profname), "{0}_{1}".format(profname, energy), 1400, 1200)
    prof.SetMarkerStyle(20); prof.SetMarkerColor(ROOT.kBlack); prof.SetMarkerSize(2); fit.SetLineColor(ROOT.kRed)
    prof.Draw()
    fit.Draw("same")
    cTestProfile.SaveAs("profile_{0}_{1}.png".format(profname, energy))
    return prof, fit




def DrawEnergyHist(dfs, energies, binning_S, binning_C, varname_S, varname_C):
    ROOT.gStyle.SetOptStat(0)
    histlistS = []; histlistC = []; histlistDR = []; histlistDRcorr = []
    for i, (df, energy) in enumerate(zip(dfs, energies)):
        histS = ROOT.TH1D("histS{0}".format(energy), "histS{0}; E [GeV]; Normalized Counts".format(energy), 160, 0., 160)
        histC = ROOT.TH1D("histC{0}".format(energy), "histC{0}; E [GeV]; Normalized Counts".format(energy), 160, 0., 160)
        histDR = ROOT.TH1D("histDR{0}".format(energy), "histDR{0}; E [GeV]; Normalized Counts".format(energy), 160, 0., 160)
        histDRcorr = ROOT.TH1D("histDRcorr{0}".format(energy), "histDRcorr{0}; E [GeV]; Normalized Counts".format(energy), 160, 0., 160)


        for varS, varC, varDR, varDRcorr in zip(df[varname_S], df[varname_C], df["totDRene"], df["totDRene_cont"]): 
            histS.Fill(varS)
            histC.Fill(varC)
            histDR.Fill(varDR)
            histDRcorr.Fill(varDRcorr)

        histS.SetLineColor(i+1)
        histC.SetLineColor(i+1)
        histDR.SetLineColor(i+1)
        histDRcorr.SetLineColor(i+1)

        histlistS.append(histS)
        histlistC.append(histC)
        histlistDR.append(histDR)
        histlistDRcorr.append(histDRcorr)        

    c_histS = ROOT.TCanvas("c_histS", "c_histS", 1400, 1200)
    c_histS.SetRightMargin(0.05)
    c_histS.SetLeftMargin(0.15)
    leg = ROOT.TLegend(0.67, 0.60, 0.92, 0.88)
    leg.SetTextSize(0.018)
    for hist, energy in zip(histlistS, energies):
        hist.SetTitle(varname_S)
        integral = hist.Integral()
        scale_factor = 1/integral
        hist.SetLineWidth(2)
        hist.Scale(scale_factor)
        hist.Draw("same hist")
        leg.AddEntry(hist,"{0}, {1} GeV".format(varname_S, energy))
        leg.Draw()
    c_histS.SaveAs("HistCollection_S_{0}.png".format(varname_S))

    c_histC = ROOT.TCanvas("c_histC", "c_histC", 1400, 1200)
    c_histC.SetRightMargin(0.05)
    c_histC.SetLeftMargin(0.15)
    leg = ROOT.TLegend(0.67, 0.60, 0.92, 0.88)
    leg.SetTextSize(0.018)   
    for hist, energy in zip(histlistC, energies):
        hist.SetTitle(varname_C)
        integral = hist.Integral()
        scale_factor = 1/integral
        hist.Scale(scale_factor)
        hist.SetLineWidth(2)
        hist.Draw("same hist")
        leg.AddEntry(hist, "{0}, {1} GeV".format(varname_C, energy))
        leg.Draw()
    c_histC.SaveAs("HistCollection_C_{0}.png".format(varname_C))
    
    c_histDR = ROOT.TCanvas("c_histDR", "c_histDR", 1400, 1200)
    c_histDR.SetRightMargin(0.05)
    c_histDR.SetLeftMargin(0.15)
    leg = ROOT.TLegend(0.67, 0.60, 0.92, 0.88)
    leg.SetTextSize(0.018)   
    for hist, energy in zip(histlistDR, energies):
        hist.SetTitle("Dual-Readout ene")
        integral = hist.Integral()
        scale_factor = 1/integral
        hist.Scale(scale_factor)
        hist.SetLineWidth(2)
        hist.Draw("same hist")
        leg.AddEntry(hist, "{0}, {1} GeV".format("Dual-Readout ene", energy))
        leg.Draw()
    c_histDR.SaveAs("HistCollection_DR_{0}.png".format("DRene"))
    
    c_histDRcorr = ROOT.TCanvas("c_histDRcorr", "c_histDRcorr", 1400, 1200)
    c_histDRcorr.SetRightMargin(0.05)
    c_histDRcorr.SetLeftMargin(0.15)
    leg = ROOT.TLegend(0.67, 0.60, 0.92, 0.88)
    leg.SetTextSize(0.018)   
    for hist, energy in zip(histlistDRcorr, energies):
        hist.SetTitle("DR ene, corrected for containment")
        integral = hist.Integral()
        scale_factor = 1/integral
        hist.Scale(scale_factor)
        hist.SetLineWidth(2)
        hist.Draw("same hist")
        leg.AddEntry(hist, "{0}, {1} GeV".format("DR ene, corrected for containment", energy))
        leg.Draw()
    c_histDRcorr.SaveAs("HistCollection_DRcorr_{0}.png".format("DReneCorr"))
    

def DrawFem(energy, eneS, eneC, labelS, labelC, title, outname):
    ROOT.gStyle.SetOptStat(0)
    cFem = ROOT.TCanvas(title, title, 1400, 1200)
    histfem = ROOT.TH2D("histfem{0}GeV".format(energy), "{3}, {0}GeV;{1}/E;{2}/E".format(energy, labelS, labelC, title), 50, 0., 1.5, 50, 0., 1.5)
    for s, c in zip(eneS, eneC):
        histfem.Fill(s/energy, c/energy)
    histfem.Draw("colz")
    cFem.SaveAs("{0}{1}.png".format(outname, energy))
    ROOT.gStyle.SetOptStat(1)











##########
## Main ##
##########
def main():
    print("Hello there")






    cutCalib = 5000


    # 20 GeV runs (new HV): 0766, 0999, 1018
    # open one run, get map of tower values
    

    runs = ["0786", "0766", "0772", "0774", "0775", "0778", "0779", "0792"]
    #runs = ["0786", "0999", "0772", "0774", "0775", "0778", "0779", "0792"]
    #runs = ["0786", "0766", "0772", "0774", "0775", "0776", "0779", "0792"]
    #runs = ["0786", "1018", "0772", "0774", "0775", "0776", "0779", "0792"] 



    run_energies = [10, 20, 30, 40, 60, 80, 100, 120]
    #myCut = "(abs(XDWC2 - XDWC1) < 5) & (abs(YDWC2 - YDWC1)<5) & (MCounter<200) & (TailC<300) & (C2>160) & ( (totLeakage - L20)<5000 ) & (PShower>550)"
    myCut = "(abs(XDWC2 - XDWC1) < 5) & (abs(YDWC2 - YDWC1)<5) & (MCounter<200) & (TailC<300) & (C2>160) & ( (totLeakage - L20)<5000 )"

    varProf = "YDWC2"
    cut_x_min = [-19.83, -16.74, -16.22, -15.95, -15.60, -16.12, -16.07, -15.50]
    cut_x_max = [23.90, 22.19, 23.27, 23.44, 24.27, 23.79, 23.63, 24.12]
    cut_y_min = [-26.54, -25.80, -26.15, -26.15, -26.39, -25.63, -25.63, -26.03]
    cut_y_max = [13.38, 10.89, 9.72, 9.50, 9.86, 10.89, 10.54, 10.17]


    #cut_y_min = [-20, -20, -20, -20, -20, -20, -20, -20]
    #cut_y_max = [4, 4, 4, 4, 4, 4, 4, 4]


    FitS = []; profS = []
    FitC = []; profC = [] 

    # Store ntuples as pandas dataframes

    #for index, (run, energy) in enumerate(zip(runs, energies)):

    run = "0999"
    calib_energy = 20
    filename = "physics_sps2024_run" + run + ".root"
    print(filename, calib_energy)
    ROOT.gStyle.SetOptFit(0)
    index=2
    print("Current cuts on DWCs: ", cut_x_min[index], cut_x_max[index], cut_y_min[index], cut_y_max[index])
    CurrentCut = myCut + " & (XDWC2 > {0}) & (XDWC2 < {1}) & (YDWC2 > {2}) & (YDWC2 < {3})".format(cut_x_min[index], cut_x_max[index], cut_y_min[index], cut_y_max[index]) 

    #df, funcS, funcC, eneSprof, eneCprof = GetDFparametrization(run, CurrentCut, filename, energy)
    root_file = uproot.open(calibfolder+filename)
    tree = root_file[treename]
    data = tree.arrays(cut=myCut, library="pd")
    data = data[:cutCalib]

    inputs_name_S = ['TS55', 'TS54', 'TS53', 'TS45', 'TS44', 'TS43',
                    'TS35', 'TS34', 'TS33', 'TS25', 'TS24', 'TS23', 'TS16',
                    'TS15', 'TS14', 'TS17', 'TS00', 'TS13', 'TS10', 'TS11',
                    'TS12', 'TS20', 'TS21', 'TS22', 'TS30', 'TS31', 'TS32',
                    'TS40', 'TS41', 'TS42', 'TS50', 'TS51', 'TS52', 'TS60',
                    'TS61', 'TS62']
    
    inputC_name_C = ['TC55', 'TC54', 'TC53', 'TC45', 'TC44', 'TC43',
                    'TC35', 'TC34', 'TC33', 'TC25', 'TC24', 'TC23', 'TC16',
                    'TC15', 'TC14', 'TC17', 'TC00', 'TC13', 'TC10', 'TC11',
                    'TC12', 'TC20', 'TC21', 'TC22', 'TC30', 'TC31', 'TC32',
                    'TC40', 'TC41', 'TC42', 'TC50', 'TC51', 'TC52', 'TC60',
                    'TC61', 'TC62']


    nBinsX = 3; nBinsY = 12
    mapSpmt = ROOT.TH2D("MapSpmt", "Map of S tower energies, run {0}".format(run), nBinsX, 0, nBinsX, nBinsY, 0, nBinsY)
    print("Tower ", "\tRow ", "\tcolumn")
    
    for index, tow in enumerate(inputs_name_S):
        column = index%3
        row = nBinsY-int(index/3)-1
        print(tow, "\t", row, "\t", column)

        for pmtS in data[tow]:
            #mapSpmt.Fill(index+1, pmtS/energy)
            mapSpmt.Fill(column+.5, row+.5, pmtS/calib_energy)
        #mapSpmt.GetXaxis().SetBinLabel(column+1,  )    
        

    # Overlay text inside each bin
    latex = ROOT.TLatex()
    latex.SetTextSize(0.04)  # Set text size
    latex.SetTextAlign(22)   # Center alignment
    cSpmt = ROOT.TCanvas(r"cSpmt", r"Map of S tower energies, run {0}; Tower; Pmt/E_beam".format(run))    
    mapSpmt.Draw("colz")
    #for i in range(1, 3+1):
    #    for j in range(1, 12+1):
    #        x = mapSpmt.GetXaxis().GetBinCenter(i)
    #        y = mapSpmt.GetYaxis().GetBinCenter(j)
    #        latex.DrawLatex(x, y, inputs_name[(i-1)*(j-1)+i-1])  # Print label in the bin
    for index, tow in enumerate(inputs_name_S):
        column = index%3
        row = nBinsY-int(index/3)-1
        print(tow, "\t", row, "\t", column)
        latex.DrawLatex(column+0.5, row+0.5, tow)  # Print label in the bin


    cSpmt.SaveAs("mapSpmt.png")

    A_S = np.array([
        data['TS55'].to_numpy(), data['TS54'].to_numpy(), data['TS53'].to_numpy(),
        data['TS45'].to_numpy(), data['TS44'].to_numpy(), data['TS43'].to_numpy(),
        data['TS35'].to_numpy(), data['TS34'].to_numpy(), data['TS33'].to_numpy(),
        data['TS25'].to_numpy(), data['TS24'].to_numpy(), data['TS23'].to_numpy(),
        data['TS16'].to_numpy(), data['TS15'].to_numpy(), data['TS14'].to_numpy(),
        data['TS17'].to_numpy(), data['TS00'].to_numpy(), data['TS13'].to_numpy(),
        data['TS10'].to_numpy(), data['TS11'].to_numpy(), data['TS12'].to_numpy(),
        data['TS20'].to_numpy(), data['TS21'].to_numpy(), data['TS22'].to_numpy(),
        data['TS30'].to_numpy(), data['TS31'].to_numpy(), data['TS32'].to_numpy(),
        data['TS40'].to_numpy(), data['TS41'].to_numpy(), data['TS42'].to_numpy(),
        data['TS50'].to_numpy(), data['TS51'].to_numpy(), data['TS52'].to_numpy(),
        data['TS60'].to_numpy(), data['TS61'].to_numpy(), data['TS62'].to_numpy()
        ])

    A_C = np.array([
        data['TC55'].to_numpy(), data['TC54'].to_numpy(), data['TC53'].to_numpy(),
        data['TC45'].to_numpy(), data['TC44'].to_numpy(), data['TC43'].to_numpy(),
        data['TC35'].to_numpy(), data['TC34'].to_numpy(), data['TC33'].to_numpy(),
        data['TC25'].to_numpy(), data['TC24'].to_numpy(), data['TC23'].to_numpy(),
        data['TC16'].to_numpy(), data['TC15'].to_numpy(), data['TC14'].to_numpy(),
        data['TC17'].to_numpy(), data['TC00'].to_numpy(), data['TC13'].to_numpy(),
        data['TC10'].to_numpy(), data['TC11'].to_numpy(), data['TC12'].to_numpy(),
        data['TC20'].to_numpy(), data['TC21'].to_numpy(), data['TC22'].to_numpy(),
        data['TC30'].to_numpy(), data['TC31'].to_numpy(), data['TC32'].to_numpy(),
        data['TC40'].to_numpy(), data['TC41'].to_numpy(), data['TC42'].to_numpy(),
        data['TC50'].to_numpy(), data['TC51'].to_numpy(), data['TC52'].to_numpy(),
        data['TC60'].to_numpy(), data['TC61'].to_numpy(), data['TC62'].to_numpy()
        ])

    print(A_S)
    print(A_S.T)
    print(data.shape)
    print(data.shape[0])
    target = np.full(data.shape[0], calib_energy)
    sol_Sci = optimize.lsq_linear(A_S.T, target)
    weights_Sci = sol_Sci.x

    sol_Cer = optimize.lsq_linear(A_C.T, target)
    weights_Cer = sol_Cer.x


    print(weights_Sci)
    print(weights_Cer)    


    pmtShist = ROOT.TH1D("pmtSraw", "pmtSraw run {0}".format(run), 36, 1, 36)
    pmtSweights = ROOT.TH1D("pmtSweights", "pmtSweights run {0}".format(run), 36, 1, 36)
    for index, tow in enumerate(inputs_name_S):
        print(index, tow, data[tow].mean(), weights_Sci[index])
        pmtShist.Fill(index+1, data[tow].mean())
        pmtShist.GetXaxis().SetBinLabel(index+1, tow)
        pmtSweights.SetBinContent(index+1, weights_Sci[index])
        pmtSweights.GetXaxis().SetBinLabel(index+1, tow)

    pmtChist = ROOT.TH1D("pmtCraw", "pmtCraw run {0}".format(run), 36, 1, 36)
    pmtCweights = ROOT.TH1D("pmtCweightC", "pmtCweightC run {0}".format(run), 36, 1, 36)
    for index, tow in enumerate(inputC_name_C):
        print(index, tow, data[tow].mean(), weights_Cer[index])
        pmtChist.Fill(index+1, data[tow].mean())
        pmtChist.GetXaxis().SetBinLabel(index+1, tow)
        pmtCweights.SetBinContent(index+1, weights_Cer[index])
        pmtCweights.GetXaxis().SetBinLabel(index+1, tow)



    cCanvaSci = ROOT.TCanvas("cSci", "cSci", 1600, 600)
    pmtSweights.SetFillColor(ROOT.kRed)
    pmtShist.Scale(1/pmtSweights.GetMaximum())
    #pmtSweights.SetMaximum( max(pmtSweights.GetMaximum(), pmtShist.GetMaximum())*1.15 )
    pmtSweights.SetMaximum( 7 )
    print("Max weight S: ", pmtSweights.GetMaximum())
    line = ROOT.TLine(pmtShist.GetXaxis().GetXmin(), 1.0, pmtShist.GetXaxis().GetXmax(), 1.0)
    line.SetLineColor(ROOT.kGreen)  # Set line color to red
    line.SetLineWidth(2); line.SetLineStyle(2)          # Set line thickness    
    #pmtSweights.SetMaximum( pmtSweights.GetMaximum()*1.5 )
    pmtSweights.Draw("hist")
    pmtShist.Draw("sameHIST")
    line.Draw("same")
    cCanvaSci.SaveAs("testS.png")

    cCanvaCer = ROOT.TCanvas("cCer", "cCer", 1600, 600)
    pmtCweights.SetFillColor(ROOT.kRed)
    pmtChist.Scale(1/pmtCweights.GetMaximum())
    #pmtCweights.SetMaximum( max(pmtCweights.GetMaximum(), pmtChist.GetMaximum())*1.15 )
    pmtCweights.SetMaximum( 7 )
    line = ROOT.TLine(pmtShist.GetXaxis().GetXmin(), 1.0, pmtShist.GetXaxis().GetXmax(), 1.0)
    line.SetLineColor(ROOT.kGreen)  # Set line color to red
    line.SetLineWidth(2); line.SetLineStyle(2)          # Set line thickness    
    print("Max weight C: ", pmtCweights.GetMaximum())
    #pmtCweights.SetMaximum( pmtCweights.GetMaximum()*1.5 )
    pmtCweights.Draw("hist")
    pmtChist.Draw("sameHIST")
    line.Draw("same")
    cCanvaCer.SaveAs("testC.png")


    # Apply and test recalibration
    correctedSciene = np.dot(A_S.T, weights_Sci)
    correctedCerene = np.dot(A_C.T, weights_Cer)

    print( correctedSciene.shape )
    #print(data["totPMTSene"].mean(), data["totPMTSene"].rms())
    #print( correctedene.mean(), correctedene.rms() )

    hist1Sci = ROOT.TH1D("totPMTSene", "totPMTSene", 100, 0., calib_energy*1.5)
    hist2Sci = ROOT.TH1D("eneScorrected", "eneScorrected", 100, 0., calib_energy*1.5)
    for pmtS, eneS in zip(data["totPMTSene"], correctedSciene):
        hist1Sci.Fill(pmtS)
        hist2Sci.Fill(eneS)

    print("totPMTSene: ", hist1Sci.GetMean(), " +- ", hist1Sci.GetRMS())
    print("After Calib: ", hist2Sci.GetMean(), " +- ", hist2Sci.GetRMS())

    hist1Cer = ROOT.TH1D("totPMTCene", "totPMTCene", 100, 0., calib_energy*1.5)
    hist2Cer = ROOT.TH1D("eneCcorrected", "eneCcorrected", 100, 0.,calib_energy*1.5)
    for pmtC, eneC in zip(data["totPMTCene"], correctedCerene):
        hist1Cer.Fill(pmtC)
        hist2Cer.Fill(eneC)

    print("totPMTCene: ", hist1Cer.GetMean(), " +- ", hist1Cer.GetRMS())
    print("After Calib: ", hist2Cer.GetMean(), " +- ", hist2Cer.GetRMS())



    c2Sci = ROOT.TCanvas("c2Sci", "c2Sci", 700, 600)
    hist1Sci.SetLineColor(ROOT.kBlack)
    hist2Sci.SetLineColor(ROOT.kRed)
    hist1Sci.SetMaximum( max(hist1Sci.GetMaximum(), hist2Sci.GetMaximum())*1.1 )
    hist1Sci.Draw()
    hist2Sci.Draw("same")
    c2Sci.SaveAs("CalibTestSci.png")

    c2Cer = ROOT.TCanvas("c2Cer", "c2Cer", 700, 600)
    hist1Cer.SetLineColor(ROOT.kBlack)
    hist2Cer.SetLineColor(ROOT.kRed)
    hist1Cer.SetMaximum( max(hist1Cer.GetMaximum(), hist2Cer.GetMaximum())*1.1 )
    hist1Cer.Draw()
    hist2Cer.Draw("same")
    c2Cer.SaveAs("CalibTestCer.png")





    ##################################
    ####### Start running on energy scan ###
    ######################################







    #runs = ["0714", "0715", "0716", "0717", "0718", "0721"]
    #runs = ["0968", "0967", "0966", "0965", "0963", "0962"]
    #runs = ["0972", "0967", "0966", "0965", "0963", "0962"]
    runs = ["1000", "0967", "0966", "0965", "0963", "0962"]

    energies = [20, 40, 60, 80, 100, 120]
    exp_containment = [0.865, 0.87, 0.875, 0.88, 0.885, 0.89]
    #exp_containment = [0.875, 0.875, 0.875, 0.875, 0.875, 0.875]
   
    #myCut = "(abs(XDWC2 - XDWC1) < 5) & (abs(YDWC2 - YDWC1)<5) & (totPMTSene>0) & (PShower<500) & (TailC<400) & (totLeakage<7000) & (MCounter<150)"
    #myCut = "(abs(XDWC2 - XDWC1) < 5) & (abs(YDWC2 - YDWC1)<5) & (totPMTSene>0) & (PShower<450) & (TailC<400) & (totLeakage<6500) & (totLeakage>4700) & (MCounter<150) & (PShower>350) & (YDWC2>-20) & (YDWC2<5) & (XDWC2>-20) & (XDWC2<20)"
    #myCut = "(abs(XDWC2 - XDWC1) < 5) & (abs(YDWC2 - YDWC1)<5) & (totPMTSene>0) & (PShower<450) & (TailC<400) & (totLeakage<6500) & (totLeakage>4500) & (MCounter<150) & (PShower>350) & (YDWC2>-20) & (YDWC2<5) & (XDWC2>-20) & (XDWC2<20)"
    
    # Final one
    #myCut = "(abs(XDWC2 - XDWC1) < 5) & (abs(YDWC2 - YDWC1)<5) & (totPMTSene>0) & (PShower<450) & (TailC<400) & (TailC>170) & (totLeakage<6500) & (MCounter<160) & (PShower>350) & (YDWC2>-20) & (YDWC2<5) & (XDWC2>-20) & (XDWC2<20)"
    
    # test
    #myCut = "(abs(XDWC2 - XDWC1) < 5) & (abs(YDWC2 - YDWC1)<5) & (totPMTSene>0) & (PShower<450) & (TailC<400) & (TailC>170) & (totLeakage<6500) & (MCounter<160) & (PShower>350) & (YDWC2>-20) & (YDWC2<5) & (XDWC2>-20) & (XDWC2<20) & (totLeakage>4500)"
    #myCut = "(TDC_TS00>600)"
    #myCut = "(abs(XDWC2 - XDWC1) < 5) & (abs(YDWC2 - YDWC1)<5) & (totPMTSene>0) & (PShower<450) & (TailC<400) & (TailC>170) & (totLeakage<5500) & (MCounter<160) & (PShower>350) & (YDWC2>-20) & (YDWC2<5) & (XDWC2>-20) & (XDWC2<20)"
    myCut = "(abs(XDWC2 - XDWC1) < 5) & (abs(YDWC2 - YDWC1)<5) & (totPMTSene>0) & (PShower<450) & (TailC<400) & (TailC>170) & (totLeakage<6500) & (MCounter<160) & (PShower>350) & (YDWC2>-20) & (YDWC2<5) & (XDWC2>-20) & (XDWC2<20) & (L02+L03+L04)<1100"



    varProf = "YDWC2"
    cut_x_min = [-19.83, -16.74, -16.22, -15.95, -15.60, -16.12, -16.07, -15.50]
    cut_x_max = [23.90, 22.19, 23.27, 23.44, 24.27, 23.79, 23.63, 24.12]
    cut_y_min = [-26.54, -25.80, -26.15, -26.15, -26.39, -25.63, -25.63, -26.03]
    cut_y_max = [13.38, 10.89, 9.72, 9.50, 9.86, 10.89, 10.54, 10.17]

    # Asymmetry profiles
    DReneAsymFitSvec = []; DReneProfSvec = []
    DReneAsymFitCvec = []; DReneProfCvec = [] 
    SciAsymFitvec = []; SciAsymProfvec = []
    CerAsymFitvec = []; CerAsymProfvec = [] 

    # TDC profiles
    DReneTdcFitSvec = []; DReneProfSvec = []
    DReneTdcFitCvec = []; DReneProfCvec = [] 
    SciTdcFitvec = []; SciTdcProfvec = []
    CerTdcFitvec = []; CerTdcProfvec = [] 

    # store profiles for energy corrected for asymmetry and profiled over TS11 TDCs
    DReneAsym_TdcFitVec = []; DReneAsym_TdcProfVec = []

    # Store ntuples as pandas dataframes
    dfs = []

    # array of dataframe with corrected energies
    dfCorrected_array = []

    # Arrays to store reco energy parameters
    MeanVec_S=[]; RmsVec_S=[]; MeanErrVec_S=[]; RmsErrVec_S=[] 
    MeanVec_C=[]; RmsVec_C=[]; MeanErrVec_C=[]; RmsErrVec_C=[] 
    MeanVec_Comb=[]; RmsVec_Comb=[]; MeanErrVec_Comb=[]; RmsErrVec_Comb=[] 
    MeanVec_Comb_Asym=[]; RmsVec_Comb_Asym=[]; MeanErrVec_Comb_Asym=[]; RmsErrVec_Comb_Asym=[] 

    # First loop: prepare data and extract parametrizations for asymmetry and TDC information    
    for index, (run, energy, cont) in enumerate(zip(runs, energies, exp_containment)):
        print("\n#######################\nSTARTING RUN {}\n#########################\n".format(run))
        filename = "physics_sps2024_run" + run + ".root"
        #print(filename, energy, cont)
        print("file: ", filename, "\tEnergy: ", energy, "\tExpected containment: ", cont)
        #ROOT.gStyle.SetOptFit(0)
        print("Current cuts on DWCs: ", cut_x_min[index], cut_x_max[index], cut_y_min[index], cut_y_max[index])
        #CurrentCut = myCut + " & (XDWC2 > {0}) & (XDWC2 < {1}) & (YDWC2 > {2}) & (YDWC2 < {3})".format(cut_x_min[index], cut_x_max[index], cut_y_min[index], cut_y_max[index]) 
        CurrentCut = myCut

        #df = GetDFparametrization(run, CurrentCut, filename, energy)
        # read data and add asymmetry variable

        df = GetDF(run, CurrentCut, filename, energy, cont)



        # Apply and test recalibration
        A_S = np.array([
            df['TS55'].to_numpy(), df['TS54'].to_numpy(), df['TS53'].to_numpy(),
            df['TS45'].to_numpy(), df['TS44'].to_numpy(), df['TS43'].to_numpy(),
            df['TS35'].to_numpy(), df['TS34'].to_numpy(), df['TS33'].to_numpy(),
            df['TS25'].to_numpy(), df['TS24'].to_numpy(), df['TS23'].to_numpy(),
            df['TS16'].to_numpy(), df['TS15'].to_numpy(), df['TS14'].to_numpy(),
            df['TS17'].to_numpy(), df['TS00'].to_numpy(), df['TS13'].to_numpy(),
            df['TS10'].to_numpy(), df['TS11'].to_numpy(), df['TS12'].to_numpy(),
            df['TS20'].to_numpy(), df['TS21'].to_numpy(), df['TS22'].to_numpy(),
            df['TS30'].to_numpy(), df['TS31'].to_numpy(), df['TS32'].to_numpy(),
            df['TS40'].to_numpy(), df['TS41'].to_numpy(), df['TS42'].to_numpy(),
            df['TS50'].to_numpy(), df['TS51'].to_numpy(), df['TS52'].to_numpy(),
            df['TS60'].to_numpy(), df['TS61'].to_numpy(), df['TS62'].to_numpy()
            ])

        A_C = np.array([
            df['TC55'].to_numpy(), df['TC54'].to_numpy(), df['TC53'].to_numpy(),
            df['TC45'].to_numpy(), df['TC44'].to_numpy(), df['TC43'].to_numpy(),
            df['TC35'].to_numpy(), df['TC34'].to_numpy(), df['TC33'].to_numpy(),
            df['TC25'].to_numpy(), df['TC24'].to_numpy(), df['TC23'].to_numpy(),
            df['TC16'].to_numpy(), df['TC15'].to_numpy(), df['TC14'].to_numpy(),
            df['TC17'].to_numpy(), df['TC00'].to_numpy(), df['TC13'].to_numpy(),
            df['TC10'].to_numpy(), df['TC11'].to_numpy(), df['TC12'].to_numpy(),
            df['TC20'].to_numpy(), df['TC21'].to_numpy(), df['TC22'].to_numpy(),
            df['TC30'].to_numpy(), df['TC31'].to_numpy(), df['TC32'].to_numpy(),
            df['TC40'].to_numpy(), df['TC41'].to_numpy(), df['TC42'].to_numpy(),
            df['TC50'].to_numpy(), df['TC51'].to_numpy(), df['TC52'].to_numpy(),
            df['TC60'].to_numpy(), df['TC61'].to_numpy(), df['TC62'].to_numpy()
            ])


        #df["totPMTSene"] = np.dot(A_S.T, weights_Sci)
        #df["totPMTCene"] = np.dot(A_C.T, weights_Cer)

        # Get the bin with the highest content (most probable bin)
        scienehist = ROOT.TH1D("scienehist{0}".format(energy), "totPMTSene after cuts, containment correction applied ({0} GeV); E [GeV]; Counts".format(energy), 100, energy-0.9*energy, energy+1.2*energy)
        cerenehist = ROOT.TH1D("cerenehist{0}".format(energy), "totPMTCene after cuts, containment correction applied ({0} GeV); E [GeV]; Counts".format(energy), 100, energy-0.9*energy, energy+1.2*energy)
        drenehist = ROOT.TH1D("drenehist{0}".format(energy), "DRene after cuts, containment correction applied ({0} GeV); E [GeV]; Counts".format(energy), 100, energy-0.9*energy, energy+1.2*energy)
        for eneS, eneC, eneDR in zip(df["pmtS_cont"], df["pmtC_cont"], df["totDRene_cont"]):
            scienehist.Fill(eneS)
            cerenehist.Fill(eneC)
            drenehist.Fill(eneDR)

        ROOT.gStyle.SetOptStat(111111)
        cS = ROOT.TCanvas("cS{0}".format(energy), "cS{0}".format(energy), 1400, 1200)
        cS.SetLeftMargin(0.14); cS.SetRightMargin(0.12)
        scienehist.SetLineColor(ROOT.kRed); scienehist.SetLineWidth(2)
        scienehist.Draw()
        cS.SaveAs("scieneTest{0}.png".format(energy))    

        cC = ROOT.TCanvas("cC{0}".format(energy), "cC{0}".format(energy), 1400, 1200)
        cC.SetLeftMargin(0.14); cC.SetRightMargin(0.12)
        cerenehist.SetLineColor(ROOT.kBlue); cerenehist.SetLineWidth(2)
        cerenehist.Draw()
        cC.SaveAs("cereneTest{0}.png".format(energy))    
        cDR = ROOT.TCanvas("cDR{0}".format(energy), "cDR{0}".format(energy), 1400, 1200)
        cDR.SetLeftMargin(0.14); cDR.SetRightMargin(0.12)
        drenehist.SetLineColor(ROOT.kGreen+1); drenehist.SetLineWidth(2)
        drenehist.Draw()
        cDR.SaveAs("dreneTest{0}.png".format(energy))    

        # Get the bin with the highest content (most probable bin)
        mp_bin_sci = scienehist.GetMaximumBin()
        mpv_sci = scienehist.GetXaxis().GetBinCenter(mp_bin_sci)
        mp_bin_cer = cerenehist.GetMaximumBin()
        mpv_cer = cerenehist.GetXaxis().GetBinCenter(mp_bin_cer)
        mp_bin_dr = drenehist.GetMaximumBin()
        mpv_dr = drenehist.GetXaxis().GetBinCenter(mp_bin_dr)


        print(f"Sciene MPV: {mpv_sci:.2f}")
        print(f"Cerene MPV: {mpv_cer:.2f}")
        print(f"DRene MPV: {mpv_dr:.2f}")


        #df, fDReneAsymS, fDReneAsymC, DReneAsymSprof, DReneAsymCprof,  fSciAsymProf, fCerAsymProf, SciAsymProf, CerAsymProf = GetDFparametrization(run, CurrentCut, filename, energy, cont)
        
        # simplify
        #df, fDReneAsymS, fDReneAsymC, DReneAsymSprof, DReneAsymCprof,  fSciAsymProf, fCerAsymProf, SciAsymProf, CerAsymProf = GetAsymProfiles(df, energy)
        # profile over most probable value instead of truth energy
        df, fDReneAsymS, fDReneAsymC, DReneAsymSprof, DReneAsymCprof,  fSciAsymProf, fCerAsymProf, SciAsymProf, CerAsymProf = GetAsymProfilesMPV(df, energy, mpv_sci, mpv_cer, mpv_dr)

        # clean events with TS11<100
        # cleanup for TDCs
        if(run=="1000"): 
            df["TDC_TS11"] = df["TDC_TS11"]-60; df["TDC_TC11"] = df["TDC_TC11"]-60
            df["TDC_TS15"] = df["TDC_TS15"]-60; df["TDC_TC15"] = df["TDC_TC15"]-60
            df["TDC_TS00"] = df["TDC_TS00"]-60; df["TDC_TC00"] = df["TDC_TC00"]-60


        min_tdc_sci = 600; max_tdc_sci = 720
        min_tdc_cer = 600; max_tdc_cer = 675

        df = df[df["TDC_TS00"]>600]
        df = df[df["TDC_TS11"]>600]
        df = df[df["TDC_TS15"]>600]
        df = df[df["TDC_TC11"]>600]
        df = df[df["TDC_TC15"]>600]
        df = df[df["TDC_TC00"]>600]
        df = df[df["TDC_TC11"]<675]
        df = df[df["TDC_TS11"]<720]


        # plot TDC profiles and return
        # profile using truth energy
        #DReneStdcProf, DReneCtdcProf, ScieneStdcProf, CereneCtdcProf, fDReneStdc, fDReneCtdc, fScieneStdc, fCereneCtdc = GetTDCProfiles(df, energy)
        # profile using most probable values of the three channels
        DReneStdcProf, DReneCtdcProf, ScieneStdcProf, CereneCtdcProf, fDReneStdc, fDReneCtdc, fScieneStdc, fCereneCtdc = GetTDCProfilesMPV(df, energy, mpv_sci, mpv_cer, mpv_dr)





        dfs.append(df)
        DrawFem(energy, df["totPMTSene"], df["totPMTCene"], "totPMTSene", "totPMTCene", "Electromagnetic Fraction, raw variables", "HistFemPMT")
        
        # append DRene profile and fits over S or C asymmetry
        DReneAsymFitSvec.append(fDReneAsymS); DReneAsymFitCvec.append(fDReneAsymC);  DReneProfSvec.append(DReneAsymSprof); DReneProfCvec.append(DReneAsymCprof)
        # append S/C energy profile and fits over S/C asymmetry
        SciAsymFitvec.append(fSciAsymProf); CerAsymFitvec.append(fCerAsymProf); SciAsymProfvec.append(SciAsymProf); CerAsymProfvec.append(CerAsymProf)
        # append DRene profile and fits over S and C tdcs
        DReneTdcFitSvec.append(fDReneStdc); DReneTdcFitCvec.append(fDReneCtdc);  DReneProfSvec.append(DReneStdcProf); DReneProfCvec.append(DReneCtdcProf)
        # append S/C energy over S/C tdcs
        SciTdcFitvec.append(fScieneStdc); CerTdcFitvec.append(fCereneCtdc); SciTdcProfvec.append(ScieneStdcProf); CerTdcProfvec.append(CereneCtdcProf)


        # Plot dual-readout energy profile over asymmetry variable
        #ROOT.gStyle.SetOptStat(0)
        ctest = ROOT.TCanvas("cDREneProfAsymmetry{0}".format(energy), "cDREneProfAsymmetry{0}".format(energy), 1400, 1200)
        ctest.SetLeftMargin(0.15)
        DReneAsymSprof.GetYaxis().SetTitle("PMT dual-readout energy / Beam nominal energy")
        DReneAsymSprof.SetLineWidth(2); DReneAsymSprof.SetLineColor(ROOT.kRed); DReneAsymSprof.SetMarkerStyle(ROOT.kFullCircle); DReneAsymSprof.SetMarkerColor(ROOT.kRed)
        DReneAsymSprof.SetMinimum(0.5); DReneAsymSprof.SetMaximum(1.5); DReneAsymSprof.Draw()
        DReneAsymCprof.SetLineWidth(2); DReneAsymCprof.SetLineColor(ROOT.kBlue); DReneAsymCprof.SetMarkerStyle(ROOT.kFullCircle); DReneAsymCprof.SetMarkerColor(ROOT.kBlue)
        DReneAsymCprof.SetMinimum(0.8); DReneAsymCprof.SetMaximum(1.2); DReneAsymCprof.Draw("same")
        leg = ROOT.TLegend(0.53, 0.75, 0.97, 0.92)
        leg.AddEntry(DReneAsymSprof, "DR ene profiled over S tower asymmetry", "L")
        leg.AddEntry(DReneAsymCprof, "DR ene profiled over C tower asymmetry", "L")
        leg.Draw()
        ctest.SaveAs("DREnergyProfOverAsymmetry_{0}GeV.png".format(energy))  

        # Plot S/C energy profile over asymmetry variable
        #ROOT.gStyle.SetOptStat(0)
        ctest = ROOT.TCanvas("cSciCerEneProfAsymmetry{0}".format(energy), "cSciCerEneProfAsymmetry{0}".format(energy), 1400, 1200)
        ctest.SetLeftMargin(0.15)
        SciAsymProf.GetYaxis().SetTitle("PMT dual-readout energy / Beam nominal energy")
        SciAsymProf.SetLineWidth(2); SciAsymProf.SetLineColor(ROOT.kRed); SciAsymProf.SetMarkerStyle(ROOT.kFullCircle); SciAsymProf.SetMarkerColor(ROOT.kRed)
        SciAsymProf.SetMinimum(0.5); SciAsymProf.SetMaximum(1.5); SciAsymProf.Draw()
        CerAsymProf.SetLineWidth(2); CerAsymProf.SetLineColor(ROOT.kBlue); CerAsymProf.SetMarkerStyle(ROOT.kFullCircle); CerAsymProf.SetMarkerColor(ROOT.kBlue)
        CerAsymProf.SetMinimum(0.8); CerAsymProf.SetMaximum(1.2); CerAsymProf.Draw("same")
        
        #leg = ROOT.TLegend(0.65, 0.15, 0.85, 0.25)
        leg = ROOT.TLegend(0.53, 0.75, 0.97, 0.92)
        leg.AddEntry(SciAsymProf, "Sciene profiled over S tower asymmetry", "L")
        leg.AddEntry(CerAsymProf, "Cerene profiled over C tower asymmetry", "L")
        leg.Draw()
        ctest.SaveAs("SciCerProfOverAsymmetry_{0}GeV.png".format(energy))  

        # Plot dual-readout energy profile over TDC_TS11 variable
        #ROOT.gStyle.SetOptStat(0)
        ctest = ROOT.TCanvas("cDREneProfAsymmetry{0}".format(energy), "cDREneProfAsymmetry{0}".format(energy), 1400, 1200)
        ctest.SetLeftMargin(0.15)
        DReneStdcProf.GetYaxis().SetTitle("PMT dual-readout energy / Beam nominal energy")
        DReneStdcProf.SetLineWidth(2); DReneStdcProf.SetLineColor(ROOT.kRed); DReneStdcProf.SetMarkerStyle(ROOT.kFullCircle); DReneStdcProf.SetMarkerColor(ROOT.kRed)
        DReneStdcProf.SetMinimum(0.5); DReneStdcProf.SetMaximum(1.5); DReneStdcProf.Draw()
        DReneCtdcProf.SetLineWidth(2); DReneCtdcProf.SetLineColor(ROOT.kBlue); DReneCtdcProf.SetMarkerStyle(ROOT.kFullCircle); DReneCtdcProf.SetMarkerColor(ROOT.kBlue)
        DReneCtdcProf.SetMinimum(0.8); DReneCtdcProf.SetMaximum(1.2); DReneCtdcProf.Draw("same")
        leg = ROOT.TLegend(0.53, 0.75, 0.97, 0.92)
        leg.AddEntry(DReneStdcProf, "DR ene profiled over S tower asymmetry", "L")
        leg.AddEntry(DReneCtdcProf, "DR ene profiled over C tower asymmetry", "L")
        leg.Draw()
        ctest.SaveAs("DREnergyProfOverTDC_{0}GeV.png".format(energy))  

        # Plot S/C energy profile over TDC variable
        #ROOT.gStyle.SetOptStat(0)
        ctest = ROOT.TCanvas("cSciCerEneProfTdc{0}".format(energy), "cSciCerEneProfTdc{0}".format(energy), 1400, 1200)
        ctest.SetLeftMargin(0.15)
        ScieneStdcProf.GetYaxis().SetTitle("PMT dual-readout energy / Beam nominal energy")
        ScieneStdcProf.SetLineWidth(2); ScieneStdcProf.SetLineColor(ROOT.kRed); ScieneStdcProf.SetMarkerStyle(ROOT.kFullCircle); ScieneStdcProf.SetMarkerColor(ROOT.kRed)
        ScieneStdcProf.SetMinimum(0.5); ScieneStdcProf.SetMaximum(1.5); ScieneStdcProf.Draw()
        CereneCtdcProf.SetLineWidth(2); CereneCtdcProf.SetLineColor(ROOT.kBlue); CereneCtdcProf.SetMarkerStyle(ROOT.kFullCircle); CereneCtdcProf.SetMarkerColor(ROOT.kBlue)
        CereneCtdcProf.SetMinimum(0.8); CereneCtdcProf.SetMaximum(1.2); CereneCtdcProf.Draw("same")
        
        #leg = ROOT.TLegend(0.65, 0.15, 0.85, 0.25)
        leg = ROOT.TLegend(0.53, 0.75, 0.97, 0.92)
        leg.AddEntry(ScieneStdcProf, "Sciene profiled over TDC_TS11", "L")
        leg.AddEntry(CereneCtdcProf, "Cerene profiled over TDC_TC11", "L")
        leg.Draw()
        ctest.SaveAs("SciCerProfOverTdc_{0}GeV.png".format(energy))  

        # END first loop with parametrization extraction 
        

    # comment since I'm only using 40 GeV parametrisation
    # extract asymmetry parametrisation for dual-readout energy at 40 GeV
    # with respect to S or C asymmetry, respectively
    fDRS40Asym  = np.vectorize(DReneAsymFitSvec[1].Eval); fDRC40Asym  = np.vectorize(DReneAsymFitCvec[1].Eval)

    # extract asymmetry parametrisation for Sci and Cer energy at 40 GeV
    fSci40Asym  = np.vectorize(SciAsymFitvec[1].Eval); fCer40Asym  = np.vectorize(CerAsymFitvec[1].Eval)

    # extract TDC parametrisation for dual-readout energy at 40 GeV
    fDReneStdc40  = np.vectorize(DReneTdcFitSvec[1].Eval); fDReneCtdc40  = np.vectorize(DReneTdcFitCvec[1].Eval)

    # extract asymmetry parametrisation for Sci and Cer energy at 40 GeV
    fSci40Tdc  = np.vectorize(SciTdcFitvec[1].Eval); fCer40Tdc  = np.vectorize(CerTdcFitvec[1].Eval)


    binning_S = [0, 30, 60, 90]
    binning_C = [0, 30, 60, 90]
    DrawEnergyHist(dfs, energies, binning_S, binning_C, "totPMTSene", "totPMTCene")

    
    # Read PMT noise
    noiseS, noiseC = GetPMTnoise()

    ROOT.gStyle.SetOptStat(1)

    # START LOOP 
    # Correct reco energy given parametrizations obtain with previous loop
        
    for index, (energy, cont) in enumerate(zip(energies, exp_containment)):
        print("Index ", index, "\tUsing file with energy: ", energy)

        data = dfs[index]

        # Define binnings on PMT energy
        conditions_DR = [data["totDRene_cont"] >= 0]
        choices_DR_asym = [data["totDRene"] / fDRS40Asym(data["AsymS"])/cont]
        choices_DR_tdc = [data["totDRene"] / fDReneStdc40(data["TDC_TS11"])/cont]


        # Use np.select to apply the conditions
        data["energyDR_AsymCorrected"] = np.select(conditions_DR, choices_DR_asym, default=data["totDRene_cont"])
        data["energyS_AsymCorrected"] = data["totPMTSene"]/fSci40Asym(data["AsymS"])/cont
        data["energyC_AsymCorrected"] = data["totPMTCene"]/fCer40Asym(data["AsymC"])/cont
        #data["energyC"] = np.select(conditions_C, choices_C, default=data["totPMTCene"])

        # Evaluate function on Asym variable and use it to correct energy
        DrawFem(energy, data["energyS_AsymCorrected"], data["energyC_AsymCorrected"], "S (corrected)", "C (corrected)", "Electromagnetic Fraction, Corrected for asymmetry", "HistFemCorrected")

        # Use np.select to apply the conditions
        data["energyDR_TdcCorrected"] = np.select(conditions_DR, choices_DR_tdc, default=data["totDRene_cont"])
        data["energyS_TdcCorrected"] = data["totPMTSene"]/fSci40Tdc(data["TDC_TS11"])/cont
        data["energyC_TdcCorrected"] = data["totPMTCene"]/fCer40Tdc(data["TDC_TC11"])/cont

        

        funcname = "pol1"; fitmin = 600; fitmax=750
        xmin = 620; xmax=750; nbinsX = 50
        varx = "TDC_TS11"; vary = "energyDR_AsymCorrected"
        profname = "DRene_AsymCorrectedOverTdc"; labelX = "TDC_TS11"; labelY = "DR reco E / Nominal Beam E"
        DReneAsym_TdcProfile, fDReneAsym_TdcFitProf = GetProfileFit(data, energy, funcname, varx, vary, nbinsX, xmin, xmax, profname, labelX, labelY, fitmin, fitmax)

        DReneAsym_TdcFitVec.append(fDReneAsym_TdcFitProf)
        DReneAsym_TdcProfVec.append(DReneAsym_TdcProfile)

        dfs[index] = data

    # END LOOP 


    # extract TDC parametrisation for dual-readout energy at 40 GeV corrected for asymmetry
    fDReneAsymTdc40  = np.vectorize(DReneAsym_TdcFitVec[1].Eval)

    # extract TDC parametrisation for dual-readout energy at 40 GeV corrected for asymmetry
    #fDReneAsymTdc40  = np.vectorize(DReneAsym_TdcFitVec[1].Eval)


    # START LOOP
    # Now All considered corrections have been added to dataframe
    # Plot before/after correction distributions

    for index, (energy, cont) in enumerate(zip(energies, exp_containment)):
        print("Index ", index, "\tUsing file with energy: ", energy)

        data = dfs[index]
        #data["energyDR_TdcAsym"] = np.select(conditions_DR_TdcAsym, choices_DR_TdcAsym, default=data["energyDR_AsymCorrected"])

        # Use TDC parametrisation over asymmetry-corrected DR energy
        conditions_DR_TdcAsym = [data["energyDR_AsymCorrected"]>=0]
        choices_DR_TdcAsym = [data["energyDR_AsymCorrected"] / fDReneStdc40(data["TDC_TS11"])]
        data["energyDR_TdcAsym"] = np.select(conditions_DR_TdcAsym, choices_DR_TdcAsym, default=data["energyDR_AsymCorrected"])


        # filter tree using only events with an asymmetry within range
        #data = data[ (np.abs(data["AsymS"]<AsymCut) ) & (np.abs(data["AsymC"]<AsymCut) ) ]
        #data_filtered = data[ (np.abs(data["AsymS"])<AsymCut ) & (np.abs(data["AsymC"])<AsymCut ) & (np.abs(data["BaryS"])<4) & (np.abs(data["BaryC"])<4) ]
        #data = data[ (np.abs(data["TDC_TS11"]-660)<TdcCut )]



        ########### Colz Plots  ###############
        # colz Sci energy plot (over asymmetry), before correction
        myOutfile = "ScieneColz_NoAsymCorrection_{0}GeV.png".format(energy); labelPlot= "S energy, before asymmetry correction {0} GeV".format(energy)
        varX = data["AsymS"]; varY = data["pmtS_cont"]; title = "cColzAsymS_pre"; labelX = "TS24-TS21 / TS24+TS21"; labelY = "totPMTSene/containment [GeV]"
        nbinX = 50; xmin = -1.4; xmax = 1.4; nbinY = 50; ymin = 0.; ymax = energy*1.7
        DrawColzPlot(myOutfile, title, varX, varY, labelPlot, labelX, labelY, nbinX, xmin, xmax, nbinY, ymin, ymax)

        # colz Sci energy plot, after correction with asymmetry
        myOutfile = "ScieneColz_AsymCorrection_{0}GeV.png".format(energy); labelPlot= "S energy, Asymmetry correction {0} GeV".format(energy)
        varX = data["AsymS"]; varY = data["energyS_AsymCorrected"]; title = "cColzAsymS_post"; labelX = "TS24-TS21 / TS24+TS21"; labelY = "totPMTSene/containment [GeV]"
        nbinX = 50; xmin = -1.4; xmax = 1.4; nbinY = 50; ymin = 0.; ymax = energy*1.7
        DrawColzPlot(myOutfile, title, varX, varY, labelPlot, labelX, labelY, nbinX, xmin, xmax, nbinY, ymin, ymax)

        # colz Cer energy plot (over asymmetry), before correction
        myOutfile = "CereneColz_NoAsymCorrection_{0}GeV.png".format(energy); labelPlot= "C energy, before asymmetry correction {0} GeV".format(energy)
        varX = data["AsymC"]; varY = data["pmtC_cont"]; title = "cColzAsymC_pre"; labelX = "TC24-TC21 / TC24+TC21"; labelY = "totPMTCene/containment [GeV]"
        nbinX = 50; xmin = -1.4; xmax = 1.4; nbinY = 50; ymin = 0.; ymax = energy*1.7
        DrawColzPlot(myOutfile, title, varX, varY, labelPlot, labelX, labelY, nbinX, xmin, xmax, nbinY, ymin, ymax)

        # colz Cer energy plot, after correction with asymmetry
        myOutfile = "CereneColz_AsymCorrection_{0}GeV.png".format(energy); labelPlot= "C energy, Asymmetry correction {0} GeV".format(energy)
        varX = data["AsymC"]; varY = data["energyC_AsymCorrected"]; title = "cColzAsymC_post"; labelX = "TC24-TC21 / TC24+TC21"; labelY = "totPMTCene/containment [GeV]"
        nbinX = 50; xmin = -1.4; xmax = 1.4; nbinY = 50; ymin = 0.; ymax = energy*1.7
        DrawColzPlot(myOutfile, title, varX, varY, labelPlot, labelX, labelY, nbinX, xmin, xmax, nbinY, ymin, ymax)


        # colz DR energy plot (over asymmetry), before correction
        myOutfile = "DReneColz_NoAsymCorrection_{0}GeV.png".format(energy); labelPlot= "DR energy, before asymmetry correction {0} GeV".format(energy)
        varX = data["AsymS"]; varY = data["totDRene_cont"]; title = "cColzAsymS"; labelX = "TS24-TS21 / TS24+TS21"; labelY = "Reco E_DR [GeV]"
        nbinX = 50; xmin = -1.4; xmax = 1.4; nbinY = 50; ymin = 0.; ymax = energy*1.7
        DrawColzPlot(myOutfile, title, varX, varY, labelPlot, labelX, labelY, nbinX, xmin, xmax, nbinY, ymin, ymax)

        # colz DR energy plot, after correction with asymmetry
        myOutfile = "DReneColz_AsymCorrection_{0}GeV.png".format(energy); labelPlot= "DR energy, Asymmetry correction {0} GeV".format(energy)
        varX = data["AsymS"]; varY = data["energyDR_AsymCorrected"]; title = "cColzAsymS"; labelX = "TS24-TS21 / TS24+TS21"; labelY = "Reco E_DR [GeV]"
        nbinX = 50; xmin = -1.4; xmax = 1.4; nbinY = 50; ymin = 0.; ymax = energy*1.7
        DrawColzPlot(myOutfile, title, varX, varY, labelPlot, labelX, labelY, nbinX, xmin, xmax, nbinY, ymin, ymax)


        # colz DR energy plot (over TDC 11), before any correction
        myOutfile = "DReneColz_NoTdcCorrection_{0}GeV.png".format(energy); labelPlot= "DR energy, before Tdc TS11 correction {0} GeV".format(energy)
        varX = data["TDC_TS11"]; varY = data["totDRene_cont"]; title = "cColzTdcTS11"; labelX = "TS11 Tdc"; labelY = "Reco E_DR [GeV]"
        nbinX = 50; xmin = 600; xmax = 720; nbinY = 50; ymin = 0.; ymax = energy*1.7
        DrawColzPlot(myOutfile, title, varX, varY, labelPlot, labelX, labelY, nbinX, xmin, xmax, nbinY, ymin, ymax)


        # colz DR energy plot, after correction with TDCs
        myOutfile = "DReneColz_TdcCorrection_{0}GeV.png".format(energy); labelPlot= "DR energy, Tdc TS11 correction {0} GeV".format(energy)
        varX = data["TDC_TS11"]; varY = data["energyDR_TdcCorrected"]; title = "cColzTdcTS11"; labelX = "TS11 Tdc"; labelY = "Reco E_DR [GeV]"
        nbinX = 50; xmin = 600; xmax = 720; nbinY = 50; ymin = 0.; ymax = energy*1.7
        DrawColzPlot(myOutfile, title, varX, varY, labelPlot, labelX, labelY, nbinX, xmin, xmax, nbinY, ymin, ymax)

        # colz DR energy plot, energy corrected with asymmetry, over TDCs
        myOutfile = "DReneColz_AsymCorrectedOverTdc_{0}GeV.png".format(energy); labelPlot= "DRene, asymmetry correction over Tdc TS11 {0} GeV".format(energy)
        varX = data["TDC_TS11"]; varY = data["energyDR_AsymCorrected"]; title = "cColzTdcTS11"; labelX = "TS11 Tdc"; labelY = "Reco E_DR [GeV]"
        nbinX = 50; xmin = 600; xmax = 720; nbinY = 50; ymin = 0.; ymax = energy*1.7
        DrawColzPlot(myOutfile, title, varX, varY, labelPlot, labelX, labelY, nbinX, xmin, xmax, nbinY, ymin, ymax)

        #######
        # colz DR energy plot, energy corrected with both asymmetry and TDCs, over TDCs
        myOutfile = "DReneColz_AsymTdcCorrectionOverTdc_{0}GeV.png".format(energy); labelPlot= "DRene, asymmetry&TDC correction over Tdc TS11 {0} GeV".format(energy)
        varX = data["TDC_TS11"]; varY = data["energyDR_TdcAsym"]; title = "cColzTdcTS11"; labelX = "TS11 Tdc"; labelY = "Reco E_DR [GeV]"
        nbinX = 50; xmin = 600; xmax = 720; nbinY = 50; ymin = 0.; ymax = energy*1.7
        DrawColzPlot(myOutfile, title, varX, varY, labelPlot, labelX, labelY, nbinX, xmin, xmax, nbinY, ymin, ymax)


        # Profile plot, energy corrected with asymmetry, over TDCs
        myOutfile = "DReneProf_AsymCorrectedOverTdc_{0}GeV.png".format(energy); labelPlot= "DRene, asymmetry correction over Tdc TS11 {0} GeV".format(energy)
        varX = data["TDC_TS11"]; varY = data["energyDR_AsymCorrected"]/energy; title = "cProfTdcTS11"; labelX = "TS11 Tdc"; labelY = "Reco E_DR / E_beam"
        nbinX = 50; xmin = 600; xmax = 720; nbinY = 50; ymin = energy*0.5; ymax = energy*1.5
        DrawProfPlot(myOutfile, title, varX, varY, labelPlot, labelX, labelY, nbinX, xmin, xmax)

        # Profile plot, energy corrected with asymmetry, over TDCs
        myOutfile = "DReneProf_AsymTdcCorrectedOverTdc_{0}GeV.png".format(energy); labelPlot= "DRene, asymmetry&TDC correction over Tdc TS11 {0} GeV".format(energy)
        varX = data["TDC_TS11"]; varY = data["energyDR_TdcAsym"]/energy; title = "cProfTdcTS11"; labelX = "TS11 Tdc"; labelY = "Reco E_DR / E_beam"
        nbinX = 50; xmin = 600; xmax = 720; nbinY = 50; ymin = energy*0.5; ymax = energy*1.5
        DrawProfPlot(myOutfile, title, varX, varY, labelPlot, labelX, labelY, nbinX, xmin, xmax)



        # colz Sci energy plot (over TDC_TS11), before correction
        myOutfile = "ScieneColz_NoTdcCorrection_{0}GeV.png".format(energy); labelPlot= "S energy, before timing correction {0} GeV".format(energy)
        varX = data["TDC_TS11"]; varY = data["pmtS_cont"]; title = "cColzScieneTdcTS11_pre"; labelX = "TDC_TS11"; labelY = "totPMTSene/containment [GeV]"
        nbinX = 50; xmin = min_tdc_sci; xmax = max_tdc_sci; nbinY = 50; ymin = 0.; ymax = energy*1.7
        DrawColzPlot(myOutfile, title, varX, varY, labelPlot, labelX, labelY, nbinX, xmin, xmax, nbinY, ymin, ymax)

        # colz Sci energy plot, after correction with tdc
        myOutfile = "ScieneColz_TdcCorrection_{0}GeV.png".format(energy); labelPlot= "S energy, after timing correction {0} GeV".format(energy)
        varX = data["TDC_TS11"]; varY = data["energyS_TdcCorrected"]; title = "cColzScieneTdcTS11_post"; labelX = "TDC_TS11"; labelY = "totPMTSene/containment [GeV]"
        nbinX = 50; xmin = min_tdc_sci; xmax = max_tdc_sci; nbinY = 50; ymin = 0.; ymax = energy*1.7
        DrawColzPlot(myOutfile, title, varX, varY, labelPlot, labelX, labelY, nbinX, xmin, xmax, nbinY, ymin, ymax)

        # colz Cer energy plot (over TDC_TC11), before correction
        myOutfile = "CereneColz_NoTdcCorrection_{0}GeV.png".format(energy); labelPlot= "C energy, before timing correction {0} GeV".format(energy)
        varX = data["TDC_TC11"]; varY = data["pmtC_cont"]; title = "cColzCereneTdcTC11_pre"; labelX = "TDC_TC11"; labelY = "totPMTCene/containment [GeV]"
        nbinX = 50; xmin = min_tdc_cer; xmax = max_tdc_cer; nbinY = 50; ymin = 0.; ymax = energy*1.7
        DrawColzPlot(myOutfile, title, varX, varY, labelPlot, labelX, labelY, nbinX, xmin, xmax, nbinY, ymin, ymax)

        # colz Cer energy plot, after correction with tdc
        myOutfile = "CereneColz_TdcCorrection_{0}GeV.png".format(energy); labelPlot= "C energy, after timing correction {0} GeV".format(energy)
        varX = data["TDC_TC11"]; varY = data["energyC_TdcCorrected"]; title = "cColzCereneTdcTC11_post"; labelX = "TDC_TC11"; labelY = "totPMTCene/containment [GeV]"
        nbinX = 50; xmin = min_tdc_cer; xmax = max_tdc_cer; nbinY = 50; ymin = 0.; ymax = energy*1.7
        DrawColzPlot(myOutfile, title, varX, varY, labelPlot, labelX, labelY, nbinX, xmin, xmax, nbinY, ymin, ymax)





        data["DRasym"] = (data["energyS_AsymCorrected"]-chi_value*data["energyC_AsymCorrected"])/(1-chi_value)
        data["DRtdc"] = (data["energyS_TdcCorrected"]-chi_value*data["energyC_TdcCorrected"])/(1-chi_value)


        # colz DR energy plot (over asymmetry), before correction
        myOutfile = "DReneColzAfterAsymCorrections_OverAsym_{0}GeV.png".format(energy); labelPlot= "DR energy built after S&C asymmetry corrections, over asymmetry {0} GeV".format(energy)
        varX = data["AsymS"]; varY = data["DRasym"]; title = "cColzDReneAsymS_OverAsymS"; labelX = "TS24-TS21 / TS24+TS21"; labelY = "Reco E_DR [GeV]"
        nbinX = 50; xmin = -1.4; xmax = 1.4; nbinY = 50; ymin = 0.; ymax = energy*1.7
        DrawColzPlot(myOutfile, title, varX, varY, labelPlot, labelX, labelY, nbinX, xmin, xmax, nbinY, ymin, ymax)
        myOutfile = "DReneProfAfterAsymCorrections_OverAsym_{0}GeV.png".format(energy); title = "cProfDReneAsymS_OverAsymS"
        DrawProfPlot(myOutfile, title, varX, varY, labelPlot, labelX, labelY, nbinX, xmin, xmax)


        # colz DR energy plot (over asymmetry), before correction
        myOutfile = "DReneColzAfterAsymCorrections_OverTdc_{0}GeV.png".format(energy); labelPlot= "DR energy built after S&C asymmetry corrections, over timing {0} GeV".format(energy)
        varX = data["TDC_TS11"]; varY = data["DRasym"]; title = "cColzDReneAsymSOverTdc"; labelX = "TDC_TS11"; labelY = "Reco E_DR [GeV]"
        nbinX = 50; xmin = min_tdc_sci; xmax = max_tdc_sci; nbinY = 50; ymin = 0.; ymax = energy*1.7
        DrawColzPlot(myOutfile, title, varX, varY, labelPlot, labelX, labelY, nbinX, xmin, xmax, nbinY, ymin, ymax)
        myOutfile = "DReneProfAfterAsymCorrections_OverTdc_{0}GeV.png".format(energy); title = "cProfDReneAsymSOverTdc"
        DrawProfPlot(myOutfile, title, varX, varY, labelPlot, labelX, labelY, nbinX, xmin, xmax)


        # colz DR energy plot, after correction with asymmetry
        myOutfile = "DReneColzAfterTdcCorrections_OverTdc_{0}GeV.png".format(energy); labelPlot= "DR energy built after S&C TDC corrections, over timing {0} GeV".format(energy)
        varX = data["TDC_TS11"]; varY = data["DRtdc"]; title = "cColzDReneTdc_OverTDC"; labelX = "TDC_TS11"; labelY = "Reco E_DR [GeV]"
        nbinX = 50; xmin = min_tdc_sci; xmax = max_tdc_sci; nbinY = 50; ymin = 0.; ymax = energy*1.7
        DrawColzPlot(myOutfile, title, varX, varY, labelPlot, labelX, labelY, nbinX, xmin, xmax, nbinY, ymin, ymax)
        myOutfile = "DReneProfAfterTdcCorrections_OverTdc_{0}GeV.png".format(energy); title = "cProfDReneTdc_OverTDC"
        DrawProfPlot(myOutfile, title, varX, varY, labelPlot, labelX, labelY, nbinX, xmin, xmax)


       # colz DR energy plot, after correction with asymmetry
        myOutfile = "DReneColzAfterTdcCorrections_OverAsym_{0}GeV.png".format(energy); labelPlot= "DR energy built after S&C TDC corrections, over asymmetry {0} GeV".format(energy)
        varX = data["AsymS"]; varY = data["DRtdc"]; title = "cColzDReneTdc_OverAsymS"; labelX = "TS24-TS21 / TS24+TS21"; labelY = "Reco E_DR [GeV]"
        nbinX = 50; xmin = -1.4; xmax = 1.4; nbinY = 50; ymin = 0.; ymax = energy*1.7
        DrawColzPlot(myOutfile, title, varX, varY, labelPlot, labelX, labelY, nbinX, xmin, xmax, nbinY, ymin, ymax)
        myOutfile = "DReneProfAfterTdcCorrections_OverAsym_{0}GeV.png".format(energy); itle = "cProfDReneTdc_OverAsymS"
        DrawProfPlot(myOutfile, title, varX, varY, labelPlot, labelX, labelY, nbinX, xmin, xmax)






        # filter tree using only events with an asymmetry within range
        AsymCut = 0.9
        #data = data[ (np.abs(data["AsymS"]<AsymCut) ) & (np.abs(data["AsymC"]<AsymCut) ) ]
        #data_filtered = data[ (np.abs(data["AsymS"])<AsymCut ) & (np.abs(data["AsymC"])<AsymCut ) & (np.abs(data["BaryS"])<4) & (np.abs(data["BaryC"])<4) ]
        data_filtered = data[ (np.abs(data["AsymS"])<AsymCut ) & (np.abs(data["AsymC"])<AsymCut ) ]


        # Fill histograms with corrected energy values
        HistScorrected = ROOT.TH1D("HistScorrected_{0}GeV".format(energy), "totPMTSene/avg_containment ({0}GeV); E [GeV]; Counts".format(energy), 100, energy-0.7*energy, energy+0.7*energy)
        #HistScorrected = ROOT.TH1D("HistSraw_{0}GeV".format(energy), "totPMTSene ({0}GeV); E [GeV]; Counts".format(energy), 100, energy-0.7*energy, energy+0.7*energy)
        HistScorrected_Asymcut = ROOT.TH1D("HistSAsymcorrected_{0}GeV".format(energy), "totPMTSene/avg_containment (Corrected with asymmetry, {0}GeV); E [GeV]; Counts".format(energy), 100, energy-0.7*energy, energy+0.7*energy)
        HistScorrected_Tdc = ROOT.TH1D("HistStdccorrected_{0}GeV".format(energy), "totPMTSene/avg_containment (Corrected with timing, {0}GeV); E [GeV]; Counts".format(energy), 100, energy-0.7*energy, energy+0.7*energy)


        # Fill histograms with corrected energy values
        HistCcorrected = ROOT.TH1D("HistCcorrected_{0}GeV".format(energy), "totPMTCene/avg_containment ({0}GeV); E [GeV]; Counts".format(energy), 100, energy-0.7*energy, energy+0.7*energy)
        #HistCcorrected = ROOT.TH1D("HistCraw_{0}GeV".format(energy), "totPMTCene ({0}GeV); E [GeV]; Counts".format(energy), 100, energy-0.7*energy, energy+0.7*energy)
        HistCcorrected_Asymcut = ROOT.TH1D("HistCAsymcorrected_{0}GeV".format(energy), "totPMTCene/avg_containment (Corrected with asymmetry, {0}GeV); E [GeV]; Counts".format(energy), 100, energy-0.7*energy, energy+0.7*energy)
        HistCcorrected_Tdc = ROOT.TH1D("HistCtdccorrected_{0}GeV".format(energy), "totPMTCene/avg_containment (Corrected with timing, {0}GeV); E [GeV]; Counts".format(energy), 100, energy-0.7*energy, energy+0.7*energy)

        # Fill histograms with corrected energy values
        HistCombCorrected = ROOT.TH1D("HistCombcorrected_{0}GeV".format(energy), r"(S-#chi C)/(1-#chi)/avg_containment  ({0}GeV); E [GeV]; Counts".format(energy), 100, energy-0.7*energy, energy+0.7*energy)
        HistComb = ROOT.TH1D("HistComb_{0}GeV".format(energy), r"(S-#chi C)/(1-#chi) ({0}GeV); E [GeV]; Counts".format(energy), 100, energy-0.7*energy, energy+0.7*energy)
        HistCombCorrected_Asymcut = ROOT.TH1D("HistDRene_AsymCorrected_{0}GeV".format(energy), r"(S-#chi C)/(1-#chi)/avg_containment (Corrected, {0}GeV); E [GeV]; Counts".format(energy), 100, energy-0.7*energy, energy+0.7*energy)

        # corrected energy for TDCs
        HistCombCorrected_tdc = ROOT.TH1D("HistDRene_TdcCorrected_{0}GeV".format(energy), r"(S-#chi C)/(1-#chi)/avg_containment (CorrectedTDC, {0}GeV); E [GeV]; Counts".format(energy), 100, energy-0.7*energy, energy+0.7*energy)


        TdcCut = 20
        # Filter only events with time close to mean 
        data = data[ (np.abs(data["TDC_TS11"]-666)<TdcCut )]
        data = data[ (np.abs(data["TDC_TC11"]-640)<TdcCut )]


 

        # Combining TDC and Asymmetry corrections
        HistDR = ROOT.TH1D("HistDRene_TdcAsym_{0}GeV".format(energy), r"(S-#chi C)/(1-#chi)/avg_containment (TDC&Asym corrections, {0}GeV); E [GeV]; Counts".format(energy), 100, energy-0.7*energy, energy+0.7*energy)

        HistDRAsym2 = ROOT.TH1D("HistDRAsym2_{0}GeV".format(energy), r"(S-#chi C)/(1-#chi)/avg_containment (Asym corrections, {0}GeV); E [GeV]; Counts".format(energy), 100, energy-0.7*energy, energy+0.7*energy)
        HistDRTdc2 = ROOT.TH1D("HistDRTdc2_{0}GeV".format(energy), r"reco E_DR (built with S and C corrected with timing, {0}GeV); E [GeV]; Counts".format(energy), 100, energy-0.7*energy, energy+0.7*energy)


        for pmtS, pmtC, eneDR_cont, eneDR_corrected in zip(data["pmtS_cont"].values, data["pmtC_cont"].values, data["totDRene_cont"].values, data["energyDR_TdcCorrected"]):
            HistScorrected.Fill(pmtS)
            HistCcorrected.Fill(pmtC)
            HistComb.Fill(eneDR_cont)
            HistCombCorrected.Fill(eneDR_cont)
            HistCombCorrected_tdc.Fill(eneDR_corrected)

        for eneS, eneC, eneDR_corrected in zip(data_filtered["energyS_AsymCorrected"], data_filtered["energyC_AsymCorrected"], data_filtered["energyDR_AsymCorrected"]):
            HistScorrected_Asymcut.Fill(eneS)
            HistCcorrected_Asymcut.Fill(eneC)
            HistCombCorrected_Asymcut.Fill(eneDR_corrected)
            HistDRAsym2.Fill( (eneS-chi_value*eneC)/(1-chi_value)   )


        for eneS, eneC in zip(data["energyS_TdcCorrected"], data["energyC_TdcCorrected"]):
            HistScorrected_Tdc.Fill(eneS)
            HistCcorrected_Tdc.Fill(eneC)
            HistDRTdc2.Fill((eneS-chi_value*eneC)/(1-chi_value) )

        # Show reco energy Vs Asymmetry
        #RecoEvsAsymSHist = ROOT.TH2D("RecoEvsAsymSHist_{0}GeV".format(energy), "RecoEvsAsymSHist_{0}GeV; TS24-TS21/TS24+TS21; Reco E_dr".format(energy), 50, -1.5, 1.5, 50, 0., 2*energy)
        #RecoEvsAsymCHist = ROOT.TH2D("RecoEvsAsymCHist_{0}GeV".format(energy), "RecoEvsAsymCHist_{0}GeV; TSC4-TC21/TC24+TC21; Reco E_dr".format(energy), 50, -1.5, 1.5, 50, 0., 2*energy)

        #for recoE, ts24, ts21, tc24, tc21 in zip(data["totDRene_cont"].values, data["TS24"].values, data["TS21"].values, data["TC24"].values, data["TC21"].values):
        #    RecoEvsAsymSHist.Fill(((ts24-ts21)/(ts24+ts21)), recoE)
        #    RecoEvsAsymCHist.Fill(((ts24-ts21)/(tc24+tc21)), recoE)
        #cAsymS = ROOT.TCanvas("cAsymS{}".format(energy), "cAsymS{}".format(energy), 1400, 1200)
        #RecoEvsAsymSHist.Draw("colz")
        #cAsymS.SaveAs("RecoEvsAysmS{0}GeV.png".format(energy))    
        #cAsymC = ROOT.TCanvas("cAsymC{}".format(energy), "cAsymC{}".format(energy), 1400, 1200)
        #RecoEvsAsymCHist.Draw("colz")
        #cAsymC.SaveAs("RecoEvsAysmC{0}GeV.png".format(energy))    



        # Filter only events with time close to mean
        #data = data[ (np.abs(data["TDC_TS11"]-666)<10 )]
        for ene in data["energyDR_TdcAsym"]:
            HistDR.Fill(ene)    


            
        # Fit S energy before and after asymmetry correction
        HistScorrected.Fit("gaus", "Q")
        ftotPMTSeneFit = HistScorrected.GetFunction("gaus") 
        HistScorrected_Asymcut.Fit("gaus", "Q")
        fSeneAsymFit = HistScorrected_Asymcut.GetFunction("gaus")
        # Fit C energy before and after asymmetry correction
        HistCcorrected.Fit("gaus", "Q")
        ftotPMTCeneFit = HistCcorrected.GetFunction("gaus")
        HistCcorrected_Asymcut.Fit("gaus", "Q")
        fCeneAsymFit = HistCcorrected_Asymcut.GetFunction("gaus")

        HistComb.Fit("gaus","Q")
        fDReneFit = HistComb.GetFunction("gaus")
        HistCombCorrected.Fit("gaus","Q")
        fDReneFitCorrected = HistCombCorrected.GetFunction("gaus")

        HistCombCorrected_Asymcut.Fit("gaus","Q")
        fDReneFitCorrected_Asymcut = HistCombCorrected_Asymcut.GetFunction("gaus")

        HistCombCorrected_tdc.Fit("gaus","Q")
        fDReneFitCorrected_tdc = HistCombCorrected_tdc.GetFunction("gaus")

        HistScorrected_Tdc.Fit("gaus", "Q")
        fScieneFitCorrected_tdc = HistScorrected_Tdc.GetFunction("gaus")
        
        HistCcorrected_Tdc.Fit("gaus", "Q")
        fCereneFitCorrected_tdc = HistCcorrected_Tdc.GetFunction("gaus")

        # Fit gaussian over both corrections
        HistDR.Fit("gaus", "Q")
        fHistDRFit = HistDR.GetFunction("gaus")

        # Fit gaussians after both corrections
        HistDRAsym2.Fit("gaus", "Q")
        fHistDRAsym2Fit = HistDRAsym2.GetFunction("gaus")

        HistDRTdc2.Fit("gaus", "Q")
        fHistDRTdc2Fit = HistDRTdc2.GetFunction("gaus")



        #ROOT.gStyle.SetOptFit(111)
        ROOT.gStyle.SetOptFit(0)
        ROOT.gStyle.SetOptStat(0)

        # refit within -1.5 and +3 sigma
        HistScorrected.Fit("gaus", "Q", "", ftotPMTSeneFit.GetParameter(1)-1.5*ftotPMTSeneFit.GetParameter(2), ftotPMTSeneFit.GetParameter(1)+1.5*ftotPMTSeneFit.GetParameter(2)) 
        HistCcorrected.Fit("gaus", "Q", "", ftotPMTCeneFit.GetParameter(1)-1.5*ftotPMTCeneFit.GetParameter(2), ftotPMTCeneFit.GetParameter(1)+1.5*ftotPMTCeneFit.GetParameter(2)) 
        HistComb.Fit("gaus", "Q", "", fDReneFit.GetParameter(1)-1.5*fDReneFit.GetParameter(2), fDReneFit.GetParameter(1)+1.5*fDReneFit.GetParameter(2)) 
        HistCombCorrected.Fit("gaus", "Q", "", fDReneFitCorrected.GetParameter(1)-1.5*fDReneFitCorrected.GetParameter(2), fDReneFitCorrected.GetParameter(1)+1.5*fDReneFitCorrected.GetParameter(2)) 
        HistCombCorrected_Asymcut.Fit("gaus", "Q", "", fDReneFitCorrected_Asymcut.GetParameter(1)-1.5*fDReneFitCorrected_Asymcut.GetParameter(2), fDReneFitCorrected_Asymcut.GetParameter(1)+1.5*fDReneFitCorrected_Asymcut.GetParameter(2)) 
        HistScorrected_Asymcut.Fit("gaus", "Q", "", fSeneAsymFit.GetParameter(1)-1.5*fSeneAsymFit.GetParameter(2), fSeneAsymFit.GetParameter(1)+1.5*fSeneAsymFit.GetParameter(2)) 
        HistCcorrected_Asymcut.Fit("gaus", "Q", "", fCeneAsymFit.GetParameter(1)-1.5*fCeneAsymFit.GetParameter(2), fCeneAsymFit.GetParameter(1)+1.5*fCeneAsymFit.GetParameter(2)) 

        HistScorrected_Tdc.Fit("gaus", "Q", "", fScieneFitCorrected_tdc.GetParameter(1)-1.5*fScieneFitCorrected_tdc.GetParameter(2), fScieneFitCorrected_tdc.GetParameter(1)+1.5*fScieneFitCorrected_tdc.GetParameter(2))
        HistCcorrected_Tdc.Fit("gaus", "Q", "", fCereneFitCorrected_tdc.GetParameter(1)-1.5*fCereneFitCorrected_tdc.GetParameter(2), fCereneFitCorrected_tdc.GetParameter(1)+1.5*fCereneFitCorrected_tdc.GetParameter(2))

        HistCombCorrected_tdc.Fit("gaus", "Q", "", fDReneFitCorrected_tdc.GetParameter(1)-1.5*fDReneFitCorrected_tdc.GetParameter(2), fDReneFitCorrected_tdc.GetParameter(1)+1.5*fDReneFitCorrected_tdc.GetParameter(2)) 

        # Fit histogram with both corrections
        HistDR.Fit("gaus", "Q", "", fHistDRFit.GetParameter(1)-1.5*fHistDRFit.GetParameter(2), fHistDRFit.GetParameter(1)+1.5*fHistDRFit.GetParameter(2))
        HistDRAsym2.Fit("gaus", "Q", "", fHistDRAsym2Fit.GetParameter(1)-1.5*fHistDRAsym2Fit.GetParameter(2), fHistDRAsym2Fit.GetParameter(1)+1.5*fHistDRAsym2Fit.GetParameter(2))
        HistDRTdc2.Fit("gaus", "Q", "", fHistDRTdc2Fit.GetParameter(1)-1.5*fHistDRTdc2Fit.GetParameter(2), fHistDRTdc2Fit.GetParameter(1)+1.5*fHistDRTdc2Fit.GetParameter(2))



        # Plot Sci Energy
        scihist = ROOT.TCanvas("cSciEnergy{0}".format(energy), "cSciEnergy{0}".format(energy), 1400, 1200)
        HistScorrected.SetLineColor(ROOT.kBlack); HistScorrected_Asymcut.SetLineColor(ROOT.kRed); HistScorrected.SetLineWidth(2); HistScorrected_Asymcut.SetLineWidth(2)
        HistScorrected.GetFunction("gaus").SetLineColor(ROOT.kBlack)
        HistScorrected_Asymcut.GetFunction("gaus").SetLineColor(ROOT.kRed); HistScorrected_Asymcut.GetFunction("gaus").SetLineStyle(2); HistScorrected.GetFunction("gaus").SetLineColor(ROOT.kBlack)
        HistScorrected_Tdc.GetFunction("gaus").SetLineColor(ROOT.kRed+1); HistScorrected_Tdc.GetFunction("gaus").SetLineWidth(2); HistScorrected.GetFunction("gaus").SetLineColor(ROOT.kBlack)
        HistScorrected_Tdc.SetLineWidth(2); HistScorrected_Tdc.SetLineColor(ROOT.kRed)

        maximum = max(HistScorrected.GetMaximum(), HistScorrected_Asymcut.GetMaximum(), HistScorrected_Tdc.GetMaximum())
        HistScorrected.SetMaximum(maximum*1.2)
        HistScorrected.Draw()
        #HistScorrected_Asymcut.Draw("same")
        HistScorrected_Tdc.Draw("same")
        str_rawS = r"pmtS/containment (raw): #mu = " + str( round(HistScorrected.GetFunction("gaus").GetParameter(1), 3) ) + r" #sigma = " + str( round(HistScorrected.GetFunction("gaus").GetParameter(2), 3) ) 
        str_asymS = r"S energy (corrected with asym): #mu = " + str( round(HistScorrected_Asymcut.GetFunction("gaus").GetParameter(1), 3) ) + r" #sigma = " + str( round(HistScorrected_Asymcut.GetFunction("gaus").GetParameter(2), 3) ) 
        str_tdcS = r"S energy (corrected with timing): #mu = " + str( round(HistScorrected_Tdc.GetFunction("gaus").GetParameter(1), 3) ) + r" #sigma = " + str( round(HistScorrected_Tdc.GetFunction("gaus").GetParameter(2), 3) )         
        #legend = ROOT.TLegend(0.15, 0.8, 0.4, 0.9); legend.AddEntry(HistScorrected, "SPMT/containment", "l"); legend.AddEntry(HistScorrected_Asymcut, "Asym Correction", "l"); 
        legend = ROOT.TLegend(0.15, 0.75, 0.6, 0.9); legend.AddEntry(HistScorrected, str_rawS, "l")
        #legend.AddEntry(HistScorrected_Asymcut, str_asymS, "l") 
        legend.AddEntry(HistScorrected_Tdc, str_tdcS, "l")
        # Make legend background transparent
        legend.SetFillStyle(0)  # 0 means fully transparent
        legend.SetFillColor(0)   # Ensure no background color
        legend.SetBorderSize(0)   # Thin border around the legend
        legend.SetTextSize(0.03); 
        legend.Draw("same")
        scihist.SaveAs("Scienehist{0}.png".format(energy))
        
        # Plot Cer Energy
        cerhist = ROOT.TCanvas("cCerEnergy{0}".format(energy), "cCerEnergy{0}".format(energy), 1400, 1200)
        HistCcorrected.SetLineColor(ROOT.kBlack); HistCcorrected_Asymcut.SetLineColor(ROOT.kBlue+1); HistCcorrected.SetLineWidth(2); HistCcorrected_Asymcut.SetLineWidth(2)
        HistCcorrected_Asymcut.GetFunction("gaus").SetLineColor(ROOT.kBlue); HistCcorrected_Asymcut.GetFunction("gaus").SetLineStyle(2); HistCcorrected.GetFunction("gaus").SetLineColor(ROOT.kBlack)
        maximum = max(HistCcorrected.GetMaximum(), HistCcorrected_Asymcut.GetMaximum())
        HistCcorrected.GetFunction("gaus").SetLineColor(ROOT.kBlack); HistCcorrected.GetFunction("gaus").SetLineWidth(2)
        HistCcorrected_Tdc.SetLineWidth(2); HistCcorrected_Tdc.SetLineColor(ROOT.kAzure+1); HistCcorrected_Tdc.GetFunction("gaus").SetLineWidth(2); HistCcorrected_Tdc.GetFunction("gaus").SetLineColor(ROOT.kBlue)

        HistCcorrected.SetMaximum(maximum*1.2)
        HistCcorrected.Draw()
        #HistCcorrected_Asymcut.Draw("same")
        HistCcorrected_Tdc.Draw("same")
        str_rawC = r"pmtC/containment (raw): #mu = " + str( round(HistCcorrected.GetFunction("gaus").GetParameter(1), 3) ) + r" #sigma = " + str( round(HistCcorrected.GetFunction("gaus").GetParameter(2), 3) ) 
        str_asymC = r"C energy (corrected with asym): #mu = " + str( round(HistCcorrected_Asymcut.GetFunction("gaus").GetParameter(1), 3) ) + r" #sigma = " + str( round(HistCcorrected_Asymcut.GetFunction("gaus").GetParameter(2), 3) ) 
        str_tdcC = r"C energy (corrected with timing): #mu = " + str( round(HistCcorrected_Tdc.GetFunction("gaus").GetParameter(1), 3) ) + r" #sigma = " + str( round(HistCcorrected_Tdc.GetFunction("gaus").GetParameter(2), 3) ) 
        #legend = ROOT.TLegend(0.15, 0.8, 0.4, 0.9); legend.AddEntry(HistCcorrected, "CPMT/containment", "l"); legend.AddEntry(HistCcorrected_Asymcut, "Asym Correction", "l");
        legend = ROOT.TLegend(0.15, 0.75, 0.6, 0.9);
        legend.AddEntry(HistCcorrected, str_rawC, "l"); 
        #legend.AddEntry(HistCcorrected_Asymcut, str_asymC, "l")
        legend.AddEntry(HistCcorrected_Tdc, str_tdcC, "l")
        legend.SetFillStyle(0)  # 0 means fully transparent
        legend.SetFillColor(0)   # Ensure no background color
        legend.SetBorderSize(0)   # Thin border around the legend
        legend.SetTextSize(0.03); legend.Draw("same")
        cerhist.SaveAs("Cerenehist{0}.png".format(energy))
  
        
        # Plot DR energy with asymmetry correction
        combcorrhist_Asymcut = ROOT.TCanvas("cCombEnergyCorrectedAsymmetry{0}".format(energy), "cCombEnergyCorrectedAsymmetry{0}".format(energy), 1400, 1200)
        HistCombCorrected.SetLineColor(ROOT.kBlack); HistCombCorrected_Asymcut.SetLineColor(ROOT.kGreen+1); HistCombCorrected.SetLineWidth(2); HistCombCorrected_Asymcut.SetLineWidth(2)
        HistCombCorrected_Asymcut.GetFunction("gaus").SetLineColor(ROOT.kBlue); HistCombCorrected_Asymcut.GetFunction("gaus").SetLineStyle(2); HistCombCorrected.GetFunction("gaus").SetLineColor(ROOT.kBlack)
        maximum = max(HistCombCorrected.GetMaximum(), HistCombCorrected_Asymcut.GetMaximum())
        HistCombCorrected.SetMaximum(maximum*1.2)
        HistCombCorrected.Draw()
        HistCombCorrected_Asymcut.Draw("same")
        str_rawDR = r"DR ene/containment (raw): #mu = " + str( round(HistCombCorrected.GetFunction("gaus").GetParameter(1), 3) ) + r" #sigma = " + str( round(HistCombCorrected.GetFunction("gaus").GetParameter(2), 3) ) 
        str_asymDR = r"DR energy (corrected with asym): #mu = " + str( round(HistCombCorrected_Asymcut.GetFunction("gaus").GetParameter(1), 3) ) + r" #sigma = " + str( round(HistCombCorrected_Asymcut.GetFunction("gaus").GetParameter(2), 3) ) 
        #legend = ROOT.TLegend(0.15, 0.8, 0.4, 0.9); legend.AddEntry(HistCombCorrected, "Raw PMT (DR), containment correction", "l"); legend.AddEntry(HistCombCorrected_Asymcut, "with asym Correction", "l"); legend.SetTextSize(0.015); legend.Draw("same")
        legend = ROOT.TLegend(0.15, 0.75, 0.6, 0.9); legend.AddEntry(HistCombCorrected, str_rawDR, "l"); legend.AddEntry(HistCombCorrected_Asymcut, str_asymDR, "l")
        legend.SetFillStyle(0)  # 0 means fully transparent
        legend.SetFillColor(0)   # Ensure no background color
        legend.SetBorderSize(0)   # Thin border around the legend
        legend.SetTextSize(0.03)
        legend.Draw("same")
        combcorrhist_Asymcut.SaveAs("DRcorrectedAsymmetryenehist{0}.png".format(energy))
        
        # Plot DR energy with TDC correction
        combcorrhist_tdc = ROOT.TCanvas("cCombEnergyCorrectedTdc{0}".format(energy), "cCombEnergyCorrectedTdc{0}".format(energy), 1400, 1200)
        HistCombCorrected.SetLineColor(ROOT.kBlack); HistCombCorrected_tdc.SetLineColor(ROOT.kGreen+2); HistCombCorrected.SetLineWidth(2); HistCombCorrected_tdc.SetLineWidth(2)
        HistCombCorrected_tdc.GetFunction("gaus").SetLineColor(ROOT.kGreen+2); HistCombCorrected_tdc.GetFunction("gaus").SetLineStyle(2); HistCombCorrected.GetFunction("gaus").SetLineColor(ROOT.kBlack)
        maximum = max(HistCombCorrected.GetMaximum(), HistCombCorrected_tdc.GetMaximum())
        HistCombCorrected.SetMaximum(maximum*1.2)     
        HistCombCorrected.Draw()
        HistCombCorrected_tdc.Draw("same")
        str_rawDR = r"DR ene/containment (raw): #mu = " + str( round(HistCombCorrected.GetFunction("gaus").GetParameter(1), 3) ) + r" #sigma = " + str( round(HistCombCorrected.GetFunction("gaus").GetParameter(2), 3) ) 
        str_tdcDR = r"DR energy (corrected with timing): #mu = " + str( round(HistCombCorrected_tdc.GetFunction("gaus").GetParameter(1), 3) ) + r" #sigma = " + str( round(HistCombCorrected_tdc.GetFunction("gaus").GetParameter(2), 3) ) 
        #legend = ROOT.TLegend(0.15, 0.8, 0.4, 0.9); legend.AddEntry(HistCombCorrected, "Raw PMT (DR), containment correction", "l"); legend.AddEntry(HistCombCorrected_Asymcut, "with tdc Correction", "l")
        legend = ROOT.TLegend(0.15, 0.75, 0.6, 0.9); legend.AddEntry(HistCombCorrected, str_rawDR, "l"); legend.AddEntry(HistCombCorrected_tdc, str_tdcDR, "l")
        legend.SetFillStyle(0)  # 0 means fully transparent
        legend.SetFillColor(0)   # Ensure no background color
        legend.SetBorderSize(0)   # Thin border around the legend
        legend.SetTextSize(0.03); legend.Draw("same")
        combcorrhist_tdc.SaveAs("DRcorrectedTdcenehist{0}.png".format(energy))



        # Plot Histgram with both corrections over the one with asymmetry only
        cBothCorrections = ROOT.TCanvas("cBothCorrections{0}".format(energy), "cBothCorrections{0}".format(energy), 1400, 1200)
        HistCombCorrected_Asymcut.GetFunction("gaus").SetLineColor(ROOT.kBlack); HistCombCorrected_Asymcut.GetFunction("gaus").SetLineWidth(2); HistCombCorrected_Asymcut.SetLineColor(ROOT.kBlack); HistCombCorrected_Asymcut.SetLineWidth(2)
        HistDR.GetFunction("gaus").SetLineColor(ROOT.kGreen+2); HistDR.GetFunction("gaus").SetLineWidth(2); HistDR.SetLineColor(ROOT.kGreen+1); HistDR.SetLineWidth(2)

        #HistCombCorrected_tdc.GetFunction("gaus").SetLineColor(ROOT.kGreen+2); HistCombCorrected_tdc.GetFunction("gaus").SetLineStyle(2); HistCombCorrected.GetFunction("gaus").SetLineColor(ROOT.kBlack)
        #HistDR.GetFunction("gaus").SetLineColor(ROOT.kGreen+2); HistDR.GetFunction("gaus").SetLineStyle(2); HistCombCorrected_tdc.GetFunction("gaus").SetLineColor(ROOT.kBlack)
        maximum = max(HistCombCorrected_Asymcut.GetMaximum(), HistDR.GetMaximum())
        HistCombCorrected_Asymcut.SetMaximum(maximum*1.15)
        #HistCombCorrected.SetMaximum(maximum*1.2)     
        HistCombCorrected_Asymcut.Draw()
        #HistCombCorrected_tdc.Draw("same")
        HistDR.Draw("same")
        legend = ROOT.TLegend(0.15, 0.8, 0.4, 0.9); legend.AddEntry(HistCombCorrected_Asymcut, "Asymmetry correction", "l"); legend.AddEntry(HistDR, "Asymmetry & TDC Correction", "l"); legend.SetTextSize(0.013); legend.Draw("same")
        cBothCorrections.SaveAs("DReneBothCorrections{0}.png".format(energy))

    
        #BestFitS_Asym = HistScorrected.GetFunction("gaus")
        #BestFitC_Asym = HistCcorrected.GetFunction("gaus")
        BestFitS_Asym = HistScorrected_Asymcut.GetFunction("gaus")
        BestFitC_Asym = HistCcorrected_Asymcut.GetFunction("gaus")
        BestFitComb = HistComb.GetFunction("gaus")
        BestFitCombCorrected = HistCombCorrected.GetFunction("gaus")
        BestFitCombCorrected_Asymcut = HistCombCorrected_Asymcut.GetFunction("gaus")
        BestFitCombCorrected_tdc = HistCombCorrected_tdc.GetFunction("gaus")
        BestFitS_Tdc = HistScorrected_Tdc.GetFunction("gaus")
        BestFitC_Tdc = HistCcorrected_Tdc.GetFunction("gaus")

        BestFitCombCorrected_AsymTdc = HistDR.GetFunction("gaus")
        BestFitCombAsym2 = HistDRAsym2.GetFunction("gaus")
        BestFitCombTdc2 = HistDRTdc2.GetFunction("gaus")



        ROOT.gStyle.SetOptStat(11111)
        ROOT.gStyle.SetOptFit(111)
        # Plot Histgram with dual-readout energies built after their corrections
        cDRAfterCorrection = ROOT.TCanvas("cDRAfterCorrections{0}".format(energy), "cDRAfterCorrections{0}".format(energy), 1400, 1200)
        #HistCombCorrected_Asymcut.GetFunction("gaus").SetLineColor(ROOT.kBlack); HistCombCorrected_Asymcut.GetFunction("gaus").SetLineWidth(2); HistCombCorrected_Asymcut.SetLineColor(ROOT.kBlack); HistCombCorrected_Asymcut.SetLineWidth(2)
        #HistDR.GetFunction("gaus").SetLineColor(ROOT.kGreen+2); HistDR.GetFunction("gaus").SetLineWidth(2); HistDR.SetLineColor(ROOT.kGreen+1); HistDR.SetLineWidth(2)
        #HistCombCorrected_tdc.GetFunction("gaus").SetLineColor(ROOT.kGreen+2); HistCombCorrected_tdc.GetFunction("gaus").SetLineStyle(2); HistCombCorrected.GetFunction("gaus").SetLineColor(ROOT.kBlack)
        #HistDR.GetFunction("gaus").SetLineColor(ROOT.kGreen+2); HistDR.GetFunction("gaus").SetLineStyle(2); HistCombCorrected_tdc.GetFunction("gaus").SetLineColor(ROOT.kBlack)
        #maximum = max(HistCombCorrected_Asymcut.GetMaximum(), HistDR.GetMaximum())
        #HistCombCorrected_Asymcut.SetMaximum(maximum*1.15)
        #HistCombCorrected.SetMaximum(maximum*1.2)  
        HistDRTdc2.SetLineColor(ROOT.kBlack); HistDRTdc2.SetLineWidth(2)
        HistDRTdc2.GetFunction("gaus").SetLineColor(ROOT.kGreen+1); HistDRTdc2.GetFunction("gaus").SetLineWidth(2)   
        HistDRTdc2.Draw()
        #HistCombCorrected_tdc.Draw("same")
        #HistDR.Draw("same")
        #legend = ROOT.TLegend(0.15, 0.8, 0.4, 0.9); legend.AddEntry(HistCombCorrected_Asymcut, "Asymmetry correction", "l"); legend.AddEntry(HistDR, "Asymmetry & TDC Correction", "l"); legend.SetTextSize(0.013); legend.Draw("same")
        cDRAfterCorrection.SaveAs("DReneAfterCorrections{0}.png".format(energy))



        # Uncomment to store resolution with Asymmetry Correction
        #MeanVec_S.append(BestFitS_Asym.GetParameter(1)); MeanErrVec_S.append(BestFitS_Asym.GetParError(1)); RmsVec_S.append(BestFitS_Asym.GetParameter(2)); RmsErrVec_S.append(BestFitS_Asym.GetParError(2))
        #MeanVec_C.append(BestFitC_Asym.GetParameter(1)); MeanErrVec_C.append(BestFitC_Asym.GetParError(1)); RmsVec_C.append(BestFitC_Asym.GetParameter(2)); RmsErrVec_C.append(BestFitC_Asym.GetParError(2))
        #MeanVec_Comb.append(BestFitCombCorrected_Asymcut.GetParameter(1)); MeanErrVec_Comb.append(BestFitCombCorrected_Asymcut.GetParError(1)); RmsVec_Comb.append(BestFitCombCorrected_Asymcut.GetParameter(2)); RmsErrVec_Comb.append(BestFitCombCorrected_Asymcut.GetParError(2))

        #MeanVec_Comb.append(BestFitCombCorrected.GetParameter(1)); MeanErrVec_Comb.append(BestFitCombCorrected.GetParError(1)); RmsVec_Comb.append(BestFitCombCorrected.GetParameter(2)); RmsErrVec_Comb.append(BestFitCombCorrected.GetParError(2))
        #MeanVec_Comb.append(BestFitComb.GetParameter(1)); MeanErrVec_Comb.append(BestFitComb.GetParError(1)); RmsVec_Comb.append(BestFitComb.GetParameter(2)); RmsErrVec_Comb.append(BestFitComb.GetParError(2))
        
        # Only append DR energy to be plotted
        #MeanVec_Comb.append(BestFitCombCorrected_Asymcut.GetParameter(1)); MeanErrVec_Comb.append(BestFitCombCorrected_Asymcut.GetParError(1)); RmsVec_Comb.append(BestFitCombCorrected_Asymcut.GetParameter(2)); RmsErrVec_Comb.append(BestFitCombCorrected_Asymcut.GetParError(2))


        # Uncomment to print TDC-corrected energy parameters
        MeanVec_S.append(BestFitS_Tdc.GetParameter(1)); MeanErrVec_S.append(BestFitS_Tdc.GetParError(1)); RmsVec_S.append(BestFitS_Tdc.GetParameter(2)); RmsErrVec_S.append(BestFitS_Tdc.GetParError(2))
        MeanVec_C.append(BestFitC_Tdc.GetParameter(1)); MeanErrVec_C.append(BestFitC_Tdc.GetParError(1)); RmsVec_C.append(BestFitC_Tdc.GetParameter(2)); RmsErrVec_C.append(BestFitC_Tdc.GetParError(2))
        #MeanVec_Comb.append(BestFitCombCorrected_tdc.GetParameter(1)); MeanErrVec_Comb.append(BestFitCombCorrected_tdc.GetParError(1)); RmsVec_Comb.append(BestFitCombCorrected_tdc.GetParameter(2)); RmsErrVec_Comb.append(BestFitCombCorrected_tdc.GetParError(2))
        # print energy parameters with both corrections
        #MeanVec_Comb.append(BestFitCombCorrected_AsymTdc.GetParameter(1)); MeanErrVec_Comb.append(BestFitCombCorrected_AsymTdc.GetParError(1)); RmsVec_Comb.append(BestFitCombCorrected_AsymTdc.GetParameter(2)); RmsErrVec_Comb.append(BestFitCombCorrected_AsymTdc.GetParError(2))
 
        #MeanVec_Comb.append(BestFitCombAsym2.GetParameter(1)); MeanErrVec_Comb.append(BestFitCombAsym2.GetParError(1)); RmsVec_Comb.append(BestFitCombAsym2.GetParameter(2)); RmsErrVec_Comb.append(BestFitCombAsym2.GetParError(2))
        MeanVec_Comb.append(BestFitCombTdc2.GetParameter(1)); MeanErrVec_Comb.append(BestFitCombTdc2.GetParError(1)); RmsVec_Comb.append(BestFitCombTdc2.GetParameter(2)); RmsErrVec_Comb.append(BestFitCombTdc2.GetParError(2))



        # Show effect of correction
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptFit(0)
        cEnergyPlot = ROOT.TCanvas("cEnePlot{0}".format(energy), "cEnePlot{0}".format(energy), 1400, 1200)
        cEnergyPlot.cd()
        cEnergyPlot.SetLeftMargin(0.12)
        integral = HistComb.Integral()
        scale_factor = 1/integral
        HistComb.Scale(scale_factor)
        integral = HistCombCorrected.Integral()
        scale_factor = 1/integral
        HistCombCorrected.Scale(scale_factor)
        integral = HistCombCorrected_Asymcut.Integral()
        scale_factor = 1/integral
        HistCombCorrected_Asymcut.Scale(scale_factor)
        integral = HistCombCorrected_tdc.Integral()
        scale_factor = 1/integral
        HistCombCorrected_tdc.Scale(scale_factor)
        # Get Maximum
        maximum = max(HistComb.GetMaximum(), HistCombCorrected.GetMaximum(), HistCombCorrected_Asymcut.GetMaximum(), HistCombCorrected_tdc.GetMaximum())
        # Plot raw DR energy
        HistComb.SetLineColor(ROOT.kBlue+1); HistComb.SetLineWidth(2)
        HistComb.SetLineWidth(2)     
        HistComb.SetMaximum(maximum*1.2)
        HistComb.SetTitle(r"Reco Energy ({0} GeV)".format(energy))     
        HistComb.GetYaxis().SetTitle("Normalised Counts")     
        HistComb.Draw("hist")
        # Plot DR energy / average containment
        HistCombCorrected.SetLineColor(ROOT.kRed+1); HistCombCorrected.SetLineWidth(2)
        HistCombCorrected.GetYaxis().SetTitle("Normalised Counts")     
        HistCombCorrected.SetLineWidth(2)
        HistCombCorrected.Draw("same hist")
        # Plot DR energy / average containment, with asymmetry correction -> attenuation
        HistCombCorrected_Asymcut.SetLineColor(ROOT.kGreen+2); HistCombCorrected_Asymcut.SetLineWidth(2)
        HistCombCorrected_Asymcut.SetLineWidth(2)
        HistCombCorrected_Asymcut.GetYaxis().SetTitle("Normalised Counts")     
        HistCombCorrected_Asymcut.Draw("same hist")
        # Plot DR energy / average containment, with tdc correction -> attenuation
        HistCombCorrected_tdc.SetLineColor(ROOT.kBlack); HistCombCorrected_Asymcut.SetLineWidth(2)
        HistCombCorrected_tdc.SetLineWidth(2)
        HistCombCorrected_tdc.GetYaxis().SetTitle("Normalised Counts")     
        HistCombCorrected_tdc.Draw("same hist")

        leg = ROOT.TLegend(0.53, 0.75, 0.97, 0.92)
        leg.SetTextSize(0.018)   
        #leg.AddEntry(HistComb, r"(S-#chi C)/(1-#chi)")
        leg.AddEntry(HistCombCorrected, r"(S-#chi C)/(1-#chi)/avg_cont")
        leg.AddEntry(HistCombCorrected_Asymcut, r"(S-#chi C)/(1-#chi)/avg_cont (Corrected for Asymmetry)")
        leg.AddEntry(HistCombCorrected_tdc, r"(S-#chi C)/(1-#chi)/avg_cont (Corrected for Timing)")
        leg.Draw() 
        cEnergyPlot.SaveAs("EnergyComparison_{0}GeV.png".format(energy))
        ROOT.gStyle.SetOptStat(1)
        ROOT.gStyle.SetOptFit(1)






    #DrawEnergyHist(dfCorrected_array, energies, binning_S, binning_C, "energyS", "energyC")

    # Get Noise value and make it an array
    noiseSvec = np.full(len(RmsVec_S), noiseS)
    RmsVec_S_corrected = np.sqrt(np.asarray(RmsVec_S)**2 - noiseSvec**2)
    print("S resolution before noise subtraction: ", RmsVec_S)
    print("S resolution after noise subtraction: ", RmsVec_S_corrected)

    noiseCvec = np.full(len(RmsVec_C), noiseC)
    RmsVec_C_corrected = np.sqrt(np.asarray(RmsVec_C)**2 - noiseCvec**2)
    print("C resolution before noise subtraction: ", RmsVec_C)
    print("C resolution after noise subtraction: ", RmsVec_C_corrected)

    # take pedestal S+C/2 and subtract it from combined energy
    noiseCombVec = np.full(len(RmsErrVec_C), (noiseS+noiseC)/2)
    RmsVec_Comb_corrected = np.sqrt(np.asarray(RmsVec_Comb)**2 - noiseCombVec**2)
    print("Combined resolution before noise subtraction: ", RmsVec_Comb)
    print("Combined resolution after noise subtraction: ", RmsVec_Comb_corrected)


    '''
    df = pd.DataFrame({
        "Energy" : energies,
        "Mean_S" : MeanVec_S,
        "RMS_S" : RmsVec_S_corrected,  
        "MeanErr_S" : MeanErrVec_S,  
        "RMSErr_S" : RmsErrVec_S,      
        "Mean_C" : MeanVec_C,     
        "RMS_C" : RmsVec_C_corrected,  
        "MeanErr_C" : MeanErrVec_C,  
        "RMSErr_C" : RmsErrVec_C,  
        "Mean_Comb" : MeanVec_Comb,
        "RMS_Comb" :  RmsVec_Comb_corrected,
        "MeanErr_Comb" : MeanErrVec_Comb,
        "RMSErr_Comb" : RmsErrVec_Comb
    })
    '''
    
        
    df = pd.DataFrame({
        "Energy" : energies,
        "Mean_S" : MeanVec_S,
        "RMS_S" : RmsVec_S,  
        "MeanErr_S" : MeanErrVec_S,  
        "RMSErr_S" : RmsErrVec_S,      
        "Mean_C" : MeanVec_C,     
        "RMS_C" : RmsVec_C,  
        "MeanErr_C" : MeanErrVec_C,  
        "RMSErr_C" : RmsErrVec_C,  
        "Mean_Comb" : MeanVec_Comb,
        "RMS_Comb" :  RmsVec_Comb,
        "MeanErr_Comb" : MeanErrVec_Comb,
        "RMSErr_Comb" : RmsErrVec_Comb
    })
    
                
        

    df.to_csv("DRcaloData.csv", index=False, sep="\t")

                
    
    








if __name__ == "__main__":
    main()
