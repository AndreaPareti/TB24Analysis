#############################################
### Analysis script for 2024 Test Beam  #####
### of the dual-readout prototype DRAGO #####
### A. Pareti, 14 March 2025            #####
### andrea.pareti01@ateneo.pv.it        #####
#############################################
# Simpy plot some dependencies of prototype responses to pions


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
chi_value = 0.35
#chi_value = 0.75

#containment = 0.875







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

    #data["pmtS_cont_att"] = data["totPMTSene"]/containment/S_attenuation_correction
    #data["pmtC_cont_att"] = data["totPMTCene"]/containment/C_attenuation_correction

    data["totDRene"] = (data["totPMTSene"]-chi_value*data["totPMTCene"]) / (1-chi_value)
    data["totDRene_cont"] = data["totDRene"]/containment
    #data["totDRene_cont_att"] = (data["totPMTSene"]/S_attenuation_correction-chi_value*data["totPMTCene"]/C_attenuation_correction) / (1-chi_value) / containment

    # input truth energy to dataframe
    data["TruthE"] = energy

    return data


def GetAsymProfiles(data, energy):

    # Profile normalised energy Vs Asymmetry
    DReneAsymSprof = ROOT.TProfile("DReneAsymSprof_{0}GeV".format(energy), r" DR Energy profile over Asymmetry(S) {0}GeV; TS24-TS21 / TS24 + TS21; E_DR/E".format(energy), 100, -1, 1)
    DReneAsymCprof = ROOT.TProfile("DReneAsymCprof_{0}GeV".format(energy), r"DR Energy profile over Asymmetry(C) {0}GeV; TC24-TC21 / TC24 + TC21; E_DR/E".format(energy), 100, -1, 1)

    # Profile normalised energy Vs Asymmetry
    SciAsymprof = ROOT.TProfile("sciprof_{0}GeV".format(energy), "totPMTSene profile over Asymmetry(S) {0}GeV; TS24-TS21 / TS24 + TS21; totPMTSene/E".format(energy), 100, -1, 1)
    CerAsymprof = ROOT.TProfile("cerprof_{0}GeV".format(energy), "totPMTCene profile over Asymmetry(C) {0}GeV; TC24-TC21 / TC24 + TC21; totPMTCene/E".format(energy), 100, -1, 1)

    for DRene, asymS,  asymC in zip(data["totDRene_cont"].values, data["AsymS"].values, data["AsymC"].values):
        #if(asymS < 0.9 and asymS > -0.9 and asymC < 0.9 and asymC > -0.9):
        DReneAsymSprof.Fill(asymS, DRene/energy)
        DReneAsymCprof.Fill(asymC, DRene/energy)

    for pmtS, pmtC, asymS, asymC, tdcS, tdcC  in zip(data["pmtS_cont"], data["pmtC_cont"], data["AsymS"], data["AsymC"], data["TDC_TS11"].values, data["TDC_TC11"].values):
        #if(asymS < 0.9 and asymS > -0.9 and asymC < 0.9 and asymC > -0.9):
        #SciAsymHist.Fill(asymS, pmtS)
        #CerAsymHist.Fill(asymC, pmtC)
        SciAsymprof.Fill(asymS, pmtS/energy)
        CerAsymprof.Fill(asymC, pmtC/energy)
        #if(tdcS > 512):
        #    ScieneStdcProf.Fill(tdcS, pmtS/energy)
        #if(tdcC>512):
        #    CereneCtdcProf.Fill(tdcC, pmtC/energy)

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





    '''
    cutCalib = 5000



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


    '''








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
    myCut = "(abs(XDWC2 - XDWC1) < 5) & (abs(YDWC2 - YDWC1)<5) & (totPMTSene>0) & (PShower<450) & (TailC<400) & (TailC>170) & (totLeakage<6500) & (MCounter<160) & (PShower>350) & (YDWC2>-20) & (YDWC2<5) & (XDWC2>-20) & (XDWC2<20)"
    # test
    #myCut = "(abs(XDWC2 - XDWC1) < 5) & (abs(YDWC2 - YDWC1)<5) & (totPMTSene>0) & (PShower<450) & (TailC<400) & (TailC>170) & (totLeakage<6500) & (MCounter<160) & (PShower>350) & (YDWC2>-20) & (YDWC2<5) & (XDWC2>-20) & (XDWC2<20) & (totLeakage>4500)"
    #myCut = "(TDC_TS00>600)"



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

 


    # First loop: prepare data and extract parametrizations for asymmetry and TDC information    
    for index, (run, energy, cont) in enumerate(zip(runs, energies, exp_containment)):
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


        #df, fDReneAsymS, fDReneAsymC, DReneAsymSprof, DReneAsymCprof,  fSciAsymProf, fCerAsymProf, SciAsymProf, CerAsymProf = GetDFparametrization(run, CurrentCut, filename, energy, cont)
        df, fDReneAsymS, fDReneAsymC, DReneAsymSprof, DReneAsymCprof,  fSciAsymProf, fCerAsymProf, SciAsymProf, CerAsymProf = GetAsymProfiles(df, energy)

        # clean events with TS11<100
        # cleanup for TDCs
        if(run=="1000"): 
            df["TDC_TS11"] = df["TDC_TS11"]-60; df["TDC_TC11"] = df["TDC_TC11"]-60
            df["TDC_TS15"] = df["TDC_TS15"]-60; df["TDC_TC15"] = df["TDC_TC15"]-60
            df["TDC_TS00"] = df["TDC_TS00"]-60; df["TDC_TC00"] = df["TDC_TC00"]-60

        # some cleanup on timing values
        df = df[df["TDC_TS00"]>600]
        df = df[df["TDC_TS11"]>600]
        df = df[df["TDC_TS15"]>600]
        #df = df[df["TDC_TC11"]>600]
        #df = df[df["TDC_TC15"]>600]
        #df = df[df["TDC_TC00"]>600]
        #df = df[df["TDC_TC11"]<675]
        df = df[df["TDC_TS11"]<720]
        #df = df[df["TDC_TC15"]<675]
        df = df[df["TDC_TS15"]<720]
        #df = df[df["TDC_TC00"]<675]
        df = df[df["TDC_TS00"]<720]


        # plot TDC profiles and return
        DReneStdcProf, DReneCtdcProf, ScieneStdcProf, CereneCtdcProf, fDReneStdc, fDReneCtdc, fScieneStdc, fCereneCtdc = GetTDCProfiles(df, energy)

        # Define some variable of interest 
        # Plot dependency of "leakage rings" with timing 
        df["leakRing1"] = df["L02"]+df["L04"]+df["L03"]
        df["leakRing2"] = df["L05"]+df["L07"]+df["L08"]+df["L09"]
        df["leakRing3"] = df["L10"]+df["L11"]+df["L12"]+df["L13"]
        df["leakRing4"] = df["L14"]+df["L15"]+df["L16"]+df["L20"]

        leakrings = ["leakRing1", "leakRing2", "leakRing3", "leakRing4"]
        for i, lRing in enumerate(leakrings):
            print("LeakRing{0}VsTdcTS11_{1}GeV.png".format(i+1, energy))
            # colz DR energy plot, energy corrected with both asymmetry and TDCs, over TDCs
            myOutfile = "LeakRing{0}ColzVsTdcTS11_{1}GeV.png".format(i+1, energy); labelPlot= "Leakage Ring{0} over TdcTS11, {1}GeV".format(i+1,energy)
            varX = df["TDC_TS11"]; varY = df[lRing]; title = "cColzLringTdcTS11"; labelX = "TS11 Tdc"; labelY = lRing
            nbinX = 50; xmin = 600; xmax = 720; nbinY = 50; ymin = df[lRing].min()*0.90; ymax = df[lRing].max()*1.1
            print(ymin, ymax, labelPlot, labelX, labelY)
            DrawColzPlot(myOutfile, title, varX, varY, labelPlot, labelX, labelY, nbinX, xmin, xmax, nbinY, ymin, ymax)
            myOutfile = "LeakRing{0}ProfVsTdcTS11_{1}GeV.png".format(i+1, energy); labelPlot= "Leakage Ring{0} over TdcTS11, {1}GeV".format(i+1,energy)
            varX = df["TDC_TS11"]; varY = df[lRing]; title = "cProfLringTdcTS11"; labelX = "TS11 Tdc"; labelY = lRing
            DrawProfPlot(myOutfile, title, varX, varY, labelPlot, labelX, labelY, nbinX, xmin, xmax)

            # Profile plot, energy corrected with asymmetry, over TDCs
            #myOutfile = "DReneProf_AsymCorrectedOverTdc_{0}GeV.png".format(energy); labelPlot= "DRene, asymmetry correction over Tdc TS11 {0} GeV".format(energy)
            ##varX = df["TDC_TS11"]; varY = df["energyDR_AsymCorrected"]/energy; title = "cProfTdcTS11"; labelX = "TS11 Tdc"; labelY = "Reco E_DR / E_beam"
            #nbinX = 50; xmin = 600; xmax = 720; nbinY = 50; ymin = energy*0.5; ymax = energy*1.5
            #DrawProfPlot(myOutfile, title, varX, varY, labelPlot, labelX, labelY, nbinX, xmin, xmax)


            #superimpose profiles on the same plot

        tdcs = ["TDC_TS11", "TDC_TS00", "TDC_TS15", "TDC_TC00", "TDC_TC11", "TDC_TC15"]
        ROOT.gStyle.SetOptStat(0)

        for tdcvar in tdcs:    

            ctitle = "cLeakProfiles_{1}_{0}".format(energy, tdcvar) 
            outfile = "LeakProfilesVs{1}_{0}GeV.png".format(energy, tdcvar)
            labelPlot = "LeakRings Profiles over {1}, {0}GeV".format(energy, tdcvar); labelX = "{0}".format(tdcvar); labelY = "Leakage Ring signal [ADC]"

            cCanva = ROOT.TCanvas(ctitle, "{0};{1};{2}".format(labelPlot, labelX, labelY), 1400, 1200)
            nbinX = 50; xmin = 600; xmax = 720
            Ring1Prof = ROOT.TProfile("Ring1Prof{1}{0}".format(energy, tdcvar), "{0}, {3}GeV;{1};{2}".format("Ring1Prof", labelX, labelY, energy), nbinX, xmin, xmax)
            Ring2Prof = ROOT.TProfile("Ring2Prof{1}{0}".format(energy, tdcvar), "{0}, {3}GeV;{1};{2}".format("Ring2Prof", labelX, labelY, energy), nbinX, xmin, xmax)
            Ring3Prof = ROOT.TProfile("Ring3Prof{1}{0}".format(energy, tdcvar), "{0}, {3}GeV;{1};{2}".format("Ring3Prof", labelX, labelY, energy), nbinX, xmin, xmax)
            Ring4Prof = ROOT.TProfile("Ring4Prof{1}{0}".format(energy, tdcvar), "{0}, {3}GeV;{1};{2}".format("Ring4Prof", labelX, labelY, energy), nbinX, xmin, xmax)

            for tdc, lring1, lring2, lring3, lring4 in zip(df[tdcvar], df["leakRing1"], df["leakRing2"], df["leakRing3"], df["leakRing4"]):
                Ring1Prof.Fill(tdc, lring1)
                Ring2Prof.Fill(tdc, lring2)
                Ring3Prof.Fill(tdc, lring3)
                Ring4Prof.Fill(tdc, lring4)

            #myProf.SetMarkerStyle(markerStyle); myProf.SetMarkerColor(markerColor)
            #myProf.Draw("colz")

            Ring1Prof.SetMarkerColor(ROOT.kBlack); Ring1Prof.SetMarkerStyle(20); Ring1Prof.SetMarkerSize(2)
            Ring2Prof.SetMarkerColor(ROOT.kRed); Ring2Prof.SetMarkerStyle(21); Ring2Prof.SetMarkerSize(2)
            Ring3Prof.SetMarkerColor(ROOT.kBlue); Ring3Prof.SetMarkerStyle(22); Ring3Prof.SetMarkerSize(2)
            Ring4Prof.SetMarkerColor(ROOT.kGreen+1); Ring4Prof.SetMarkerStyle(23); Ring3Prof.SetMarkerSize(2)

            myMin = min(Ring1Prof.GetMinimum(), Ring2Prof.GetMinimum(), Ring3Prof.GetMinimum(), Ring4Prof.GetMinimum())
            myMax = max(Ring1Prof.GetMaximum(), Ring2Prof.GetMaximum(), Ring3Prof.GetMaximum(), Ring4Prof.GetMaximum())
            Ring1Prof.Scale(1/myMax); Ring2Prof.Scale(1/myMax), Ring3Prof.Scale(1/myMax); Ring4Prof.Scale(1/myMax)

            #Ring1Prof.SetMinimum(600); Ring1Prof.SetMaximum(myMax*1.05)
            Ring1Prof.SetMinimum(0.3); Ring1Prof.SetMaximum(1.1)
            cCanva.SetLeftMargin(0.13)

            Ring1Prof.Draw(); Ring2Prof.Draw("same"); Ring3Prof.Draw("same"); Ring4Prof.Draw("same")
            legend = ROOT.TLegend(0.75, 0.75, 0.90, 0.90)
            legend.AddEntry(Ring1Prof, "Leakage Ring 1", "P")
            legend.AddEntry(Ring2Prof, "Leakage Ring 2", "P")
            legend.AddEntry(Ring3Prof, "Leakage Ring 3", "P")
            legend.AddEntry(Ring4Prof, "Leakage Ring 4", "P")
            legend.Draw()

            cCanva.SaveAs(outfile)









if __name__ == "__main__":
    main()
