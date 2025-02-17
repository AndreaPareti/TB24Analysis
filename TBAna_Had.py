import pandas as pd
import numpy as np
import ROOT 
import uproot


#infolder = "/home/storage/data_apareti/TB24/PionScan_OldHVsaturated/"
infolder = "/home/storage/data_apareti/TB24/PionScan/"
#infolder = "/home/storage/data_apareti/TB24/PionScanCorrected/pion24/"
pedfolder = "/home/storage/data_apareti/TB24/PedestalRuns/"
treename = "Ftree"
chi_value = 0.35
#chi_value = 0.75

containment = 0.875


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




def GetDFparametrization(run, Cut, filename, energy): 
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


    #data["totPMTDRene"] = (data["totPMTSene"]/S_attenuation_correction-chi_value*data["totPMTCene"]/C_attenuation_correction) / (1-chi_value)
    #data["totDRene"] = (data["pmtScorr"]/S_attenuation_correction-chi_value*data["totPMTCene"]/C_attenuation_correction) / (1-chi_value)
    #data["totDRene_corrected"] =  data["totDRene"]/containment

    # input truth energy to dataframe
    data["TruthE"] = energy


    # Profile normalised energy Vs Asymmetry
    DReneAsymSprof = ROOT.TProfile("DReneAsymSprof_{0}GeV".format(energy), r" DR Energy profile over Asymmetry(S) {0}GeV; TS24-TS21 / TS24 + TS21; E_DR/E".format(energy), 100, -1, 1)
    DReneAsymCprof = ROOT.TProfile("DReneAsymCprof_{0}GeV".format(energy), r"DR Energy profile over Asymmetry(C) {0}GeV; TC24-TC21 / TC24 + TC21; E_DR/E".format(energy), 100, -1, 1)

    # Profile normalised energy Vs Asymmetry
    SciAsymprof = ROOT.TProfile("sciprof_{0}GeV".format(energy), "totPMTSene profile over Asymmetry(S) {0}GeV; TS24-TS21 / TS24 + TS21; totPMTSene/E".format(energy), 100, -1, 1)
    CerAsymprof = ROOT.TProfile("cerprof_{0}GeV".format(energy), "totPMTCene profile over Asymmetry(C) {0}GeV; TC24-TC21 / TC24 + TC21; totPMTCene/E".format(energy), 100, -1, 1)

    # Show reco energy Vs Asymmetry
    #SciAsymHist = ROOT.TH2D("ScivsAsymHist_{0}GeV".format(energy), "ScivsAsymHist_{0}GeV; TS24-TS21/TS24+TS21; totPMTSene".format(energy), 50, -1.5, 1.5, 50, 0., 2*energy)
    #CerAsymHist = ROOT.TH2D("CervsAsymHist_{0}GeV".format(energy), "CervsAsymHist_{0}GeV; TSC4-TC21/TC24+TC21; totPMTCene".format(energy), 50, -1.5, 1.5, 50, 0., 2*energy)


    for DRene, asymS,  asymC in zip(data["totDRene_cont"].values, data["AsymS"].values, data["AsymC"].values):
        #if(asymS < 0.9 and asymS > -0.9 and asymC < 0.9 and asymC > -0.9):
        DReneAsymSprof.Fill(asymS, DRene/energy)
        DReneAsymCprof.Fill(asymC, DRene/energy)

    for pmtS, pmtC, asymS, asymC in zip(data["pmtS_cont"], data["pmtC_cont"], data["AsymS"], data["AsymC"]):
        #if(asymS < 0.9 and asymS > -0.9 and asymC < 0.9 and asymC > -0.9):
        #SciAsymHist.Fill(asymS, pmtS)
        #CerAsymHist.Fill(asymC, pmtC)
        SciAsymprof.Fill(asymS, pmtS/energy)
        CerAsymprof.Fill(asymC, pmtC/energy)

    #cSciAsym = ROOT.TCanvas("cSciAsym{}".format(energy), "cSciAsymS{}".format(energy), 1400, 1200)
    #SciAsymHist.Draw("colz")
    #cSciAsym.SaveAs("ScivsAysm{0}GeV.png".format(energy))    
    #cCerAsym = ROOT.TCanvas("cCerAsym{}".format(energy), "cCerAsym{}".format(energy), 1400, 1200)
    #CerAsymHist.Draw("colz")
    #cCerAsym.SaveAs("CervsAysm{0}GeV.png".format(energy))  


    # Fit with 5 degree polynomial
    DReneAsymSprof.Fit("pol5", "", "", -0.9, 0.9)
    DReneAsymCprof.Fit("pol5", "", "", -0.9, 0.9)
    # Get fitted function
    fDReneAsymS = DReneAsymSprof.GetFunction("pol5")
    fDReneAsymC = DReneAsymCprof.GetFunction("pol5")

    SciAsymprof.Fit("pol5", "", "", -0.9, 0.9)
    CerAsymprof.Fit("pol5", "", "", -0.9, 0.9) 
    fPMTSeneAsym = SciAsymprof.GetFunction("pol5")
    fPMTCeneAsym = CerAsymprof.GetFunction("pol5")
    ctotPMTeneAsymProf = ROOT.TCanvas("ctotPMTeneAsymProf{0}".format(energy),"PMT signal over Asymmetry ({0} GeV)".format(energy), 1400, 1200)
    SciAsymprof.SetMarkerColor(ROOT.kRed); SciAsymprof.SetMarkerStyle(20); fPMTSeneAsym.SetLineColor(ROOT.kRed)
    CerAsymprof.SetMarkerColor(ROOT.kBlue); CerAsymprof.SetMarkerStyle(20); fPMTCeneAsym.SetLineColor(ROOT.kBlue)
    SciAsymprof.Draw()
    CerAsymprof.Draw("same")
    fPMTSeneAsym.Draw("same"); fPMTCeneAsym.Draw("same")
    ctotPMTeneAsymProf.SaveAs("PMTsignalAsym_{0}GeV.png".format(energy))



    # vectorize it
    #fDReneAsymSvector = np.vectorize(fDReneAsymS.Eval)
    #fDReneAsymCvector = np.vectorize(fDReneAsymC.Eval)
    #return data
    return data, fDReneAsymS, fDReneAsymC, DReneAsymSprof, DReneAsymCprof, fPMTSeneAsym, fPMTCeneAsym, SciAsymprof, CerAsymprof



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


def main():
    print("Hello there")
    #runs = ["0714", "0715", "0716", "0717", "0718", "0721"]
    runs = ["0968", "0967", "0966", "0965", "0963", "0962"]

    energies = [20, 40, 60, 80, 100, 120]
    exp_containment = [0.865, 0.87, 0.875, 0.88, 0.885, 0.89]
    #exp_containment = [0.875, 0.875, 0.875, 0.875, 0.875, 0.875]
    myCut = "(abs(XDWC2 - XDWC1) < 5) & (abs(YDWC2 - YDWC1)<5) & (totPMTSene>0) & (PShower<500) & (TailC<400) & (totLeakage<7000) & (MCounter<150)"

    varProf = "YDWC2"
    cut_x_min = [-19.83, -16.74, -16.22, -15.95, -15.60, -16.12, -16.07, -15.50]
    cut_x_max = [23.90, 22.19, 23.27, 23.44, 24.27, 23.79, 23.63, 24.12]
    cut_y_min = [-26.54, -25.80, -26.15, -26.15, -26.39, -25.63, -25.63, -26.03]
    cut_y_max = [13.38, 10.89, 9.72, 9.50, 9.86, 10.89, 10.54, 10.17]

    # Declare Parametrization functions to fill in a first loop
    #DReneAsymFitS10 = ROOT.TF1(); DReneAsymFitS20 = ROOT.TF1(); DReneAsymFitS30 = ROOT.TF1(); DReneAsymFitS40 = ROOT.TF1(); DReneAsymFitS60 = ROOT.TF1(); DReneAsymFitS80 = ROOT.TF1(); DReneAsymFitS100 = ROOT.TF1(); DReneAsymFitS120 = ROOT.TF1()
    #DReneAsymFitC10 = ROOT.TF1(); DReneAsymFitC20 = ROOT.TF1(); DReneAsymFitC30 = ROOT.TF1(); DReneAsymFitC40 = ROOT.TF1(); DReneAsymFitC60 = ROOT.TF1(); DReneAsymFitC80 = ROOT.TF1(); DReneAsymFitC100 = ROOT.TF1(); DReneAsymFitC120 = ROOT.TF1()


    DReneAsymFitSvec = []; DReneProfSvec = []
    DReneAsymFitCvec = []; DReneProfCvec = [] 

    SciAsymFitvec = []; SciAsymProfvec = []
    CerAsymFitvec = []; CerAsymProfvec = [] 

    # Store ntuples as pandas dataframes
    #df10, df20, df30, df40, df60, df80, df100, df120
    dfs = []

    # array of dataframe with corrected energies
    dfCorrected_array = []

    # Arrays to store reco energy parameters
    MeanVec_S=[]; RmsVec_S=[]; MeanErrVec_S=[]; RmsErrVec_S=[] 
    MeanVec_C=[]; RmsVec_C=[]; MeanErrVec_C=[]; RmsErrVec_C=[] 
    MeanVec_Comb=[]; RmsVec_Comb=[]; MeanErrVec_Comb=[]; RmsErrVec_Comb=[] 
    MeanVec_Comb_Asym=[]; RmsVec_Comb_Asym=[]; MeanErrVec_Comb_Asym=[]; RmsErrVec_Comb_Asym=[] 
    
    for index, (run, energy, cont) in enumerate(zip(runs, energies, exp_containment)):
        filename = "physics_sps2024_run" + run + ".root"
        #print(filename, energy, cont)
        print("file: ", filename, "\tEnergy: ", energy, "\tExpected containment: ", cont)
        #ROOT.gStyle.SetOptFit(0)
        print("Current cuts on DWCs: ", cut_x_min[index], cut_x_max[index], cut_y_min[index], cut_y_max[index])
        #CurrentCut = myCut + " & (XDWC2 > {0}) & (XDWC2 < {1}) & (YDWC2 > {2}) & (YDWC2 < {3})".format(cut_x_min[index], cut_x_max[index], cut_y_min[index], cut_y_max[index]) 
        CurrentCut = myCut


        #df = GetDFparametrization(run, CurrentCut, filename, energy)
        df, fDReneAsymS, fDReneAsymC, DReneAsymSprof, DReneAsymCprof,  fSciAsymProf, fCerAsymProf, SciAsymProf, CerAsymProf = GetDFparametrization(run, CurrentCut, filename, energy)


        dfs.append(df)
        DrawFem(energy, df["totPMTSene"], df["totPMTCene"], "totPMTSene", "totPMTCene", "Electromagnetic Fraction, raw variables", "HistFemPMT")
        #DrawFem(energy, df["totDRene"], df["totPMTCene"], "DRene", "totPMTCene", "Electromagnetic Fraction, raw variables", "HistFemDR")

        
        DReneAsymFitSvec.append(fDReneAsymS); DReneAsymFitCvec.append(fDReneAsymC);  DReneProfSvec.append(DReneAsymSprof); DReneProfCvec.append(DReneAsymCprof)
        SciAsymFitvec.append(fSciAsymProf); CerAsymFitvec.append(fCerAsymProf); SciAsymProfvec.append(SciAsymProf); CerAsymProfvec.append(CerAsymProf)

        ROOT.gStyle.SetOptStat(0)
        ctest = ROOT.TCanvas("cEneProfAsymmetry{0}".format(energy), "cEneProfAsymmetry{0}".format(energy), 1400, 1200)
        ctest.SetLeftMargin(0.15)
        DReneAsymSprof.GetYaxis().SetTitle("PMT energy / Beam nominal energy")
        DReneAsymSprof.SetLineWidth(2); DReneAsymSprof.SetLineColor(ROOT.kRed); DReneAsymSprof.SetMarkerStyle(ROOT.kFullCircle); DReneAsymSprof.SetMarkerColor(ROOT.kRed)
        DReneAsymSprof.SetMinimum(0.5); DReneAsymSprof.SetMaximum(1.5); DReneAsymSprof.Draw()
        #ctestS.SaveAs("testS{0}GeV.png".format(energy))

        #ctestC = ROOT.TCanvas("c{0}".format(energy), "c{0}".format(energy), 1400, 1200)
        DReneAsymCprof.SetLineWidth(2); DReneAsymCprof.SetLineColor(ROOT.kBlue); DReneAsymCprof.SetMarkerStyle(ROOT.kFullCircle); DReneAsymCprof.SetMarkerColor(ROOT.kBlue)
        DReneAsymCprof.SetMinimum(0.8); DReneAsymCprof.SetMaximum(1.2); DReneAsymCprof.Draw("same")
        

        #leg = ROOT.TLegend(0.65, 0.15, 0.85, 0.25)
        leg = ROOT.TLegend(0.53, 0.75, 0.97, 0.92)

        leg.AddEntry(DReneAsymSprof, "DR ene profiled over S tower asymmetry", "L")
        leg.AddEntry(DReneAsymCprof, "DR ene profiled over C tower asymmetry", "L")
        leg.Draw()
        ctest.SaveAs("EnergyProfOverAsymmetry_{0}GeV.png".format(energy))   
        

    # comment since I'm only using 40 GeV parametrisation
    #fDRS20  = np.vectorize(DReneAsymFitS[0].Eval); fDRC20  = np.vectorize(DReneAsymFitC[0].Eval)
    fDRS40  = np.vectorize(DReneAsymFitSvec[1].Eval); fDRC40  = np.vectorize(DReneAsymFitCvec[1].Eval)
    #fDRS60  = np.vectorize(DReneAsymFitS[2].Eval); fDRC60  = np.vectorize(DReneAsymFitC[2].Eval)
    #fDRS80  = np.vectorize(DReneAsymFitS[3].Eval); fDRC80  = np.vectorize(DReneAsymFitC[3].Eval)
    #fDRS100  = np.vectorize(DReneAsymFitS[4].Eval); fDRC100  = np.vectorize(DReneAsymFitC[4].Eval)
    #fDRS120  = np.vectorize(DReneAsymFitS[5].Eval); fDRC120  = np.vectorize(DReneAsymFitC[5].Eval)
   
    fSci40  = np.vectorize(SciAsymFitvec[1].Eval); fCer40  = np.vectorize(CerAsymFitvec[1].Eval)

    binning_S = [0, 30, 60, 90]
    binning_C = [0, 30, 60, 90]
    DrawEnergyHist(dfs, energies, binning_S, binning_C, "totPMTSene", "totPMTCene")

    noiseS, noiseC = GetPMTnoise()
    ROOT.gStyle.SetOptStat(1)
        
    for index, (energy, cont) in enumerate(zip(energies, exp_containment)):
        print("Index ", index, "\tUsing file with energy: ", energy)

        data = dfs[index]

        # Define conditions
        #binning_S = [0, 15, 25, 35, 50, 75, 100, 120]
        #binning_C = [0, 15, 25, 35, 50, 75, 100, 120]
        binning_S = [0, 30, 60, 90]
        binning_C = [0, 30, 60, 90]


        # Define binnings on PMT energy
        conditions_DR = [
            (data["totDRene_cont"] >= binning_S[0]) & (data["totDRene_cont"] < binning_S[1]),
            (data["totDRene_cont"] >= binning_S[1]) & (data["totDRene_cont"] < binning_S[2]),
            (data["totDRene_cont"] >= binning_S[2]) & (data["totDRene_cont"] < binning_S[3]),
            (data["totDRene_cont"] >= binning_S[3])
        ]

        #conditions_C = [
        #    (data["totPMTCene"] >= binning_C[0]) & (data["totPMTCene"] < binning_C[1]),
        #    (data["totPMTCene"] >= binning_C[1]) & (data["totPMTCene"] < binning_C[2]),
        #    (data["totPMTCene"] >= binning_C[2]) & (data["totPMTCene"] < binning_C[3]),
        #    (data["totPMTCene"] >= binning_C[3])
        #]

        # Associate a parametrisation to each energy point
        # A parametrisation at each energy is extracted, it is possible to use only one for all
        # or one per bin
        choices_DR = [
            data["totDRene"] / fDRS40(data["AsymS"])/cont,
            data["totDRene"] / fDRS40(data["AsymS"])/cont,
            data["totDRene"] / fDRS40(data["AsymS"])/cont,
            data["totDRene"] / fDRS40(data["AsymS"])/cont
        ]
        #choices_C = [
        #    data["totPMTCene"] / fC60(data["AsymC"]),
        #    data["totPMTCene"] / fC60(data["AsymC"]),
        #    data["totPMTCene"] / fC60(data["AsymC"]),
        #    data["totPMTCene"] / fC60(data["AsymC"])
        #]


        # Use np.select to apply the conditions
        data["energyDR_AsymCorrected"] = np.select(conditions_DR, choices_DR, default=data["totDRene_cont"])
        data["energyS_AsymCorrected"] = data["totPMTSene"]/fSci40(data["AsymS"])/cont
        data["energyC_AsymCorrected"] = data["totPMTCene"]/fCer40(data["AsymC"])/cont
        #data["energyC"] = np.select(conditions_C, choices_C, default=data["totPMTCene"])

        # Evaluate function on Asym variable and use it to correct energy
        #data["energyS"] = data["totPMTSene"]/fS20(data["AsymS"])
        #data["energyC"] = data["totPMTCene"]/fC20(data["AsymC"])
        DrawFem(energy, data["energyS_AsymCorrected"], data["energyC_AsymCorrected"], "S (corrected)", "C (corrected)", "Electromagnetic Fraction, Corrected for asymmetry", "HistFemCorrected")



        # Fill histograms with corrected energy values
        #HistScorrected = ROOT.TH1D("HistScorrected_{0}GeV".format(energy), "S Energy (corrected) {0}GeV; E [GeV]; Counts".format(energy), 100, energy-0.4*energy, energy+0.4*energy)
        #HistSraw = ROOT.TH1D("HistSraw_{0}GeV".format(energy), "totPMTSene {0}GeV; E [GeV]; Counts".format(energy), 100, energy-0.4*energy, energy+0.4*energy)
        #HistScorrected_Asymcut = ROOT.TH1D("HistScorrected_AsymCut_{0}GeV".format(energy), "S Energy (corrected, Asym cut) {0}GeV; E [GeV]; Counts".format(energy), 100, energy-0.4*energy, energy+0.4*energy)

        # Fill histograms with corrected energy values
        #HistCcorrected = ROOT.TH1D("HistCcorrected_{0}GeV".format(energy), "C Energy (corrected) {0}GeV; E [GeV]; Counts".format(energy), 100, energy-0.4*energy, energy+0.4*energy)
        #HistCraw = ROOT.TH1D("HistCraw_{0}GeV".format(energy), "totPMTCene {0}GeV; E [GeV]; Counts".format(energy), 100, energy-0.4*energy, energy+0.4*energy)
        #HistCcorrected_Asymcut = ROOT.TH1D("HistCcorrected_AsymCut_{0}GeV".format(energy), "C Energy (corrected, Asym cut) {0}GeV; E [GeV]; Counts".format(energy), 100, energy-0.4*energy, energy+0.4*energy)

        # Fill histograms with corrected energy values
        #HistCombCorrected = ROOT.TH1D("HistCombcorrected_{0}GeV".format(energy), "Combined Energy (corrected) {0}GeV; E [GeV]; Counts".format(energy), 100, energy-0.4*energy, energy+0.4*energy)


        #for eneS, pmtS, eneC, pmtC in zip(data["energyS"].values, data["totPMTSene"].values, data["energyC"].values, data["totPMTCene"].values):
        #    HistScorrected.Fill(eneS)
        #    HistSraw.Fill(pmtS)
        #    HistCcorrected.Fill(eneC)
        #    HistCraw.Fill(pmtC)


        #dfCorrected_array.append(data)
        # filter tree using only events with an asymmetry within range
        AsymCut = 0.5
        #data = data[ (np.abs(data["AsymS"]<AsymCut) ) & (np.abs(data["AsymC"]<AsymCut) ) ]
        #data_filtered = data[ (np.abs(data["AsymS"])<AsymCut ) & (np.abs(data["AsymC"])<AsymCut ) & (np.abs(data["BaryS"])<4) & (np.abs(data["BaryC"])<4) ]
        data_filtered = data[ (np.abs(data["AsymS"])<AsymCut ) & (np.abs(data["AsymC"])<AsymCut ) ]

        #for eneS, eneC in zip(data_filtered["energyS"].values, data_filtered["energyC"].values):
        #    HistScorrected_Asymcut.Fill(eneS)
        #    HistCcorrected_Asymcut.Fill(eneC)
        #    HistCombCorrected.Fill( (eneS+eneC)/2 )



        #for index, energy in enumerate(energies):
        #data = dfs[index]


        # Fill histograms with corrected energy values
        HistScorrected = ROOT.TH1D("HistScorrected_{0}GeV".format(energy), "totPMTSene/avg_containment ({0}GeV); E [GeV]; Counts".format(energy), 100, energy-0.7*energy, energy+0.7*energy)
        HistSraw = ROOT.TH1D("HistSraw_{0}GeV".format(energy), "totPMTSene ({0}GeV); E [GeV]; Counts".format(energy), 100, energy-0.7*energy, energy+0.7*energy)
        HistScorrected_Asymcut = ROOT.TH1D("HistSAsymcorrected_{0}GeV".format(energy), "totPMTSene/avg_containment (Corrected, {0}GeV); E [GeV]; Counts".format(energy), 100, energy-0.7*energy, energy+0.7*energy)

        # Fill histograms with corrected energy values
        HistCcorrected = ROOT.TH1D("HistCcorrected_{0}GeV".format(energy), "totPMTCene/avg_containment ({0}GeV); E [GeV]; Counts".format(energy), 100, energy-0.7*energy, energy+0.7*energy)
        HistCraw = ROOT.TH1D("HistCraw_{0}GeV".format(energy), "totPMTCene ({0}GeV); E [GeV]; Counts".format(energy), 100, energy-0.7*energy, energy+0.7*energy)
        HistCcorrected_Asymcut = ROOT.TH1D("HistCAsymcorrected_{0}GeV".format(energy), "totPMTCene/avg_containment (Corrected, {0}GeV); E [GeV]; Counts".format(energy), 100, energy-0.7*energy, energy+0.7*energy)

        # Fill histograms with corrected energy values
        HistCombCorrected = ROOT.TH1D("HistCombcorrected_{0}GeV".format(energy), r"(S-#chi C)/(1-#chi)/avg_containment  ({0}GeV); E [GeV]; Counts".format(energy), 100, energy-0.7*energy, energy+0.7*energy)
        HistComb = ROOT.TH1D("HistComb_{0}GeV".format(energy), r"(S-#chi C)/(1-#chi) ({0}GeV); E [GeV]; Counts".format(energy), 100, energy-0.7*energy, energy+0.7*energy)
        HistCombCorrected_Asymcut = ROOT.TH1D("HistDRene_AsymCorrected_{0}GeV".format(energy), r"(S-#chi C)/(1-#chi)/avg_containment (Corrected, {0}GeV); E [GeV]; Counts".format(energy), 100, energy-0.7*energy, energy+0.7*energy)


        for pmtS, pmtC, eneDR_cont, eneDR_corrected in zip(data["pmtS_cont"].values, data["pmtC_cont"].values, data["totDRene_cont"].values, data_filtered["energyDR_AsymCorrected"]):
            HistSraw.Fill(pmtS)
            HistCraw.Fill(pmtC)
            HistComb.Fill(eneDR_cont)
            HistCombCorrected.Fill(eneDR_cont)
            HistCombCorrected_Asymcut.Fill(eneDR_corrected)

        for eneS, eneC in zip(data_filtered["energyS_AsymCorrected"], data_filtered["energyC_AsymCorrected"]):
            HistScorrected_Asymcut.Fill(eneS)
            HistCcorrected_Asymcut.Fill(eneC)


        # Show reco energy Vs Asymmetry
        RecoEvsAsymSHist = ROOT.TH2D("RecoEvsAsymSHist_{0}GeV".format(energy), "RecoEvsAsymSHist_{0}GeV; TS24-TS21/TS24+TS21; Reco E_dr".format(energy), 50, -1.5, 1.5, 50, 0., 2*energy)
        RecoEvsAsymCHist = ROOT.TH2D("RecoEvsAsymCHist_{0}GeV".format(energy), "RecoEvsAsymCHist_{0}GeV; TSC4-TC21/TC24+TC21; Reco E_dr".format(energy), 50, -1.5, 1.5, 50, 0., 2*energy)

        for recoE, ts24, ts21, tc24, tc21 in zip(data["totDRene_cont"].values, data["TS24"].values, data["TS21"].values, data["TC24"].values, data["TC21"].values):
            RecoEvsAsymSHist.Fill(((ts24-ts21)/(ts24+ts21)), recoE)
            RecoEvsAsymCHist.Fill(((ts24-ts21)/(tc24+tc21)), recoE)
        cAsymS = ROOT.TCanvas("cAsymS{}".format(energy), "cAsymS{}".format(energy), 1400, 1200)
        RecoEvsAsymSHist.Draw("colz")
        cAsymS.SaveAs("RecoEvsAysmS{0}GeV.png".format(energy))    
        cAsymC = ROOT.TCanvas("cAsymC{}".format(energy), "cAsymC{}".format(energy), 1400, 1200)
        RecoEvsAsymCHist.Draw("colz")
        cAsymC.SaveAs("RecoEvsAysmC{0}GeV.png".format(energy))    


            

        HistSraw.Fit("gaus", "Q")
        totPMTSeneFit = HistSraw.GetFunction("gaus")
        HistScorrected_Asymcut.Fit("gaus", "Q")
        SeneFit = HistScorrected_Asymcut.GetFunction("gaus")

        HistCraw.Fit("gaus", "Q")
        HistCcorrected_Asymcut.Fit("gaus", "Q")
        totPMTCeneFit = HistCraw.GetFunction("gaus")
        CeneFit = HistCcorrected_Asymcut.GetFunction("gaus")

        HistComb.Fit("gaus","Q")
        DReneFit = HistComb.GetFunction("gaus")

        HistCombCorrected.Fit("gaus","Q")
        DReneFitCorrected = HistCombCorrected.GetFunction("gaus")

        HistCombCorrected_Asymcut.Fit("gaus","Q")
        DReneFitCorrected_Asymcut = HistCombCorrected_Asymcut.GetFunction("gaus")


        ROOT.gStyle.SetOptFit(111)
        # refit within -1.5 and +3 sigma
        HistSraw.Fit("gaus", "Q", "", totPMTSeneFit.GetParameter(1)-1.5*totPMTSeneFit.GetParameter(2), totPMTSeneFit.GetParameter(1)+1.5*totPMTSeneFit.GetParameter(2)) 
        HistCraw.Fit("gaus", "Q", "", totPMTCeneFit.GetParameter(1)-1.5*totPMTCeneFit.GetParameter(2), totPMTCeneFit.GetParameter(1)+1.5*totPMTCeneFit.GetParameter(2)) 
        HistComb.Fit("gaus", "Q", "", DReneFit.GetParameter(1)-1.5*DReneFit.GetParameter(2), DReneFit.GetParameter(1)+1.5*DReneFit.GetParameter(2)) 
        HistCombCorrected.Fit("gaus", "Q", "", DReneFitCorrected.GetParameter(1)-1.5*DReneFitCorrected.GetParameter(2), DReneFitCorrected.GetParameter(1)+1.5*DReneFitCorrected.GetParameter(2)) 
        HistCombCorrected_Asymcut.Fit("gaus", "Q", "", DReneFitCorrected_Asymcut.GetParameter(1)-1.5*DReneFitCorrected_Asymcut.GetParameter(2), DReneFitCorrected_Asymcut.GetParameter(1)+1.5*DReneFitCorrected_Asymcut.GetParameter(2)) 
        HistScorrected_Asymcut.Fit("gaus", "Q", "", SeneFit.GetParameter(1)-1.5*SeneFit.GetParameter(2), SeneFit.GetParameter(1)+1.5*SeneFit.GetParameter(2)) 
        HistCcorrected_Asymcut.Fit("gaus", "Q", "", CeneFit.GetParameter(1)-1.5*CeneFit.GetParameter(2), CeneFit.GetParameter(1)+1.5*CeneFit.GetParameter(2)) 



        scihist = ROOT.TCanvas("cSciEnergy{0}".format(energy), "cSciEnergy{0}".format(energy), 1400, 1200)
        HistSraw.SetLineColor(ROOT.kRed); HistScorrected_Asymcut.SetLineColor(ROOT.kRed+2)
        HistSraw.Draw()
        HistScorrected_Asymcut.Draw("same")
        scihist.SaveAs("Scienehist{0}.png".format(energy))
        cerhist = ROOT.TCanvas("cCerEnergy{0}".format(energy), "cCerEnergy{0}".format(energy), 1400, 1200)
        HistCraw.SetLineColor(ROOT.kBlue); HistCcorrected_Asymcut.SetLineColor(ROOT.kBlue+2)
        HistCraw.Draw()
        HistCcorrected_Asymcut.Draw("same")
        cerhist.SaveAs("Cerenehist{0}.png".format(energy))
        combhist = ROOT.TCanvas("cCombEnergy{0}".format(energy), "cCombEnergy{0}".format(energy), 1400, 1200)
        HistComb.Draw()
        combhist.SaveAs("Combenehist{0}.png".format(energy))
        combcorrhist = ROOT.TCanvas("cCombEnergyCorrected{0}".format(energy), "cCombEnergyCorrected{0}".format(energy), 1400, 1200)
        HistCombCorrected.Draw()
        combcorrhist.SaveAs("DRcorrectedenehist{0}.png".format(energy))
        combcorrhist_Asymcut = ROOT.TCanvas("cCombEnergyCorrectedAsymmetry{0}".format(energy), "cCombEnergyCorrectedAsymmetry{0}".format(energy), 1400, 1200)
        HistCombCorrected_Asymcut.Draw()
        combcorrhist_Asymcut.SaveAs("DRcorrectedAsymmetryenehist{0}.png".format(energy))


    
        #BestFitS = HistSraw.GetFunction("gaus")
        #BestFitC = HistCraw.GetFunction("gaus")
        BestFitS = HistScorrected_Asymcut.GetFunction("gaus")
        BestFitC = HistCcorrected_Asymcut.GetFunction("gaus")
        BestFitComb = HistComb.GetFunction("gaus")
        BestFitCombCorrected = HistCombCorrected.GetFunction("gaus")
        BestFitCombCorrected_Asymcut = HistCombCorrected_Asymcut.GetFunction("gaus")


        MeanVec_S.append(BestFitS.GetParameter(1)); MeanErrVec_S.append(BestFitS.GetParError(1)); RmsVec_S.append(BestFitS.GetParameter(2)); RmsErrVec_S.append(BestFitS.GetParError(2))
        MeanVec_C.append(BestFitC.GetParameter(1)); MeanErrVec_C.append(BestFitC.GetParError(1)); RmsVec_C.append(BestFitC.GetParameter(2)); RmsErrVec_C.append(BestFitC.GetParError(2))
        #MeanVec_Comb.append(BestFitCombCorrected.GetParameter(1)); MeanErrVec_Comb.append(BestFitCombCorrected.GetParError(1)); RmsVec_Comb.append(BestFitCombCorrected.GetParameter(2)); RmsErrVec_Comb.append(BestFitCombCorrected.GetParError(2))
        #MeanVec_Comb.append(BestFitComb.GetParameter(1)); MeanErrVec_Comb.append(BestFitComb.GetParError(1)); RmsVec_Comb.append(BestFitComb.GetParameter(2)); RmsErrVec_Comb.append(BestFitComb.GetParError(2))
        #MeanVec_Comb_Asym.append(BestFitCombCorrected_Asymcut.GetParameter(1)); MeanErrVec_Comb_Asym.append(BestFitCombCorrected_Asymcut.GetParError(1)); RmsVec_Comb_Asym.append(BestFitCombCorrected_Asymcut.GetParameter(2)); RmsErrVec_Comb_Asym.append(BestFitCombCorrected_Asymcut.GetParError(2))
        # Only append DR energy to be plotted
        MeanVec_Comb.append(BestFitCombCorrected_Asymcut.GetParameter(1)); MeanErrVec_Comb.append(BestFitCombCorrected_Asymcut.GetParError(1)); RmsVec_Comb.append(BestFitCombCorrected_Asymcut.GetParameter(2)); RmsErrVec_Comb.append(BestFitCombCorrected_Asymcut.GetParError(2))


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
        # Get Maximum
        maximum = max(HistComb.GetMaximum(), HistCombCorrected.GetMaximum(), HistCombCorrected_Asymcut.GetMaximum())
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
        leg = ROOT.TLegend(0.53, 0.75, 0.97, 0.92)
        leg.SetTextSize(0.018)   
        #leg.AddEntry(HistComb, r"(S-#chi C)/(1-#chi)")
        leg.AddEntry(HistCombCorrected, r"(S-#chi C)/(1-#chi)/avg_cont")
        leg.AddEntry(HistCombCorrected_Asymcut, r"(S-#chi C)/(1-#chi)/avg_cont (Corrected for Asymmetry)")
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
    
        

    df.to_csv("DRcaloData.csv", index=False, sep="\t")

                
    
    








if __name__ == "__main__":
    main()
