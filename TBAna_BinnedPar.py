import pandas as pd
import numpy as np
import ROOT 
import uproot


infolder = "/home/storage/data_apareti/TB24/EnergyScan/"
treename = "Ftree"



def GetDFparametrization(run, Cut, filename, energy): 
    print("Using file: ", filename, energy)
    #ROOT.gStyle.SetOptFit(0)
    #print("Current cuts on DWCs: ", cut_x_min[index], cut_x_max[index], cut_y_min[index], cut_y_max[index])
    #CurrentCut = myCut + " & (XDWC2 > {0}) & (XDWC2 < {1}) & (YDWC2 > {2}) & (YDWC2 < {3})".format(cut_x_min[index], cut_x_max[index], cut_y_min[index], cut_y_max[index]) 

    root_file = uproot.open(infolder+filename)
    tree = root_file[treename]

    data = tree.arrays(cut=Cut, library="pd")

    # define Asymmetry variable 
    data["AsymS"] = (data["TS11"] - data["TS15"]) / (data["TS11"] + data["TS15"] )
    data["AsymC"] = (data["TC11"] - data["TC15"]) / (data["TC11"] + data["TC15"] )

    data["BaryS"] = (data["TS00"]-28.3*data["TS11"]+28.3*data["TS15"])/(data["TS00"]+data["TS11"]+data["TS15"])
    data["BaryC"] = (data["TC00"]-28.3*data["TC11"]+28.3*data["TC15"])/(data["TC00"]+data["TC11"]+data["TC15"])

    # input truth energy to dataframe
    data["TruthE"] = energy


    # Profile normalised energy Vs Asymmetry
    eneSprof = ROOT.TProfile("eneSprof_{0}GeV".format(energy), "S Energy profile over Asymmetry {0}GeV; TS11-TS15 / TS11 + TS15; totPMTSene/E".format(energy), 100, -1, 1)
    eneCprof = ROOT.TProfile("eneCprof_{0}GeV".format(energy), "C Energy profile over Asymmetry {0}GeV; TC11-TC15 / TC11 + TC15; totPMTCene/E".format(energy), 100, -1, 1)
    eneSprof.SetDirectory(0)
    eneCprof.SetDirectory(0)

    for eneS, asymS, eneC, asymC in zip(data["totPMTSene"].values, data["AsymS"].values, data["totPMTCene"].values, data["AsymC"].values):
        eneSprof.Fill(asymS, eneS/energy)
        eneCprof.Fill(asymC, eneC/energy)

    # Fit with 5 degree polynomial
    eneSprof.Fit("pol5", "", "", -0.9, 0.9)
    eneCprof.Fit("pol5", "", "", -0.9, 0.9)

    # Get fitted function
    funcS = eneSprof.GetFunction("pol5")
    funcC = eneCprof.GetFunction("pol5")


    # vectorize it
    funcSvector = np.vectorize(funcS.Eval)
    funcCvector = np.vectorize(funcC.Eval)


    #return data, funcSvector, funcCvector
    return data, funcS, funcC, eneSprof, eneCprof


def DrawEnergyHist(dfs, energies, binning_S, binning_C, varname_S, varname_C):
    ROOT.gStyle.SetOptStat(0)
    c_histS = ROOT.TCanvas("c_histS", "c_histS", 1400, 1200)
    c_histS.SetRightMargin(0.05)
    c_histS.SetLeftMargin(0.15)
    histlistS = []; histlistC = []
    for i, (df, energy) in enumerate(zip(dfs, energies)):
        histS = ROOT.TH1D("histS{0}".format(energy), "histS{0}; E [GeV]; Normalized Counts".format(energy), 160, 0., 160)
        histC = ROOT.TH1D("histC{0}".format(energy), "histC{0}; E [GeV]; Normalized Counts".format(energy), 160, 0., 160)

        for varS, varC in zip(df[varname_S], df[varname_C]): 
            histS.Fill(varS)
            histC.Fill(varC)
        histS.SetLineColor(i+1)
        histC.SetLineColor(i+1)

        histlistS.append(histS)
        histlistC.append(histC)


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
    #line1 = ROOT.TLine(binning_S[1], histlistS[0].GetYaxis().GetXmin(), binning_S[1], histlistS[0].GetYaxis().GetXmax())
    #line1.Draw("same")

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
    
    




def main():
    print("Hello there")
    runs = ["0786", "0766", "0772", "0774", "0775", "0778", "0779", "0792"]
    energies = [10, 20, 30, 40, 60, 80, 100, 120]
    myCut = "(abs(XDWC2 - XDWC1) < 5) & (abs(YDWC2 - YDWC1)<5) & (MCounter<200) & (TailC<300) & (C2>160) & ( (totLeakage - L20)<5000 ) & (PShower>550)"
    #myCut = "(abs(XDWC2 - XDWC1) < 5) & (abs(YDWC2 - YDWC1)<5) & (MCounter<200) & (TailC<300) & (C2>160) & ( (totLeakage - L20)<5000 )"

    varProf = "YDWC2"
    cut_x_min = [-19.83, -16.74, -16.22, -15.95, -15.60, -16.12, -16.07, -15.50]
    cut_x_max = [23.90, 22.19, 23.27, 23.44, 24.27, 23.79, 23.63, 24.12]
    cut_y_min = [-26.54, -25.80, -26.15, -26.15, -26.39, -25.63, -25.63, -26.03]
    cut_y_max = [13.38, 10.89, 9.72, 9.50, 9.86, 10.89, 10.54, 10.17]

    # Declare Parametrization functions to fill in a first loop

    FitS10 = ROOT.TF1(); FitS20 = ROOT.TF1(); FitS30 = ROOT.TF1(); FitS40 = ROOT.TF1(); FitS60 = ROOT.TF1(); FitS80 = ROOT.TF1(); FitS100 = ROOT.TF1(); FitS120 = ROOT.TF1()
    FitC10 = ROOT.TF1(); FitC20 = ROOT.TF1(); FitC30 = ROOT.TF1(); FitC40 = ROOT.TF1(); FitC60 = ROOT.TF1(); FitC80 = ROOT.TF1(); FitC100 = ROOT.TF1(); FitC120 = ROOT.TF1()


    FitS = []; profS = []
    FitC = []; profC = [] 

    # Store ntuples as pandas dataframes
    #df10, df20, df30, df40, df60, df80, df100, df120
    dfs = []

    # array of dataframe with corrected energies
    dfCorrected_array = []

    # Arrays to store reco energy parameters
    MeanVec_S=[]; RmsVec_S=[]; MeanErrVec_S=[]; RmsErrVec_S=[] 
    MeanVec_C=[]; RmsVec_C=[]; MeanErrVec_C=[]; RmsErrVec_C=[] 
    MeanVec_Comb=[]; RmsVec_Comb=[]; MeanErrVec_Comb=[]; RmsErrVec_Comb=[] 
    
    for index, (run, energy) in enumerate(zip(runs, energies)):
        filename = "physics_sps2024_run" + run + ".root"
        print(filename, energy)
        ROOT.gStyle.SetOptFit(0)
        print("Current cuts on DWCs: ", cut_x_min[index], cut_x_max[index], cut_y_min[index], cut_y_max[index])
        CurrentCut = myCut + " & (XDWC2 > {0}) & (XDWC2 < {1}) & (YDWC2 > {2}) & (YDWC2 < {3})".format(cut_x_min[index], cut_x_max[index], cut_y_min[index], cut_y_max[index]) 

        df, funcS, funcC, eneSprof, eneCprof = GetDFparametrization(run, CurrentCut, filename, energy)

        dfs.append(df)
        FitS.append(funcS)
        FitC.append(funcC)
        profS.append(eneSprof)
        profC.append(eneCprof)


    fS10  = np.vectorize(FitS[0].Eval); fC10  = np.vectorize(FitC[0].Eval)  
    fS20  = np.vectorize(FitS[1].Eval); fC20  = np.vectorize(FitC[1].Eval)
    fS30  = np.vectorize(FitS[2].Eval); fC30  = np.vectorize(FitC[2].Eval)
    fS40  = np.vectorize(FitS[3].Eval); fC40  = np.vectorize(FitC[3].Eval)
    fS60  = np.vectorize(FitS[4].Eval); fC60  = np.vectorize(FitC[4].Eval)
    fS80  = np.vectorize(FitS[5].Eval); fC80  = np.vectorize(FitC[5].Eval)
    fS100  = np.vectorize(FitS[6].Eval); fC100  = np.vectorize(FitC[6].Eval)
    fS120  = np.vectorize(FitS[7].Eval); fC120  = np.vectorize(FitC[7].Eval)

    binning_S = [0, 25, 53, 97]
    binning_C = [0, 25, 53, 97]


    DrawEnergyHist(dfs, energies, binning_S, binning_C, "totPMTSene", "totPMTCene")



    ROOT.gStyle.SetOptStat(0)
        
    for index, energy in enumerate(energies):
        print("Index ", index, "\tUsing file with energy: ", energy)

        data = dfs[index]

        # Define conditions
        #binning_S = [0, 15, 25, 35, 50, 75, 100, 120]
        #binning_C = [0, 15, 25, 35, 50, 75, 100, 120]
        binning_S = [0, 25, 53, 97]
        binning_C = [0, 25, 53, 97]

        #conditions_S = [
        #    (data["totPMTSene"] >= binning_S[0]) & (data["totPMTSene"] < binning_S[1]),
        #    (data["totPMTSene"] >= binning_S[1]) & (data["totPMTSene"] < binning_S[2]),
        #    (data["totPMTSene"] >= binning_S[2]) & (data["totPMTSene"] < binning_S[3]),
        #    (data["totPMTSene"] >= binning_S[3]) & (data["totPMTSene"] < binning_S[4]),
        #    (data["totPMTSene"] >= binning_S[4]) & (data["totPMTSene"] < binning_S[5]),
        #    (data["totPMTSene"] >= binning_S[5]) & (data["totPMTSene"] < binning_S[6]),
        #    (data["totPMTSene"] >= binning_S[6]) & (data["totPMTSene"] < binning_S[7]),
        #    (data["totPMTSene"] >= binning_S[7])
        #]
        # Define conditions (C channel)
        #conditions_C = [
        #    (data["totPMTCene"] >= binning_C[0]) & (data["totPMTCene"] < binning_C[1]),
        #    (data["totPMTCene"] >= binning_C[1]) & (data["totPMTCene"] < binning_C[2]),
        #    (data["totPMTCene"] >= binning_C[2]) & (data["totPMTCene"] < binning_C[3]),
        #    (data["totPMTCene"] >= binning_C[3]) & (data["totPMTCene"] < binning_C[4]),
        #    (data["totPMTCene"] >= binning_C[4]) & (data["totPMTCene"] < binning_C[5]),
        #    (data["totPMTCene"] >= binning_C[5]) & (data["totPMTCene"] < binning_C[6]),
        #    (data["totPMTCene"] >= binning_C[6]) & (data["totPMTCene"] < binning_C[7]),
        #    (data["totPMTCene"] >= binning_C[7])
        #]
        # Define the corresponding choices (functions to apply)
        #choices_S = [
        #    data["totPMTSene"] / fS20(data["AsymS"]),
        #    data["totPMTSene"] / fS20(data["AsymS"]),
        #    data["totPMTSene"] / fS20(data["AsymS"]),
        #    data["totPMTSene"] / fS20(data["AsymS"]),
        #    data["totPMTSene"] / fS20(data["AsymS"]),
        #    data["totPMTSene"] / fS20(data["AsymS"]),
        #    data["totPMTCene"] / fC20(data["AsymC"]),
        #    data["totPMTCene"] / fC20(data["AsymC"])
        #]
        #choices_C = [
        #    data["totPMTCene"] / fC20(data["AsymC"]),
        #    data["totPMTCene"] / fC20(data["AsymC"]),
        #    data["totPMTCene"] / fC20(data["AsymC"]),
        #    data["totPMTCene"] / fC20(data["AsymC"]),
        #    data["totPMTCene"] / fC20(data["AsymC"]),
        #    data["totPMTCene"] / fC20(data["AsymC"]),
        #    data["totPMTCene"] / fC20(data["AsymC"]),
        #    data["totPMTCene"] / fC20(data["AsymC"])
        #]


        conditions_S = [
            (data["totPMTSene"] >= binning_S[0]) & (data["totPMTSene"] < binning_S[1]),
            (data["totPMTSene"] >= binning_S[1]) & (data["totPMTSene"] < binning_S[2]),
            (data["totPMTSene"] >= binning_S[2]) & (data["totPMTSene"] < binning_S[3]),
            (data["totPMTSene"] >= binning_S[3])
        ]

        conditions_C = [
            (data["totPMTCene"] >= binning_C[0]) & (data["totPMTCene"] < binning_C[1]),
            (data["totPMTCene"] >= binning_C[1]) & (data["totPMTCene"] < binning_C[2]),
            (data["totPMTCene"] >= binning_C[2]) & (data["totPMTCene"] < binning_C[3]),
            (data["totPMTCene"] >= binning_C[3])
        ]

        # Define the corresponding choices (functions to apply)
        choices_S = [
            data["totPMTSene"] / fS20(data["AsymS"]),
            data["totPMTSene"] / fS40(data["AsymS"]),
            data["totPMTSene"] / fS80(data["AsymS"]),
            data["totPMTSene"] / fS120(data["AsymS"])
        ]
        # Define the corresponding choices (functions to apply)
        choices_C = [
            data["totPMTCene"] / fC20(data["AsymC"]),
            data["totPMTCene"] / fC40(data["AsymC"]),
            data["totPMTCene"] / fC80(data["AsymC"]),
            data["totPMTCene"] / fC120(data["AsymC"])
        ]


        # Use np.select to apply the conditions
        data["energyS"] = np.select(conditions_S, choices_S, default=data["totPMTSene"])
        data["energyC"] = np.select(conditions_C, choices_C, default=data["totPMTCene"])

        # Evaluate function on Asym variable and use it to correct energy
        #data["energyS"] = data["totPMTSene"]/fS20(data["AsymS"])
        #data["energyC"] = data["totPMTCene"]/fC20(data["AsymC"])



        # Fill histograms with corrected energy values
        HistScorrected = ROOT.TH1D("HistScorrected_{0}GeV".format(energy), "S Energy (corrected) {0}GeV; E [GeV]; Counts".format(energy), 100, energy-0.4*energy, energy+0.4*energy)
        HistSraw = ROOT.TH1D("HistSraw_{0}GeV".format(energy), "totPMTSene {0}GeV; E [GeV]; Counts".format(energy), 100, energy-0.4*energy, energy+0.4*energy)
        HistScorrected_Asymcut = ROOT.TH1D("HistScorrected_AsymCut_{0}GeV".format(energy), "S Energy (corrected, Asym cut) {0}GeV; E [GeV]; Counts".format(energy), 100, energy-0.4*energy, energy+0.4*energy)

        # Fill histograms with corrected energy values
        HistCcorrected = ROOT.TH1D("HistCcorrected_{0}GeV".format(energy), "C Energy (corrected) {0}GeV; E [GeV]; Counts".format(energy), 100, energy-0.4*energy, energy+0.4*energy)
        HistCraw = ROOT.TH1D("HistCraw_{0}GeV".format(energy), "totPMTCene {0}GeV; E [GeV]; Counts".format(energy), 100, energy-0.4*energy, energy+0.4*energy)
        HistCcorrected_Asymcut = ROOT.TH1D("HistCcorrected_AsymCut_{0}GeV".format(energy), "C Energy (corrected, Asym cut) {0}GeV; E [GeV]; Counts".format(energy), 100, energy-0.4*energy, energy+0.4*energy)

        # Fill histograms with corrected energy values
        HistCombCorrected = ROOT.TH1D("HistCombcorrected_{0}GeV".format(energy), "Combined Energy (corrected) {0}GeV; E [GeV]; Counts".format(energy), 100, energy-0.4*energy, energy+0.4*energy)


        for eneS, pmtS, eneC, pmtC in zip(data["energyS"].values, data["totPMTSene"].values, data["energyC"].values, data["totPMTCene"].values):
            HistScorrected.Fill(eneS)
            HistSraw.Fill(pmtS)
            HistCcorrected.Fill(eneC)
            HistCraw.Fill(pmtC)


        dfCorrected_array.append(data)
        # filter tree using only events with an asymmetry within range
        AsymCut = 0.7
        #data = data[ (np.abs(data["AsymS"]<AsymCut) ) & (np.abs(data["AsymC"]<AsymCut) ) ]
        #data_filtered = data[ (np.abs(data["AsymS"])<AsymCut ) & (np.abs(data["AsymC"])<AsymCut ) & (np.abs(data["BaryS"])<4) & (np.abs(data["BaryC"])<4) ]
        data_filtered = data[ (np.abs(data["AsymS"])<AsymCut ) & (np.abs(data["AsymC"])<AsymCut ) ]

        for eneS, eneC in zip(data_filtered["energyS"].values, data_filtered["energyC"].values):
            HistScorrected_Asymcut.Fill(eneS)
            HistCcorrected_Asymcut.Fill(eneC)
            HistCombCorrected.Fill( (eneS+eneC)/2 )




        '''
        binMax = 0
        if(varProf == "YDWC2"): binMax = 20
        elif(varProf == "XDWC2"): binMax = 30 
        # Test correction: plot totPMTSene and energyS over YDWC2
        SciEneVsYDWCprof = ROOT.TProfile("SciEneVsYDWCprof{0}".format(energy), "Corrected Energy (S) vs YDWC2 ({0}GeV); {1} [mm]; E [GeV]".format(energy, varProf), 30, -30, binMax)
        SpmtVsYDWCprof = ROOT.TProfile("SpmtVsYDWCprof{0}".format(energy), "totPMTSene vs YDWC2 ({0}GeV); {1} [mm]; E [GeV]".format(energy, varProf), 30, -30, binMax)
        CerEneVsYDWCprof = ROOT.TProfile("EneVsYDWCprof{0}".format(energy), "Corrected Energy (C) vs YDWC2 ({0}GeV); {1} [mm]; E [GeV]".format(energy, varProf), 30, -30, binMax)
        CpmtVsYDWCprof = ROOT.TProfile("CpmtVsYDWCprof{0}".format(energy), "totPMTCene vs YDWC2 ({0}GeV); {1} [mm]; E [GeV]".format(energy, varProf), 30, -30, binMax)



        for pmtS, eneS, pmtC, eneC, ydwc2 in zip(data["totPMTSene"], data["energyS"], data["totPMTCene"], data["energyC"], data["YDWC2"]):
            SciEneVsYDWCprof.Fill(ydwc2, eneS)
            SpmtVsYDWCprof.Fill(ydwc2, pmtS)
            CerEneVsYDWCprof.Fill(ydwc2, eneC)
            CpmtVsYDWCprof.Fill(ydwc2, pmtC)


        SciEneVsYDWCprof.SetLineColor(ROOT.kRed); SciEneVsYDWCprof.SetLineWidth(2); SciEneVsYDWCprof.SetMarkerColor(ROOT.kRed); SciEneVsYDWCprof.SetMarkerStyle(55)
        SpmtVsYDWCprof.SetLineColor(ROOT.kBlue); SpmtVsYDWCprof.SetLineWidth(2); SpmtVsYDWCprof.SetMarkerColor(ROOT.kBlue); SpmtVsYDWCprof.SetMarkerStyle(73)
        CerEneVsYDWCprof.SetLineColor(ROOT.kRed); CerEneVsYDWCprof.SetLineWidth(2); CerEneVsYDWCprof.SetMarkerColor(ROOT.kRed); CerEneVsYDWCprof.SetMarkerStyle(55)
        CpmtVsYDWCprof.SetLineColor(ROOT.kBlue); CpmtVsYDWCprof.SetLineWidth(2); CpmtVsYDWCprof.SetMarkerColor(ROOT.kBlue); CpmtVsYDWCprof.SetMarkerStyle(73)


        ROOT.gStyle.SetOptStat(0)
        ctest = ROOT.TCanvas("cEneProfAsymmetry{0}".format(energy), "cEneProfAsymmetry{0}".format(energy), 1400, 1200)
        ctest.SetLeftMargin(0.15)
        eneSprof.SetLineWidth(2); eneSprof.SetLineColor(ROOT.kRed); eneSprof.SetMarkerStyle(ROOT.kFullCircle); eneSprof.SetMarkerColor(ROOT.kRed)
        eneSprof.SetMinimum(0.8); eneSprof.SetMaximum(1.2); eneSprof.Draw()
        #ctestS.SaveAs("testS{0}GeV.png".format(energy))

        #ctestC = ROOT.TCanvas("c{0}".format(energy), "c{0}".format(energy), 1400, 1200)
        eneCprof.SetLineWidth(2); eneCprof.SetLineColor(ROOT.kBlue); eneCprof.SetMarkerStyle(ROOT.kFullCircle); eneCprof.SetMarkerColor(ROOT.kBlue)
        eneCprof.SetMinimum(0.8); eneCprof.SetMaximum(1.2); eneCprof.Draw("same")
        ctest.SaveAs("EnergyProfOverAsymmetry_{0}GeV.png".format(energy))
        


        cprofS = ROOT.TCanvas("cprofS{0}".format(energy), "cprofS{0}".format(energy), 1400, 1200)
        SciEneVsYDWCprof.SetMaximum(energy*1.2); SciEneVsYDWCprof.SetMinimum(energy*0.8)
        SciEneVsYDWCprof.Draw()
        SpmtVsYDWCprof.Draw("same")
        leg = ROOT.TLegend(0.75, 0.82, 0.95, 0.93)
        leg.AddEntry(SciEneVsYDWCprof, "S energy (corrected)", "PL"); leg.AddEntry(SpmtVsYDWCprof, "totPMTSene (Raw)", "PL")
        leg.SetTextSize(0.016)
        leg.Draw()
        cprofS.SaveAs("SciEneProfY{0}.png".format(energy))

        cprofC = ROOT.TCanvas("cprofC{0}".format(energy), "cprofC{0}".format(energy), 1400, 1200)
        CerEneVsYDWCprof.SetMaximum(energy*1.2); CerEneVsYDWCprof.SetMinimum(energy*0.8)
        CerEneVsYDWCprof.Draw()
        CpmtVsYDWCprof.Draw("same")
        leg = ROOT.TLegend(0.75, 0.82, 0.95, 0.93)
        leg.AddEntry(CerEneVsYDWCprof, "C energy (corrected)", "PL"); leg.AddEntry(CpmtVsYDWCprof, "totPMTCene (Raw)", "PL")
        leg.SetTextSize(0.016)
        leg.Draw()
        cprofC.SaveAs("EneProfY{0}.png".format(energy))
        
        



        # Profiles over barycenter variable
        BarycenterEneProfS = ROOT.TProfile("BarycenterEneProfS{0}".format(energy), "S Energy over Barycenter position ({0} GeV); Barycenter Y [mm]; E [GeV]".format(energy), 100, -25, 20)
        BarycenterEneProfC = ROOT.TProfile("BarycenterEneProfC{0}".format(energy), "C Energy over Barycenter position ({0} GeV); Barycenter Y [mm]; E [GeV]".format(energy), 100, -25, 20)

        # Barycenter profile over YDWC position
        YDWCBarycenterProfS = ROOT.TProfile("YDWCBarycenterProfS{0}".format(energy), "S Barycenter position over Y coordinate ({0} GeV); YDWC2 [mm]; Barycenter Y [mm]".format(energy), 100, -25, 20)
        YDWCBarycenterProfC = ROOT.TProfile("YDWCBarycenterProfC{0}".format(energy), "C Barycenter position over Y coordinate ({0} GeV); YDWC2 [mm]; Barycenter Y [mm]".format(energy), 100, -25, 20)

        # Histtograms with energy (raw and corrected) Vs Asymmetry and Barycenter position
        SpmtAsymBarHist = ROOT.TH2D("SpmtAsymBarHist_{0}".format(energy), "totPMTSene Vs Asymmetry Vs BarycenterY ({0} GeV); Asymmetry; Y Barycenter [mm]".format(energy), 50, -1, 1, 50, -25, 20)
        CpmtAsymBarHist = ROOT.TH2D("CpmtAsymBarHist_{0}".format(energy), "totPMTCene Vs Asymmetry Vs BarycenterY ({0} GeV); Asymmetry; Y Barycenter [mm]".format(energy), 50, -1, 1, 50, -25, 20)



        for pmtS, pmtC, ydwc2, barS, barC in zip(data["totPMTSene"], data["totPMTCene"], data["YDWC2"], data["BaryS"], data["BaryC"]):
            BarycenterEneProfS.Fill(barS, pmtS)
            YDWCBarycenterProfS.Fill(ydwc2, barS)
            BarycenterEneProfC.Fill(barC, pmtC)
            YDWCBarycenterProfC.Fill(ydwc2, barC)

        for pmtS, barS, asymS, pmtC, barC, asymC in zip(data["totPMTSene"], data["BaryS"], data["AsymS"], data["totPMTCene"], data["BaryC"], data["AsymC"]):    
            SpmtAsymBarHist.Fill(asymS, barS, pmtS)
            CpmtAsymBarHist.Fill(asymC, barC, pmtC)


        
        cEneBarProfS = ROOT.TCanvas("cEneBarProfS{0}".format(energy), "cEneBarProfS{0}".format(energy), 1400, 1200)
        BarycenterEneProfS.SetLineColor(ROOT.kRed); BarycenterEneProfS.SetMarkerStyle(23); BarycenterEneProfS.SetLineWidth(2); BarycenterEneProfS.SetMarkerSize(2); BarycenterEneProfS.SetMarkerColor(ROOT.kRed)
        BarycenterEneProfC.SetLineColor(ROOT.kBlue); BarycenterEneProfC.SetMarkerStyle(23); BarycenterEneProfC.SetLineWidth(2); BarycenterEneProfC.SetMarkerSize(2); BarycenterEneProfC.SetMarkerColor(ROOT.kBlue)
        BarycenterEneProfS.SetTitle("Energy over Barycenter Y position ({0})".format(energy))
        BarycenterEneProfS.Draw("P")
        BarycenterEneProfC.Draw("P same")
        leg = ROOT.TLegend(0.78, 0.85, 0.95, 0.93)
        leg.AddEntry(BarycenterEneProfS, "S Channel", "PL"); leg.AddEntry(BarycenterEneProfC, "C Channel", "PL")
        leg.SetTextSize(0.016)
        leg.Draw()
        cEneBarProfS.SaveAs("EneBarycenterProf{0}.png".format(energy))


        cYdwcBarProfS = ROOT.TCanvas("cYdwcBarProf{0}".format(energy), "cYdwcBarProf{0}".format(energy), 1400, 1200)
        YDWCBarycenterProfS.SetLineColor(ROOT.kRed); YDWCBarycenterProfS.SetMarkerStyle(23); YDWCBarycenterProfS.SetLineWidth(2); YDWCBarycenterProfS.SetMarkerSize(2); YDWCBarycenterProfS.SetMarkerColor(ROOT.kRed)
        YDWCBarycenterProfC.SetLineColor(ROOT.kBlue); YDWCBarycenterProfC.SetMarkerStyle(23); YDWCBarycenterProfC.SetLineWidth(2); YDWCBarycenterProfC.SetMarkerSize(2); YDWCBarycenterProfC.SetMarkerColor(ROOT.kBlue)
        YDWCBarycenterProfS.Draw("P")
        YDWCBarycenterProfC.Draw("P same")
        leg = ROOT.TLegend(0.8, 0.85, 0.995, 0.93)
        leg.AddEntry(BarycenterEneProfS, "S Channel", "PL"); leg.AddEntry(BarycenterEneProfC, "C Channel", "PL")
        leg.SetTextSize(0.016)
        leg.Draw()
        cYdwcBarProfS.SaveAs("YdwcBarycenterProf{0}.png".format(energy))


        cSpmtAsymBarHist = ROOT.TCanvas("SpmtAsymBarHist_{0}".format(energy), "SpmtAsymBarHist_{0}".format(energy), 1400, 1200)
        SpmtAsymBarHist.Draw("colz")
        cSpmtAsymBarHist.SaveAs("SpmtAsymBarHist_{0}GeV.png".format(energy))
        
        cCpmtAsymBarHist = ROOT.TCanvas("CpmtAsymBarHist_{0}".format(energy), "CpmtAsymBarHist_{0}".format(energy), 1400, 1200)
        CpmtAsymBarHist.Draw("colz")
        cCpmtAsymBarHist.SaveAs("CpmtAsymBarHist_{0}GeV.png".format(energy))
        '''






        cS = ROOT.TCanvas("cS{0}".format(energy), "cS{0}".format(energy), 1400, 1200)
        cS.SetLeftMargin(0.15)
        cS.SetRightMargin(0.07)
        HistSraw.SetLineWidth(2); HistSraw.SetLineColor(ROOT.kBlue+1); HistSraw.SetFillColorAlpha(ROOT.kAzure+10, 0.08)
        HistScorrected.SetLineWidth(2); HistScorrected.SetLineColor(ROOT.kRed)
        HistScorrected_Asymcut.SetLineWidth(2); HistScorrected_Asymcut.SetLineColor(ROOT.kGreen+2)
        max_value = max(HistScorrected.GetMaximum(), HistSraw.GetMaximum(), HistScorrected_Asymcut.GetMaximum())*1.15
        HistSraw.SetTitle("Reco S energy at {0}GeV".format(energy))
        HistSraw.SetMaximum(max_value)
        HistSraw.Draw()
        HistScorrected.Draw("same")
        HistScorrected_Asymcut.Draw("same")
        leg = ROOT.TLegend(0.62, 0.75, 0.92, 0.88)
        leg.SetTextSize(0.018)
        leg.AddEntry(HistSraw, "totPMTSene (raw)", "L")
        leg.AddEntry(HistScorrected, "S energy (corrected)", "L")
        leg.AddEntry(HistScorrected_Asymcut, "S energy (corrected, |asym|<{0})".format(AsymCut), "L")
        leg.Draw()
        cS.SaveAs("SciEnergy{0}GeV.png".format(energy))
        
        

        cC = ROOT.TCanvas("cC{0}".format(energy), "cC{0}".format(energy), 1400, 1200)
        cC.SetLeftMargin(0.15)
        cC.SetRightMargin(0.07)
        HistCraw.SetLineWidth(2); HistCraw.SetLineColor(ROOT.kBlue+1); HistCraw.SetFillColorAlpha(ROOT.kAzure+10, 0.08)
        HistCcorrected.SetLineWidth(2); HistCcorrected.SetLineColor(ROOT.kRed)
        HistCcorrected_Asymcut.SetLineWidth(2); HistCcorrected_Asymcut.SetLineColor(ROOT.kGreen+2)
        max_value = max(HistCcorrected.GetMaximum(), HistCraw.GetMaximum(), HistCcorrected_Asymcut.GetMaximum())*1.15
        HistCraw.SetTitle("Reco C energy at {0}GeV".format(energy))
        HistCraw.SetMaximum(max_value)
        HistCraw.Draw()
        HistCcorrected.Draw("same")
        HistCcorrected_Asymcut.Draw("same")
        leg = ROOT.TLegend(0.62, 0.75, 0.92, 0.88)
        leg.SetTextSize(0.018)
        leg.AddEntry(HistCraw, "totPMTCene (raw)", "L")
        leg.AddEntry(HistCcorrected, "C energy (corrected)", "L")
        leg.AddEntry(HistCcorrected_Asymcut, "C energy (corrected, |asym|<{0})".format(AsymCut), "L")
        leg.Draw()
        cC.SaveAs("CerEnergy{0}GeV.png".format(energy))


        
            

        HistScorrected_Asymcut.Fit("gaus", "Q")
        fitS = HistScorrected_Asymcut.GetFunction("gaus")

        HistCcorrected_Asymcut.Fit("gaus", "Q")
        fitC = HistCcorrected_Asymcut.GetFunction("gaus")

        HistCombCorrected.Fit("gaus","Q")
        fitComb = HistCombCorrected.GetFunction("gaus")


        ROOT.gStyle.SetOptFit(111)
        # refit within -1.5 and +3 sigma
        HistScorrected_Asymcut.Fit("gaus", "Q", "", fitS.GetParameter(1)-1.5*fitS.GetParameter(2), fitS.GetParameter(1)+3*fitS.GetParameter(2)) 
        HistCcorrected_Asymcut.Fit("gaus", "Q", "", fitC.GetParameter(1)-1.5*fitC.GetParameter(2), fitC.GetParameter(1)+3*fitC.GetParameter(2)) 
        HistCombCorrected.Fit("gaus", "Q", "", fitComb.GetParameter(1)-1.5*fitComb.GetParameter(2), fitComb.GetParameter(1)+3*fitComb.GetParameter(2)) 


        scihist = ROOT.TCanvas("cSciEnergy{0}".format(energy), "cSciEnergy{0}".format(energy), 1400, 1200)
        HistScorrected_Asymcut.Draw()
        scihist.SaveAs("Scienehist{0}.png".format(energy))
        cerhist = ROOT.TCanvas("cCerEnergy{0}".format(energy), "cCerEnergy{0}".format(energy), 1400, 1200)
        HistCcorrected_Asymcut.Draw()
        cerhist.SaveAs("Cerenehist{0}.png".format(energy))
        combhist = ROOT.TCanvas("cCombEnergy{0}".format(energy), "cCombEnergy{0}".format(energy), 1400, 1200)
        HistCombCorrected.Draw()
        combhist.SaveAs("Combenehist{0}.png".format(energy))


        BestFitS = HistScorrected_Asymcut.GetFunction("gaus")
        BestFitC = HistCcorrected_Asymcut.GetFunction("gaus")
        BestFitComb = HistCombCorrected.GetFunction("gaus")

        MeanVec_S.append(BestFitS.GetParameter(1)); MeanErrVec_S.append(BestFitS.GetParError(1)); RmsVec_S.append(BestFitS.GetParameter(2)); RmsErrVec_S.append(BestFitS.GetParError(2))
        MeanVec_C.append(BestFitC.GetParameter(1)); MeanErrVec_C.append(BestFitC.GetParError(1)); RmsVec_C.append(BestFitC.GetParameter(2)); RmsErrVec_C.append(BestFitC.GetParError(2))
        MeanVec_Comb.append(BestFitComb.GetParameter(1)); MeanErrVec_Comb.append(BestFitComb.GetParError(1)); RmsVec_Comb.append(BestFitComb.GetParameter(2)); RmsErrVec_Comb.append(BestFitComb.GetParError(2))
        


    DrawEnergyHist(dfCorrected_array, energies, binning_S, binning_C, "energyS", "energyC")


    
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

            


        



    '''


    # End loop here, store TF1 function to be used later

    # vectorize it
    funcSvector = np.vectorize(funcS.Eval)
    funcCvector = np.vectorize(funcC.Eval)

    # Evaluate function on Asym variable and use it to correct energy
    data["energyS"] = data["totPMTSene"]/funcSvector(data["AsymS"])
    data["energyC"] = data["totPMTCene"]/funcCvector(data["AsymC"])


    binMax = 0
    if(varProf == "YDWC2"): binMax = 20
    elif(varProf == "XDWC2"): binMax = 30 
    # Test correction: plot totPMTSene and energyS over YDWC2
    SciEneVsYDWCprof = ROOT.TProfile("SciEneVsYDWCprof{0}".format(energy), "Corrected Energy (S) vs YDWC2 ({0}GeV); {1} [mm]; E [GeV]".format(energy, varProf), 30, -30, binMax)
    SpmtVsYDWCprof = ROOT.TProfile("SpmtVsYDWCprof{0}".format(energy), "totPMTSene vs YDWC2 ({0}GeV); {1} [mm]; E [GeV]".format(energy, varProf), 30, -30, binMax)
    CerEneVsYDWCprof = ROOT.TProfile("EneVsYDWCprof{0}".format(energy), "Corrected Energy (C) vs YDWC2 ({0}GeV); {1} [mm]; E [GeV]".format(energy, varProf), 30, -30, binMax)
    CpmtVsYDWCprof = ROOT.TProfile("CpmtVsYDWCprof{0}".format(energy), "totPMTCene vs YDWC2 ({0}GeV); {1} [mm]; E [GeV]".format(energy, varProf), 30, -30, binMax)



    for pmtS, eneS, pmtC, eneC, ydwc2 in zip(data["totPMTSene"], data["energyS"], data["totPMTCene"], data["energyC"], data["YDWC2"]):
        SciEneVsYDWCprof.Fill(ydwc2, eneS)
        SpmtVsYDWCprof.Fill(ydwc2, pmtS)
        CerEneVsYDWCprof.Fill(ydwc2, eneC)
        CpmtVsYDWCprof.Fill(ydwc2, pmtC)


    SciEneVsYDWCprof.SetLineColor(ROOT.kRed); SciEneVsYDWCprof.SetLineWidth(2); SciEneVsYDWCprof.SetMarkerColor(ROOT.kRed); SciEneVsYDWCprof.SetMarkerStyle(55)
    SpmtVsYDWCprof.SetLineColor(ROOT.kBlue); SpmtVsYDWCprof.SetLineWidth(2); SpmtVsYDWCprof.SetMarkerColor(ROOT.kBlue); SpmtVsYDWCprof.SetMarkerStyle(73)
    CerEneVsYDWCprof.SetLineColor(ROOT.kRed); CerEneVsYDWCprof.SetLineWidth(2); CerEneVsYDWCprof.SetMarkerColor(ROOT.kRed); CerEneVsYDWCprof.SetMarkerStyle(55)
    CpmtVsYDWCprof.SetLineColor(ROOT.kBlue); CpmtVsYDWCprof.SetLineWidth(2); CpmtVsYDWCprof.SetMarkerColor(ROOT.kBlue); CpmtVsYDWCprof.SetMarkerStyle(73)


    ''''''
    ROOT.gStyle.SetOptStat(0)
    ctest = ROOT.TCanvas("cEneProfAsymmetry{0}".format(energy), "cEneProfAsymmetry{0}".format(energy), 1400, 1200)
    ctest.SetLeftMargin(0.15)
    eneSprof.SetLineWidth(2); eneSprof.SetLineColor(ROOT.kRed); eneSprof.SetMarkerStyle(ROOT.kFullCircle); eneSprof.SetMarkerColor(ROOT.kRed)
    eneSprof.SetMinimum(0.8); eneSprof.SetMaximum(1.2); eneSprof.Draw()
    #ctestS.SaveAs("testS{0}GeV.png".format(energy))

    #ctestC = ROOT.TCanvas("c{0}".format(energy), "c{0}".format(energy), 1400, 1200)
    eneCprof.SetLineWidth(2); eneCprof.SetLineColor(ROOT.kBlue); eneCprof.SetMarkerStyle(ROOT.kFullCircle); eneCprof.SetMarkerColor(ROOT.kBlue)
    eneCprof.SetMinimum(0.8); eneCprof.SetMaximum(1.2); eneCprof.Draw("same")
    ctest.SaveAs("EnergyProfOverAsymmetry_{0}GeV.png".format(energy))
    


    cprofS = ROOT.TCanvas("cprofS{0}".format(energy), "cprofS{0}".format(energy), 1400, 1200)
    SciEneVsYDWCprof.SetMaximum(energy*1.2); SciEneVsYDWCprof.SetMinimum(energy*0.8)
    SciEneVsYDWCprof.Draw()
    SpmtVsYDWCprof.Draw("same")
    leg = ROOT.TLegend(0.75, 0.82, 0.95, 0.93)
    leg.AddEntry(SciEneVsYDWCprof, "S energy (corrected)", "PL"); leg.AddEntry(SpmtVsYDWCprof, "totPMTSene (Raw)", "PL")
    leg.SetTextSize(0.016)
    leg.Draw()
    cprofS.SaveAs("SciEneProfY{0}.png".format(energy))

    cprofC = ROOT.TCanvas("cprofC{0}".format(energy), "cprofC{0}".format(energy), 1400, 1200)
    CerEneVsYDWCprof.SetMaximum(energy*1.2); CerEneVsYDWCprof.SetMinimum(energy*0.8)
    CerEneVsYDWCprof.Draw()
    CpmtVsYDWCprof.Draw("same")
    leg = ROOT.TLegend(0.75, 0.82, 0.95, 0.93)
    leg.AddEntry(CerEneVsYDWCprof, "C energy (corrected)", "PL"); leg.AddEntry(CpmtVsYDWCprof, "totPMTCene (Raw)", "PL")
    leg.SetTextSize(0.016)
    leg.Draw()
    cprofC.SaveAs("EneProfY{0}.png".format(energy))
    



    # Profiles over barycenter variable
    BarycenterEneProfS = ROOT.TProfile("BarycenterEneProfS{0}".format(energy), "S Energy over Barycenter position ({0} GeV); Barycenter Y [mm]; E [GeV]".format(energy), 100, -25, 20)
    BarycenterEneProfC = ROOT.TProfile("BarycenterEneProfC{0}".format(energy), "C Energy over Barycenter position ({0} GeV); Barycenter Y [mm]; E [GeV]".format(energy), 100, -25, 20)

    # Barycenter profile over YDWC position
    YDWCBarycenterProfS = ROOT.TProfile("YDWCBarycenterProfS{0}".format(energy), "S Barycenter position over Y coordinate ({0} GeV); YDWC2 [mm]; Barycenter Y [mm]".format(energy), 100, -25, 20)
    YDWCBarycenterProfC = ROOT.TProfile("YDWCBarycenterProfC{0}".format(energy), "C Barycenter position over Y coordinate ({0} GeV); YDWC2 [mm]; Barycenter Y [mm]".format(energy), 100, -25, 20)

    # Histtograms with energy (raw and corrected) Vs Asymmetry and Barycenter position
    SpmtAsymBarHist = ROOT.TH2D("SpmtAsymBarHist_{0}".format(energy), "totPMTSene Vs Asymmetry Vs BarycenterY ({0} GeV); Asymmetry; Y Barycenter [mm]".format(energy), 50, -1, 1, 50, -25, 20)
    CpmtAsymBarHist = ROOT.TH2D("CpmtAsymBarHist_{0}".format(energy), "totPMTCene Vs Asymmetry Vs BarycenterY ({0} GeV); Asymmetry; Y Barycenter [mm]".format(energy), 50, -1, 1, 50, -25, 20)



    for pmtS, pmtC, ydwc2, barS, barC in zip(data["totPMTSene"], data["totPMTCene"], data["YDWC2"], data["BaryS"], data["BaryC"]):
        BarycenterEneProfS.Fill(barS, pmtS)
        YDWCBarycenterProfS.Fill(ydwc2, barS)
        BarycenterEneProfC.Fill(barC, pmtC)
        YDWCBarycenterProfC.Fill(ydwc2, barC)

    for pmtS, barS, asymS, pmtC, barC, asymC in zip(data["totPMTSene"], data["BaryS"], data["AsymS"], data["totPMTCene"], data["BaryC"], data["AsymC"]):    
        SpmtAsymBarHist.Fill(asymS, barS, pmtS)
        CpmtAsymBarHist.Fill(asymC, barC, pmtC)


    
    cEneBarProfS = ROOT.TCanvas("cEneBarProfS{0}".format(energy), "cEneBarProfS{0}".format(energy), 1400, 1200)
    BarycenterEneProfS.SetLineColor(ROOT.kRed); BarycenterEneProfS.SetMarkerStyle(23); BarycenterEneProfS.SetLineWidth(2); BarycenterEneProfS.SetMarkerSize(2); BarycenterEneProfS.SetMarkerColor(ROOT.kRed)
    BarycenterEneProfC.SetLineColor(ROOT.kBlue); BarycenterEneProfC.SetMarkerStyle(23); BarycenterEneProfC.SetLineWidth(2); BarycenterEneProfC.SetMarkerSize(2); BarycenterEneProfC.SetMarkerColor(ROOT.kBlue)
    BarycenterEneProfS.SetTitle("Energy over Barycenter Y position ({0})".format(energy))
    BarycenterEneProfS.Draw("P")
    BarycenterEneProfC.Draw("P same")
    leg = ROOT.TLegend(0.78, 0.85, 0.95, 0.93)
    leg.AddEntry(BarycenterEneProfS, "S Channel", "PL"); leg.AddEntry(BarycenterEneProfC, "C Channel", "PL")
    leg.SetTextSize(0.016)
    leg.Draw()
    cEneBarProfS.SaveAs("EneBarycenterProf{0}.png".format(energy))


    cYdwcBarProfS = ROOT.TCanvas("cYdwcBarProf{0}".format(energy), "cYdwcBarProf{0}".format(energy), 1400, 1200)
    YDWCBarycenterProfS.SetLineColor(ROOT.kRed); YDWCBarycenterProfS.SetMarkerStyle(23); YDWCBarycenterProfS.SetLineWidth(2); YDWCBarycenterProfS.SetMarkerSize(2); YDWCBarycenterProfS.SetMarkerColor(ROOT.kRed)
    YDWCBarycenterProfC.SetLineColor(ROOT.kBlue); YDWCBarycenterProfC.SetMarkerStyle(23); YDWCBarycenterProfC.SetLineWidth(2); YDWCBarycenterProfC.SetMarkerSize(2); YDWCBarycenterProfC.SetMarkerColor(ROOT.kBlue)
    YDWCBarycenterProfS.Draw("P")
    YDWCBarycenterProfC.Draw("P same")
    leg = ROOT.TLegend(0.8, 0.85, 0.995, 0.93)
    leg.AddEntry(BarycenterEneProfS, "S Channel", "PL"); leg.AddEntry(BarycenterEneProfC, "C Channel", "PL")
    leg.SetTextSize(0.016)
    leg.Draw()
    cYdwcBarProfS.SaveAs("YdwcBarycenterProf{0}.png".format(energy))


    cSpmtAsymBarHist = ROOT.TCanvas("SpmtAsymBarHist_{0}".format(energy), "SpmtAsymBarHist_{0}".format(energy), 1400, 1200)
    SpmtAsymBarHist.Draw("colz")
    cSpmtAsymBarHist.SaveAs("SpmtAsymBarHist_{0}GeV.png".format(energy))
    
    cCpmtAsymBarHist = ROOT.TCanvas("CpmtAsymBarHist_{0}".format(energy), "CpmtAsymBarHist_{0}".format(energy), 1400, 1200)
    CpmtAsymBarHist.Draw("colz")
    cCpmtAsymBarHist.SaveAs("CpmtAsymBarHist_{0}GeV.png".format(energy))
    


    ################### Fill energy histograms ######################
    ROOT.gStyle.SetOptStat(1)

    # Fill histograms with corrected energy values
    EneHistS = ROOT.TH1D("EneHistS_{0}GeV".format(energy), "S Energy (corrected) {0}GeV; E [GeV]; Counts".format(energy), 100, energy-0.4*energy, energy+0.4*energy)
    EneHistC = ROOT.TH1D("EneHistC_{0}GeV".format(energy), "C Energy (corrected) {0}GeV; E [GeV]; Counts".format(energy), 100, energy-0.4*energy, energy+0.4*energy)
    EneHistComb = ROOT.TH1D("EneHistComb_{0}GeV".format(energy), "(S+C)/2 {0}GeV; E [GeV]; Counts".format(energy), 100, energy-0.4*energy, energy+0.4*energy)


    # filter tree using only events with an asymmetry within range
    AsymCut = 0.5
    print(data.shape)
    #data = data[ (np.abs(data["AsymS"]<AsymCut) ) & (np.abs(data["AsymC"]<AsymCut) ) ]
    data = data[ (np.abs(data["AsymS"])<AsymCut ) & (np.abs(data["AsymC"])<AsymCut ) & (np.abs(data["BaryS"])<4) & (np.abs(data["BaryC"])<4) ]

    print("\n After Asym Cut: ")
    print(data.shape)


    for eneS, eneC in zip(data["energyS"].values, data["energyC"].values):
        EneHistS.Fill(eneS)
        EneHistC.Fill(eneC)
        EneHistComb.Fill( (eneS+eneC)/2 )

    EneHistS.Fit("gaus", "Q")
    fitS = EneHistS.GetFunction("gaus")

    EneHistC.Fit("gaus", "Q")
    fitC = EneHistC.GetFunction("gaus")

    EneHistComb.Fit("gaus","Q")
    fitComb = EneHistComb.GetFunction("gaus")


    ROOT.gStyle.SetOptFit(111)
    # refit within -1.5 and +3 sigma
    EneHistS.Fit("gaus", "Q", "", fitS.GetParameter(1)-1.5*fitS.GetParameter(2), fitS.GetParameter(1)+3*fitS.GetParameter(2)) 
    EneHistC.Fit("gaus", "Q", "", fitC.GetParameter(1)-1.5*fitC.GetParameter(2), fitC.GetParameter(1)+3*fitC.GetParameter(2)) 
    EneHistComb.Fit("gaus", "Q", "", fitComb.GetParameter(1)-1.5*fitComb.GetParameter(2), fitComb.GetParameter(1)+3*fitComb.GetParameter(2)) 


    scihist = ROOT.TCanvas("cSciEnergy{0}".format(energy), "cSciEnergy{0}".format(energy), 1400, 1200)
    EneHistS.Draw()
    scihist.SaveAs("Scienehist{0}.png".format(energy))

    cerhist = ROOT.TCanvas("cCerEnergy{0}".format(energy), "cCerEnergy{0}".format(energy), 1400, 1200)
    EneHistC.Draw()
    cerhist.SaveAs("Cerenehist{0}.png".format(energy))

    combhist = ROOT.TCanvas("cCombEnergy{0}".format(energy), "cCombEnergy{0}".format(energy), 1400, 1200)
    EneHistComb.Draw()
    combhist.SaveAs("Combenehist{0}.png".format(energy))


    BestFitS = EneHistS.GetFunction("gaus")
    BestFitC = EneHistC.GetFunction("gaus")
    BestFitComb = EneHistComb.GetFunction("gaus")

    MeanVec_S.append(BestFitS.GetParameter(1)); MeanErrVec_S.append(BestFitS.GetParError(1)); RmsVec_S.append(BestFitS.GetParameter(2)); RmsErrVec_S.append(BestFitS.GetParError(2))
    MeanVec_C.append(BestFitC.GetParameter(1)); MeanErrVec_C.append(BestFitC.GetParError(1)); RmsVec_C.append(BestFitC.GetParameter(2)); RmsErrVec_C.append(BestFitC.GetParError(2))
    MeanVec_Comb.append(BestFitComb.GetParameter(1)); MeanErrVec_Comb.append(BestFitComb.GetParError(1)); RmsVec_Comb.append(BestFitComb.GetParameter(2)); RmsErrVec_Comb.append(BestFitComb.GetParError(2))
    '''



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
    '''
    

    #df.to_csv("DRcaloData.csv", index=False, sep="\t")

            

    







if __name__ == "__main__":
    main()