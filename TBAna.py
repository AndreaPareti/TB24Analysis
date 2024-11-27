import pandas as pd
import numpy as np
import ROOT 
import uproot




def main():
    print("Hello there")
    runs = ["0786", "0766", "0772", "0774", "0775", "0778", "0779", "0792"]
    energies = [10, 20, 30, 40, 60, 80, 100, 120]
    myCut = "(abs(XDWC2 - XDWC1) < 5) & (abs(YDWC2 - YDWC1)<5) & (MCounter<200) & (TailC<300) & (C2>160) & ( (totLeakage - L20)<5000 ) & (PShower>550)"
    varProf = "YDWC2"

    MeanVec_S=[]; RmsVec_S=[]; MeanErrVec_S=[]; RmsErrVec_S=[] 
    MeanVec_C=[]; RmsVec_C=[]; MeanErrVec_C=[]; RmsErrVec_C=[] 
    MeanVec_Comb=[]; RmsVec_Comb=[]; MeanErrVec_Comb=[]; RmsErrVec_Comb=[] 


    infolder = "/home/storage/data_apareti/TB24/EnergyScan/"
    treename = "Ftree"

    for run, energy in zip(runs, energies):
        filename = "physics_sps2024_run" + run + ".root"
        print(filename, energy)
        ROOT.gStyle.SetOptFit(0)

        root_file = uproot.open(infolder+filename)
        tree = root_file[treename]

        data = tree.arrays(cut=myCut, library="pd")

        # define Asymmetry variable 
        data["AsymS"] = (data["TS11"] - data["TS15"]) / (data["TS11"] + data["TS15"] )
        data["AsymC"] = (data["TC11"] - data["TC15"]) / (data["TC11"] + data["TC15"] )

        # Profile normalised energy Vs Asymmetry
        eneSprof = ROOT.TProfile("eneSprof_{0}GeV".format(energy), "eneSprof_{0}GeV".format(energy), 100, -1, 1)
        eneCprof = ROOT.TProfile("eneCprof_{0}GeV".format(energy), "eneCprof_{0}GeV".format(energy), 100, -1, 1)

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



        ctestS = ROOT.TCanvas("cS{0}".format(energy), "cS{0}".format(energy), 1400, 1200)
        eneSprof.SetMinimum(0.8); eneSprof.SetMaximum(1.2); eneSprof.Draw()
        ctestS.SaveAs("testS{0}GeV.png".format(energy))

        ctestC = ROOT.TCanvas("c{0}".format(energy), "c{0}".format(energy), 1400, 1200)
        eneCprof.SetMinimum(0.8); eneCprof.SetMaximum(1.2); eneCprof.Draw()
        ctestC.SaveAs("testC{0}GeV.png".format(energy))

        cprofS = ROOT.TCanvas("cprofS{0}".format(energy), "cprofS{0}".format(energy), 1400, 1200)
        SciEneVsYDWCprof.SetMaximum(energy*1.2); SciEneVsYDWCprof.SetMinimum(energy*0.8)
        SciEneVsYDWCprof.Draw()
        SpmtVsYDWCprof.Draw("same")
        cprofS.SaveAs("SciEneProfY{0}.png".format(energy))

        cprofC = ROOT.TCanvas("cprofC{0}".format(energy), "cprofC{0}".format(energy), 1400, 1200)
        CerEneVsYDWCprof.SetMaximum(energy*1.2); CerEneVsYDWCprof.SetMinimum(energy*0.8)
        CerEneVsYDWCprof.Draw()
        CpmtVsYDWCprof.Draw("same")
        cprofC.SaveAs("EneProfY{0}.png".format(energy))



        # Fill histograms with corrected energy values
        EneHistS = ROOT.TH1D("EneHistS_{0}GeV".format(energy), "EneHistS_{0}GeV; E [GeV]; Counts".format(energy), 100, energy-0.4*energy, energy+0.4*energy)
        EneHistC = ROOT.TH1D("EneHistC_{0}GeV".format(energy), "EneHistC_{0}GeV; E [GeV]; Counts".format(energy), 100, energy-0.4*energy, energy+0.4*energy)
        EneHistComb = ROOT.TH1D("EneHistComb_{0}GeV".format(energy), "EneHistCombined_{0}GeV; E [GeV]; Counts".format(energy), 100, energy-0.4*energy, energy+0.4*energy)


        # filter tree using only events with an asymmetry within range
        AsymCut = 0.5
        print(data.shape)
        data = data[ (np.abs(data["AsymS"]<AsymCut) ) & (np.abs(data["AsymC"]<AsymCut) ) ]
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