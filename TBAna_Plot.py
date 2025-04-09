import pandas as pd
import numpy as np
import ROOT
import argparse
import uproot

pedfolder = "/home/storage/data_apareti/TB24/PedestalRuns/"
treename = "Ftree"

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



def main():
    print("Hello there")
    # Create the argument parser
    parser = argparse.ArgumentParser(description="Input particles")
    
    # Add arguments
    parser.add_argument("--particle", type=str, required=True, help="Name of the particle.")
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Access the `particle` argument
    particle = args.particle
    print(f"Particle name provided: {particle}")

    df = pd.read_csv("DRcaloData.csv", header=0, sep="\t")  # Load CSV
    print(df)
    #columns = [df[i].to_numpy() for i in df.columns]  # Convert each column to a NumPy array

    # Dynamically create NumPy arrays with names matching DataFrame column names
    for column in df.columns:
        globals()[column] = df[column].to_numpy().astype(float)

    #column names
    #Energy      Mean_S     RMS_S  MeanErr_S  RMSErr_S      Mean_C     RMS_C  MeanErr_C  RMSErr_C   Mean_Comb  RMS_Comb  MeanErr_Comb  RMSErr_Comb

    # Read PMT noise
    noiseS, noiseC = GetPMTnoise()
    # Get Noise value and make it an array
    noiseSvec = np.full(len(RMS_S), noiseS)
    RmsVec_S_corrected = np.sqrt(np.asarray(RMS_S)**2 - noiseSvec**2)
    print("S resolution before noise subtraction: ", RMS_S)
    print("S resolution after noise subtraction: ", RmsVec_S_corrected)

    noiseCvec = np.full(len(RMS_C), noiseC)
    RmsVec_C_corrected = np.sqrt(np.asarray(RMS_C)**2 - noiseCvec**2)
    print("C resolution before noise subtraction: ", RMS_C)
    print("C resolution after noise subtraction: ", RmsVec_C_corrected)

    # take pedestal S+C/2 and subtract it from combined energy
    noiseCombVec = np.full(len(RMS_Comb), (noiseS+noiseC)/2)
    RmsVec_Comb_corrected = np.sqrt(np.asarray(RMS_Comb)**2 - noiseCombVec**2)
    print("Combined resolution before noise subtraction: ", RMS_Comb)
    print("Combined resolution after noise subtraction: ", RmsVec_Comb_corrected)


    # Linearity Plot
    c_linearity = ROOT.TCanvas("c_linearity", "c_linearity", 1400, 1200)
    c_linearity.SetLeftMargin(0.15); c_linearity.SetRightMargin(0.03)
    c_linearity.SetTopMargin(0.07)


    LinearityMultiGraph = ROOT.TMultiGraph()



    LinGraph_S = ROOT.TGraphErrors(len(Energy), Energy, (Mean_S-Energy)*100/Energy, 0, MeanErr_S/Energy)
    LinGraph_S.SetLineColor(ROOT.kRed); LinGraph_S.SetLineWidth(2)
    LinGraph_S.SetMarkerStyle(20); LinGraph_S.SetMarkerSize(2); LinGraph_S.SetMarkerColor(ROOT.kRed)

    LinGraph_C = ROOT.TGraphErrors(len(Energy), Energy, (Mean_C-Energy)*100/Energy, 0, MeanErr_C/Energy)
    LinGraph_C.SetLineColor(ROOT.kBlue); LinGraph_C.SetLineWidth(2)
    LinGraph_C.SetMarkerStyle(20); LinGraph_C.SetMarkerSize(2); LinGraph_C.SetMarkerColor(ROOT.kBlue)    

    LinGraph_Comb = ROOT.TGraphErrors(len(Energy), Energy, (Mean_Comb-Energy)*100/Energy, 0, MeanErr_Comb/Energy)
    LinGraph_Comb.SetLineColor(ROOT.kGreen+2); LinGraph_Comb.SetLineWidth(2)
    LinGraph_Comb.SetMarkerStyle(20); LinGraph_Comb.SetMarkerSize(2); LinGraph_Comb.SetMarkerColor(ROOT.kGreen+2)  

    LinearityMultiGraph.Add(LinGraph_S, "P")
    LinearityMultiGraph.Add(LinGraph_C, "P")
    LinearityMultiGraph.Add(LinGraph_Comb, "P")



    LinearityMultiGraph.SetTitle(f"Linearity ({particle})")
    LinearityMultiGraph.GetXaxis().SetTitle(r"E_{beam} [GeV]"); LinearityMultiGraph.GetYaxis().SetTitle(r"E_{beam}-E_{reco} / E_{beam} (%)")

    #LinearityMultiGraph.SetMinimum(-0.04)
    #LinearityMultiGraph.SetMaximum(0.04)
    LinearityMultiGraph.SetMinimum(-0.4*100)
    LinearityMultiGraph.SetMaximum(0.2*100)

    LinearityMultiGraph.Draw("AP")


    line1=ROOT.TLine(LinearityMultiGraph.GetXaxis().GetXmin(), 0.01*100, LinearityMultiGraph.GetXaxis().GetXmax(), 0.01*100)
    line2=ROOT.TLine(LinearityMultiGraph.GetXaxis().GetXmin(), -0.01*100, LinearityMultiGraph.GetXaxis().GetXmax(), -0.01*100)
    line3=ROOT.TLine(LinearityMultiGraph.GetXaxis().GetXmin(), 0., LinearityMultiGraph.GetXaxis().GetXmax(), 0.)
    line1.SetLineStyle(2); line1.SetLineWidth(2)
    line2.SetLineStyle(2); line2.SetLineWidth(2)
    line3.SetLineStyle(2); line3.SetLineWidth(2)
    line3.SetLineColor(ROOT.kGreen+3)
    line1.Draw("same")
    line2.Draw("same")
    line3.Draw("same")

    alignLeft = 0.2
    alignTop = 0.86
    alignRight = 0.925
    tex0 = ROOT.TLatex(alignLeft, alignTop, "TB24 DRAGO")
    tex0.SetNDC()
    tex0.SetTextFont(72)
    tex0.SetTextSize(0.042)
    tex0.SetLineWidth(2)
    tex0.Draw("same")	
    tex1 = ROOT.TLatex(alignLeft+0.25,alignTop,"Preliminary ")
    tex1.SetNDC()
    tex1.SetTextFont(42)
    tex1.SetTextSize(0.042)
    tex1.SetLineWidth(2)
    tex1.Draw("same")	

    '''tex2 = ROOT.TLatex(alignRight,alignTop+0.025,"Length: 2500 mm ");
    tex2.SetNDC();
    tex2.SetTextFont(42);
    tex2.SetTextSize(0.03);	tex2.SetTextAlign(32);
    tex2.SetLineWidth(2);
    tex2.Draw("same");	
    tex3 = ROOT.TLatex(alignRight,alignTop-0.025,"Correcting for containment");
    tex3.SetNDC();
    tex3.SetTextFont(42);	tex3.SetTextAlign(32);
    tex3.SetTextSize(0.032);
    tex3.SetLineWidth(2);
    tex3.Draw("same");	
    tex4 = ROOT.TLatex(alignRight,alignTop-0.075,"at 40 GeV");
    tex4.SetNDC();
    tex4.SetTextFont(42);	tex4.SetTextAlign(32);
    tex4.SetTextSize(0.032);
    tex4.SetLineWidth(2);
    tex4.Draw("same");	'''


    legend=ROOT.TLegend(alignLeft, alignTop-0.12, alignLeft+0.2, alignTop-0.02)
    legend.SetBorderSize(0)
    legend.SetTextSize(.035)
    legend.SetTextFont(42)
    #legend.SetHeader("HidraSim, 7m Attenuation Length")
    legend.AddEntry(LinGraph_S, "S Channel", "PL")
    legend.AddEntry(LinGraph_C, "C Channel", "PL")
    legend.AddEntry(LinGraph_Comb, "Dual-Readout, Combined", "PL")

    legend.Draw("same")
    c_linearity.SaveAs(f"LinearityPlotTest_{particle}.pdf")


    # Resolution Plot
    #########################################################
    ########### RESOLUTION PLOT #############################
    #########################################################
    ROOT.gStyle.SetOptStat(0)
    c_res=ROOT.TCanvas("c_res", "c_res", 1400, 1200)
    c_res.SetLeftMargin(0.15); c_res.SetRightMargin(0.03)
    c_res.SetTopMargin(0.07)

    SigmaOverE_S = RMS_S/Mean_S; SigmaOverE_C = RMS_C/Mean_C; SigmaOverE_Comb = RMS_Comb/Mean_Comb
    SigmaOverE_err_S = np.sqrt( (MeanErr_S/Mean_S)**2 + (RMSErr_S/RMS_S)**2 )*SigmaOverE_S
    SigmaOverE_err_C = np.sqrt( (MeanErr_C/Mean_C)**2 + (RMSErr_C/RMS_C)**2 )*SigmaOverE_C
    SigmaOverE_err_Comb = np.sqrt( (MeanErr_Comb/Mean_Comb)**2 + (RMSErr_Comb/RMS_Comb)**2 )*SigmaOverE_Comb


    ResMultiGraph = ROOT.TMultiGraph()

    #RMSGraph_S = ROOT.TGraphErrors(len(Energy), 1/np.sqrt(Mean_S), SigmaOverE_S, 0, SigmaOverE_err_S)
    RMSGraph_S = ROOT.TGraphErrors(len(Energy), 1/np.sqrt(Energy), SigmaOverE_S, 0, SigmaOverE_err_S)
    RMSGraph_S.SetMarkerColor(ROOT.kRed)
    RMSGraph_S.SetLineColor(ROOT.kRed); 
    RMSGraph_S.Fit("pol1")
    fit_S = RMSGraph_S.GetFunction("pol1")
    fit_S.SetLineWidth(2); fit_S.SetLineColor(ROOT.kRed)
    p0_S=fit_S.GetParameter(0)*100
    p1_S=fit_S.GetParameter(1)*100


    #RMSGraph_C = ROOT.TGraphErrors(len(Energy), 1/np.sqrt(Mean_C), SigmaOverE_C, 0, SigmaOverE_err_C)
    RMSGraph_C = ROOT.TGraphErrors(len(Energy), 1/np.sqrt(Energy), SigmaOverE_C, 0, SigmaOverE_err_C)
    RMSGraph_C.SetMarkerColor(ROOT.kBlue)
    RMSGraph_C.SetLineColor(ROOT.kBlue); 
    RMSGraph_C.Fit("pol1")
    fit_C = RMSGraph_C.GetFunction("pol1")
    fit_C.SetLineWidth(2); fit_C.SetLineColor(ROOT.kBlue)
    p0_C=fit_C.GetParameter(0)*100
    p1_C=fit_C.GetParameter(1)*100

    #RMSGraph_Comb = ROOT.TGraphErrors(len(Energy), 1/np.sqrt(Mean_Comb), SigmaOverE_Comb, 0, SigmaOverE_err_Comb)
    RMSGraph_Comb = ROOT.TGraphErrors(len(Energy), 1/np.sqrt(Energy), SigmaOverE_Comb, 0, SigmaOverE_err_Comb)
    RMSGraph_Comb.SetMarkerColor(ROOT.kGreen+1)
    RMSGraph_Comb.SetLineColor(ROOT.kGreen+1); 
    RMSGraph_Comb.Fit("pol1")
    fit_Comb = RMSGraph_Comb.GetFunction("pol1")
    fit_Comb.SetLineWidth(2); fit_Comb.SetLineColor(ROOT.kGreen+1)
    p0_Comb=fit_Comb.GetParameter(0)*100
    p1_Comb=fit_Comb.GetParameter(1)*100




    str_S="S Channel: #frac{#sigma(E)}{E} = #frac{"+str(round(p1_S,2))+"%}{#sqrt{E}} + "+str(round(p0_S,2))+"%"
    str_C="C Channel: #frac{#sigma(E)}{E} = #frac{"+str(round(p1_C,2))+"%}{#sqrt{E}} + "+str(round(p0_C,2))+"%"
    str_Comb="Combined: #frac{#sigma(E)}{E} = #frac{"+str(round(p1_Comb,2))+"%}{#sqrt{E}} + "+str(round(p0_Comb,2))+"%"


    #graph1.SetTitle(r"Electron Resolution in [10, 120] GeV Range 1/#sqrt{E} [GeV^{-1/2}] #sigma/E")
    ResMultiGraph.Add(RMSGraph_S, "PE")
    ResMultiGraph.Add(RMSGraph_C, "PE")
    ResMultiGraph.Add(RMSGraph_Comb, "PE")


    ResMultiGraph.SetTitle(f"Energy Resolution in [10, 120] GeV Range ({particle})")
    ResMultiGraph.GetXaxis().SetTitle(r"1/#sqrt{E_{Reco}} [GeV]^{-1/2}")
    ResMultiGraph.GetYaxis().SetTitle(r"#sigma / E_{Reco}")

    ResMultiGraph.Draw("AP")

    alignLeft = 0.20
    alignTop = 0.87
    legend=ROOT.TLegend(alignLeft-0.03, alignTop-0.25, alignLeft+0.2, alignTop-0.025)
    #legend.SetHeader("HiDRaSim, Preliminary", "C")
    legend.SetBorderSize(0)
    legend.SetTextSize(.035)
    legend.SetTextFont(42)
    #legend.AddEntry(RMSGraph_S, r"#chi={fChi:.3f}, Geometry: {fMod} mini-modules".format(fChi=chi, fMod=nMod), "ep")
    legend.AddEntry(fit_S, str_S, "l")
    legend.AddEntry(fit_C, str_C, "l")
    legend.AddEntry(fit_Comb, str_Comb, "l")

    legend.Draw("same")


    tex0 = ROOT.TLatex(alignLeft, alignTop, "TB24 DRAGO")
    tex0.SetNDC()
    tex0.SetTextFont(72)
    tex0.SetTextSize(0.042)
    tex0.SetLineWidth(2)
    tex0.Draw("same")	
    tex1 = ROOT.TLatex(alignLeft+0.25,alignTop,"Preliminary ")
    tex1.SetNDC()
    tex1.SetTextFont(42)
    tex1.SetTextSize(0.042)
    tex1.SetLineWidth(2)
    tex1.Draw("same")	

    c_res.SaveAs(f"ResolutionPlot_{particle}.png")





	############################
	### ERRORI IN QUADRATURA ###
	############################

    c_quad=ROOT.TCanvas("c_quad", "c_quad", 1400, 1200)
    c_quad.SetLeftMargin(0.15); c_quad.SetRightMargin(0.03)
    c_quad.SetTopMargin(0.07)

    SigmaOverE_S = RMS_S/Mean_S; SigmaOverE_C = RMS_C/Mean_C; SigmaOverE_Comb = RMS_Comb/Mean_Comb
    SigmaOverE_err_S = np.sqrt( (MeanErr_S/Mean_S)**2 + (RMSErr_S/RMS_S)**2 )*SigmaOverE_S
    SigmaOverE_err_C = np.sqrt( (MeanErr_C/Mean_C)**2 + (RMSErr_C/RMS_C)**2 )*SigmaOverE_C
    SigmaOverE_err_Comb = np.sqrt( (MeanErr_Comb/Mean_Comb)**2 + (RMSErr_Comb/RMS_Comb)**2 )*SigmaOverE_Comb


    QuadResMultiGraph = ROOT.TMultiGraph()

    #RMSGraph_S = ROOT.TGraphErrors(len(Energy), 1/np.sqrt(Mean_S), SigmaOverE_S, 0, SigmaOverE_err_S)
    RMSGraph_S = ROOT.TGraphErrors(len(Energy), 1/np.sqrt(Energy), SigmaOverE_S, 0, SigmaOverE_err_S)
    RMSGraph_S.SetMarkerColor(ROOT.kRed); RMSGraph_S.SetMarkerStyle(20)
    RMSGraph_S.SetLineColor(ROOT.kRed); 
    FitFunction = ROOT.TF1("FitFunc",  "sqrt( ([1]*x)*([1]*x) + [0]*[0] )")
    FitFunction.SetParameters(0.1, 0.1)

    RMSGraph_S.Fit(FitFunction)

    fit_S = RMSGraph_S.GetFunction("FitFunc")
    fit_S.SetLineWidth(2); fit_S.SetLineColor(ROOT.kRed)
    p0_S=fit_S.GetParameter(0)*100
    p1_S=fit_S.GetParameter(1)*100

    #RMSGraph_C = ROOT.TGraphErrors(len(Energy), 1/np.sqrt(Mean_C), SigmaOverE_C, 0, SigmaOverE_err_C)
    RMSGraph_C = ROOT.TGraphErrors(len(Energy), 1/np.sqrt(Energy), SigmaOverE_C, 0, SigmaOverE_err_C)
    RMSGraph_C.SetMarkerColor(ROOT.kBlue); RMSGraph_C.SetMarkerStyle(20)
    RMSGraph_C.SetLineColor(ROOT.kBlue); 
    RMSGraph_C.Fit(FitFunction)
    fit_C = RMSGraph_C.GetFunction("FitFunc")
    fit_C.SetLineWidth(2); fit_C.SetLineColor(ROOT.kBlue)
    p0_C=fit_C.GetParameter(0)*100
    p1_C=fit_C.GetParameter(1)*100

    #RMSGraph_Comb = ROOT.TGraphErrors(len(Energy), 1/np.sqrt(Mean_Comb), SigmaOverE_Comb, 0, SigmaOverE_err_Comb)
    RMSGraph_Comb = ROOT.TGraphErrors(len(Energy), 1/np.sqrt(Energy), SigmaOverE_Comb, 0, SigmaOverE_err_Comb)
    RMSGraph_Comb.SetMarkerColor(ROOT.kGreen+1); RMSGraph_Comb.SetMarkerStyle(20)
    RMSGraph_Comb.SetLineColor(ROOT.kGreen+1); 
    RMSGraph_Comb.Fit(FitFunction)
    fit_Comb = RMSGraph_Comb.GetFunction("FitFunc")
    fit_Comb.SetLineWidth(2); fit_Comb.SetLineColor(ROOT.kGreen+1)
    p0_Comb=fit_Comb.GetParameter(0)*100
    p1_Comb=fit_Comb.GetParameter(1)*100




    str_S="S Channel: #frac{#sigma(E)}{E} = #frac{"+str(round(p1_S,2))+"%}{#sqrt{E}} #oplus "+str(round(p0_S,2))+"%"
    str_C="C Channel: #frac{#sigma(E)}{E} = #frac{"+str(round(p1_C,2))+"%}{#sqrt{E}} #oplus "+str(round(p0_C,2))+"%"
    str_Comb="Combined: #frac{#sigma(E)}{E} = #frac{"+str(round(p1_Comb,2))+"%}{#sqrt{E}} #oplus "+str(round(p0_Comb,2))+"%"


    #graph1.SetTitle(r"Electron QuadResolution in [10, 120] GeV Range 1/#sqrt{E} [GeV^{-1/2}] #sigma/E")
    QuadResMultiGraph.Add(RMSGraph_S, "PE")
    QuadResMultiGraph.Add(RMSGraph_C, "PE")
    QuadResMultiGraph.Add(RMSGraph_Comb, "PE")


    QuadResMultiGraph.SetTitle(f"Energy Resolution in [10, 120] GeV Range ({particle})")
    QuadResMultiGraph.GetXaxis().SetTitle(r"1/#sqrt{E_{Beam}} [GeV]^{-1/2}")
    QuadResMultiGraph.GetYaxis().SetTitle(r"#sigma / E_{Reco}")

    QuadResMultiGraph.Draw("AP")

    alignLeft = 0.20
    alignTop = 0.87
    legend=ROOT.TLegend(alignLeft-0.03, alignTop-0.25, alignLeft+0.2, alignTop-0.025)
    #legend.SetHeader("HiDRaSim, Preliminary", "C")
    legend.SetBorderSize(0)
    legend.SetTextSize(.035)
    legend.SetTextFont(42)
    #legend.AddEntry(RMSGraph_S, r"#chi={fChi:.3f}, Geometry: {fMod} mini-modules".format(fChi=chi, fMod=nMod), "ep")
    legend.AddEntry(fit_S, str_S, "l")
    legend.AddEntry(fit_C, str_C, "l")
    legend.AddEntry(fit_Comb, str_Comb, "l")

    legend.Draw("same")


    tex0 = ROOT.TLatex(alignLeft, alignTop, "TB24 DRAGO")
    tex0.SetNDC()
    tex0.SetTextFont(72)
    tex0.SetTextSize(0.042)
    tex0.SetLineWidth(2)
    tex0.Draw("same")	
    tex1 = ROOT.TLatex(alignLeft+0.25,alignTop,"Preliminary ")
    tex1.SetNDC()
    tex1.SetTextFont(42)
    tex1.SetTextSize(0.042)
    tex1.SetLineWidth(2)
    tex1.Draw("same")	

    c_quad.SaveAs(f"QuadResolutionPlot_{particle}.png")




    ###### Add Noise term contribution ###
    pmtNoiseS, pmtNoiseC = GetPMTnoise()



	############################
	### ERRORI IN QUADRATURA ###
	############################

    c_quad=ROOT.TCanvas("c_quad_noise", "c_quad_noise", 1400, 1200)
    c_quad.SetLeftMargin(0.15); c_quad.SetRightMargin(0.03)
    c_quad.SetTopMargin(0.07)

    SigmaOverE_S = RMS_S/Mean_S; SigmaOverE_C = RMS_C/Mean_C; SigmaOverE_Comb = RMS_Comb/Mean_Comb
    SigmaOverE_err_S = np.sqrt( (MeanErr_S/Mean_S)**2 + (RMSErr_S/RMS_S)**2 )*SigmaOverE_S
    SigmaOverE_err_C = np.sqrt( (MeanErr_C/Mean_C)**2 + (RMSErr_C/RMS_C)**2 )*SigmaOverE_C
    SigmaOverE_err_Comb = np.sqrt( (MeanErr_Comb/Mean_Comb)**2 + (RMSErr_Comb/RMS_Comb)**2 )*SigmaOverE_Comb


    QuadResMultiGraph = ROOT.TMultiGraph()

    #RMSGraph_S = ROOT.TGraphErrors(len(Energy), 1/np.sqrt(Mean_S), SigmaOverE_S, 0, SigmaOverE_err_S)
    RMSGraph_S = ROOT.TGraphErrors(len(Energy), 1/np.sqrt(Energy), SigmaOverE_S, 0, SigmaOverE_err_S)
    RMSGraph_S.SetMarkerColor(ROOT.kRed); RMSGraph_S.SetMarkerStyle(20)
    RMSGraph_S.SetLineColor(ROOT.kRed); 
    FitFunctionS = ROOT.TF1("FitFuncS",  "sqrt( ([1]*x)*([1]*x) + [0]*[0] + [2]*[2]*x*x*x*x)")
    FitFunctionS.SetParameters(0.1, 0.1, 0.1)
    FitFunctionS.FixParameter(2, pmtNoiseS)

    RMSGraph_S.Fit(FitFunctionS)

    fit_S = RMSGraph_S.GetFunction("FitFuncS")
    fit_S.SetLineWidth(2); fit_S.SetLineColor(ROOT.kRed)
    p0_S=abs(fit_S.GetParameter(0)*100)
    p1_S=abs(fit_S.GetParameter(1)*100)
    p2_S=abs(fit_S.GetParameter(2)*100)

    FitFunctionC = ROOT.TF1("FitFuncC",  "sqrt( ([1]*x)*([1]*x) + [0]*[0] + [2]*[2]*x*x*x*x)")
    FitFunctionC.SetParameters(0.1, 0.1, 0.1)
    FitFunctionC.FixParameter(2, pmtNoiseC)

    #RMSGraph_C = ROOT.TGraphErrors(len(Energy), 1/np.sqrt(Mean_C), SigmaOverE_C, 0, SigmaOverE_err_C)
    RMSGraph_C = ROOT.TGraphErrors(len(Energy), 1/np.sqrt(Energy), SigmaOverE_C, 0, SigmaOverE_err_C)
    RMSGraph_C.SetMarkerColor(ROOT.kBlue); RMSGraph_C.SetMarkerStyle(20)
    RMSGraph_C.SetLineColor(ROOT.kBlue); 
    RMSGraph_C.Fit(FitFunctionC)
    fit_C = RMSGraph_C.GetFunction("FitFuncC")
    fit_C.SetLineWidth(2); fit_C.SetLineColor(ROOT.kBlue)
    p0_C=abs(fit_C.GetParameter(0)*100)
    p1_C=abs(fit_C.GetParameter(1)*100)
    p2_C=abs(fit_C.GetParameter(2)*100)

    # take pedestal S+C/2 and subtract it from combined energy
    #pmtnoiseCombVec = np.full(len(RmsErrVec_C), (noiseS+noiseC)/2)
    #RmsVec_Comb_corrected = np.sqrt(np.asarray(RmsVec_Comb)**2 - noiseCombVec**2)
    pmtNoiseComb = np.sqrt( (pmtNoiseS**2) + (pmtNoiseC**2)  )
 

    FitFunctionComb = ROOT.TF1("FitFuncComb",  "sqrt( ([1]*x)*([1]*x) + [0]*[0] + [2]*[2]*x*x*x*x)")
    FitFunctionComb.SetParameters(0.1, 0.1, 0.1)
    FitFunctionComb.FixParameter(2, pmtNoiseComb)

    #RMSGraph_Comb = ROOT.TGraphErrors(len(Energy), 1/np.sqrt(Mean_Comb), SigmaOverE_Comb, 0, SigmaOverE_err_Comb)
    RMSGraph_Comb = ROOT.TGraphErrors(len(Energy), 1/np.sqrt(Energy), SigmaOverE_Comb, 0, SigmaOverE_err_Comb)
    RMSGraph_Comb.SetMarkerColor(ROOT.kGreen+1); RMSGraph_Comb.SetMarkerStyle(20)
    RMSGraph_Comb.SetLineColor(ROOT.kGreen+1); 
    RMSGraph_Comb.Fit(FitFunctionComb)
    fit_Comb = RMSGraph_Comb.GetFunction("FitFuncComb")
    fit_Comb.SetLineWidth(2); fit_Comb.SetLineColor(ROOT.kGreen+1)
    p0_Comb=abs(fit_Comb.GetParameter(0)*100)
    p1_Comb=abs(fit_Comb.GetParameter(1)*100)
    p2_Comb=abs(fit_Comb.GetParameter(2)*100)




    str_S="S Channel: #frac{#sigma(E)}{E} = #frac{"+str(round(p1_S,2))+"%}{#sqrt{E}} #oplus "+str(round(p0_S,2))+"%"+" #oplus #frac{"+str(round(p2_S, 2))+"%}{E}"
    str_C="C Channel: #frac{#sigma(E)}{E} = #frac{"+str(round(p1_C,2))+"%}{#sqrt{E}} #oplus "+str(round(p0_C,2))+"%"+" #oplus #frac{"+str(round(p2_C, 2))+"%}{E}"
    str_Comb="Combined: #frac{#sigma(E)}{E} = #frac{"+str(round(p1_Comb,2))+"%}{#sqrt{E}} #oplus "+str(round(p0_Comb,2))+"%"+" #oplus #frac{"+str(round(p2_Comb, 2))+"%}{E}"


    #graph1.SetTitle(r"Electron QuadResolution in [10, 120] GeV Range 1/#sqrt{E} [GeV^{-1/2}] #sigma/E")
    QuadResMultiGraph.Add(RMSGraph_S, "PE")
    QuadResMultiGraph.Add(RMSGraph_C, "PE")
    QuadResMultiGraph.Add(RMSGraph_Comb, "PE")


    QuadResMultiGraph.SetTitle(f"Energy Resolution in [10, 120] GeV Range ({particle})")
    QuadResMultiGraph.GetXaxis().SetTitle(r"1/#sqrt{E_{Beam}} [GeV]^{-1/2}")
    QuadResMultiGraph.GetYaxis().SetTitle(r"#sigma / E_{Reco}")

    QuadResMultiGraph.Draw("AP")

    alignLeft = 0.20
    alignTop = 0.87
    legend=ROOT.TLegend(alignLeft-0.03, alignTop-0.25, alignLeft+0.2, alignTop-0.025)
    #legend.SetHeader("HiDRaSim, Preliminary", "C")
    legend.SetBorderSize(0)
    legend.SetTextSize(.031)
    legend.SetTextFont(42)
    #legend.AddEntry(RMSGraph_S, r"#chi={fChi:.3f}, Geometry: {fMod} mini-modules".format(fChi=chi, fMod=nMod), "ep")
    legend.AddEntry(fit_S, str_S, "l")
    legend.AddEntry(fit_C, str_C, "l")
    legend.AddEntry(fit_Comb, str_Comb, "l")

    legend.Draw("same")


    tex0 = ROOT.TLatex(alignLeft, alignTop, "TB24 DRAGO")
    tex0.SetNDC()
    tex0.SetTextFont(72)
    tex0.SetTextSize(0.042)
    tex0.SetLineWidth(2)
    tex0.Draw("same")	
    tex1 = ROOT.TLatex(alignLeft+0.25,alignTop,"Preliminary ")
    tex1.SetNDC()
    tex1.SetTextFont(42)
    tex1.SetTextSize(0.042)
    tex1.SetLineWidth(2)
    tex1.Draw("same")	

    c_quad.SaveAs(f"QuadNoiseResPlot_{particle}.png")










if __name__=="__main__":
    main()
