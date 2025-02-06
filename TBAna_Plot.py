import pandas as pd
import numpy as np
import ROOT
import argparse



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



    # Linearity Plot
    c_linearity = ROOT.TCanvas("c_linearity", "c_linearity", 1400, 1200)
    c_linearity.SetLeftMargin(0.15); c_linearity.SetRightMargin(0.03)
    c_linearity.SetTopMargin(0.07)


    LinearityMultiGraph = ROOT.TMultiGraph()



    LinGraph_S = ROOT.TGraphErrors(len(Energy), Energy, (Mean_S-Energy)/Energy, 0, MeanErr_S/Energy)
    LinGraph_S.SetLineColor(ROOT.kRed); LinGraph_S.SetLineWidth(2)

    LinGraph_C = ROOT.TGraphErrors(len(Energy), Energy, (Mean_C-Energy)/Energy, 0, MeanErr_C/Energy)
    LinGraph_C.SetLineColor(ROOT.kBlue); LinGraph_C.SetLineWidth(2)

    LinGraph_Comb = ROOT.TGraphErrors(len(Energy), Energy, (Mean_Comb-Energy)/Energy, 0, MeanErr_Comb/Energy)
    LinGraph_Comb.SetLineColor(ROOT.kGreen+1); LinGraph_Comb.SetLineWidth(2)


    LinearityMultiGraph.Add(LinGraph_S, "PL")
    LinearityMultiGraph.Add(LinGraph_C, "PL")
    LinearityMultiGraph.Add(LinGraph_Comb, "PL")



    LinearityMultiGraph.SetTitle(f"Linearity ({particle})")
    LinearityMultiGraph.GetXaxis().SetTitle(r"E_{beam} [GeV]"); LinearityMultiGraph.GetYaxis().SetTitle(r"E_{beam}-E_{reco} / E_{beam}")

    LinearityMultiGraph.SetMinimum(-0.04)
    LinearityMultiGraph.SetMaximum(0.04)
    LinearityMultiGraph.Draw("AP")


    line1=ROOT.TLine(LinearityMultiGraph.GetXaxis().GetXmin(), 0.01, LinearityMultiGraph.GetXaxis().GetXmax(), 0.01)
    line2=ROOT.TLine(LinearityMultiGraph.GetXaxis().GetXmin(), -0.01, LinearityMultiGraph.GetXaxis().GetXmax(), -0.01)
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
    legend.AddEntry(LinGraph_Comb, "Combined", "PL")

    legend.Draw("same")
    c_linearity.SaveAs(f"LinearityPlotTest_{particle}.png")


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

    RMSGraph_S = ROOT.TGraphErrors(len(Energy), 1/np.sqrt(Mean_S), SigmaOverE_S, 0, SigmaOverE_err_S)
    RMSGraph_S.SetMarkerColor(ROOT.kRed)
    RMSGraph_S.SetLineColor(ROOT.kRed); 
    RMSGraph_S.Fit("pol1")
    fit_S = RMSGraph_S.GetFunction("pol1")
    fit_S.SetLineWidth(2); fit_S.SetLineColor(ROOT.kRed)
    p0_S=fit_S.GetParameter(0)*100
    p1_S=fit_S.GetParameter(1)*100

    RMSGraph_C = ROOT.TGraphErrors(len(Energy), 1/np.sqrt(Mean_C), SigmaOverE_C, 0, SigmaOverE_err_C)
    RMSGraph_C.SetMarkerColor(ROOT.kBlue)
    RMSGraph_C.SetLineColor(ROOT.kBlue); 
    RMSGraph_C.Fit("pol1")
    fit_C = RMSGraph_C.GetFunction("pol1")
    fit_C.SetLineWidth(2); fit_C.SetLineColor(ROOT.kBlue)
    p0_C=fit_C.GetParameter(0)*100
    p1_C=fit_C.GetParameter(1)*100

    RMSGraph_Comb = ROOT.TGraphErrors(len(Energy), 1/np.sqrt(Mean_Comb), SigmaOverE_Comb, 0, SigmaOverE_err_Comb)
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
    RMSGraph_S.SetMarkerColor(ROOT.kRed)
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
    RMSGraph_C.SetMarkerColor(ROOT.kBlue)
    RMSGraph_C.SetLineColor(ROOT.kBlue); 
    RMSGraph_C.Fit(FitFunction)
    fit_C = RMSGraph_C.GetFunction("FitFunc")
    fit_C.SetLineWidth(2); fit_C.SetLineColor(ROOT.kBlue)
    p0_C=fit_C.GetParameter(0)*100
    p1_C=fit_C.GetParameter(1)*100

    #RMSGraph_Comb = ROOT.TGraphErrors(len(Energy), 1/np.sqrt(Mean_Comb), SigmaOverE_Comb, 0, SigmaOverE_err_Comb)
    RMSGraph_Comb = ROOT.TGraphErrors(len(Energy), 1/np.sqrt(Energy), SigmaOverE_Comb, 0, SigmaOverE_err_Comb)
    RMSGraph_Comb.SetMarkerColor(ROOT.kGreen+1)
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











if __name__=="__main__":
    main()
