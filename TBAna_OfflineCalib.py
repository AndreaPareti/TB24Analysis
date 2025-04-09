import ROOT
import uproot
import pandas as pd
import numpy as np
from scipy import optimize



infolder = "/home/storage/data_apareti/TB24/ElectronEnergyScan/"
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

        #print(filename, data)
        for pmtS, pmtC in zip(data["totPMTSene"], data["totPMTCene"]):
            PmtNoiseHistS.Fill(pmtS)
            PmtNoiseHistC.Fill(pmtC)

    cNoise = ROOT.TCanvas("cNoise", "cNoise", 1400, 1200)
    PmtNoiseHistS.SetLineColor(ROOT.kRed); PmtNoiseHistS.SetLineWidth(2);  PmtNoiseHistS.Draw()      
    PmtNoiseHistC.SetLineColor(ROOT.kBlue); PmtNoiseHistC.SetLineWidth(2);  PmtNoiseHistC.Draw("same")

    PmtNoiseHistS.Fit("gaus"); PmtNoiseHistC.Fit("gaus")
    RmsNoiseS = PmtNoiseHistS.GetFunction("gaus").GetParameter(2)
    RmsNoiseC = PmtNoiseHistC.GetFunction("gaus").GetParameter(2)
    #print("rms S: ", RmsNoiseS, "\trms C: ", RmsNoiseC)

    cNoise.Draw()
    cNoise.SaveAs("Noise.png")
    return(RmsNoiseS, RmsNoiseC)


# Function to get mean Y from TProfile for any X value
def get_profile_mean(profile, x):
    bin_number = profile.FindBin(x)  # Find closest bin
    return profile.GetBinContent(bin_number)  # Get mean Y


def GetDFparametrization(run, Cut, filename, energy): 
    print("Using file: ", filename, energy)
    root_file = uproot.open(infolder+filename)
    tree = root_file[treename]
    data = tree.arrays(cut=Cut, library="pd")

    # define Asymmetry variable 
    data["AsymS"] = (data["TS11"] - data["TS15"]) / (data["TS11"] + data["TS15"] )
    data["AsymC"] = (data["TC11"] - data["TC15"]) / (data["TC11"] + data["TC15"] )

    # define partial barycenter variable 
    data["BaryS"] = (data["TS00"]-28.3*data["TS11"]+28.3*data["TS15"])/(data["TS00"]+data["TS11"]+data["TS15"])
    data["BaryC"] = (data["TC00"]-28.3*data["TC11"]+28.3*data["TC15"])/(data["TC00"]+data["TC11"]+data["TC15"])

    # input truth energy to dataframe
    data["TruthE"] = energy

    # Profile normalised energy Vs Asymmetry
    eneSprof = ROOT.TProfile("eneSprof_{0}GeV".format(energy), "S Energy profile over Asymmetry {0}GeV; TS11-TS15 / TS11 + TS15; totPMTSene/E".format(energy), 100, -1, 1)
    eneCprof = ROOT.TProfile("eneCprof_{0}GeV".format(energy), "C Energy profile over Asymmetry {0}GeV; TC11-TC15 / TC11 + TC15; totPMTCene/E".format(energy), 100, -1, 1)

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
    #funcSvector = np.vectorize(funcS.Eval)
    #funcCvector = np.vectorize(funcC.Eval)
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
    print("Hello world!")

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

    run = "0766"
    calib_energy = 20
    filename = "physics_sps2024_run" + run + ".root"
    print(filename, calib_energy)
    ROOT.gStyle.SetOptFit(0)
    index=2
    print("Current cuts on DWCs: ", cut_x_min[index], cut_x_max[index], cut_y_min[index], cut_y_max[index])
    CurrentCut = myCut + " & (XDWC2 > {0}) & (XDWC2 < {1}) & (YDWC2 > {2}) & (YDWC2 < {3})".format(cut_x_min[index], cut_x_max[index], cut_y_min[index], cut_y_max[index]) 

    #df, funcS, funcC, eneSprof, eneCprof = GetDFparametrization(run, CurrentCut, filename, energy)
    root_file = uproot.open(infolder+filename)
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


    pmtShist = ROOT.TH1D("pmtSraw", "Sci PMT energy run {0}".format(run), 36, 1, 36)
    pmtSweights = ROOT.TH1D("pmtSweights", "Sci PMT weights run {0}".format(run), 36, 1, 36)
    for index, tow in enumerate(inputs_name_S):
        print(index, tow, data[tow].mean(), weights_Sci[index])
        pmtShist.Fill(index+1, data[tow].mean())
        pmtShist.GetXaxis().SetBinLabel(index+1, tow)
        pmtSweights.SetBinContent(index+1, weights_Sci[index])
        pmtSweights.GetXaxis().SetBinLabel(index+1, tow)

    pmtChist = ROOT.TH1D("pmtCraw", "Cer PMT energy run {0}".format(run), 36, 1, 36)
    pmtCweights = ROOT.TH1D("pmtCweightC", "Cer PMT weights run {0}".format(run), 36, 1, 36)
    for index, tow in enumerate(inputC_name_C):
        print(index, tow, data[tow].mean(), weights_Cer[index])
        pmtChist.Fill(index+1, data[tow].mean())
        pmtChist.GetXaxis().SetBinLabel(index+1, tow)
        pmtCweights.SetBinContent(index+1, weights_Cer[index])
        pmtCweights.GetXaxis().SetBinLabel(index+1, tow)


    ROOT.gStyle.SetOptStat(0)
    cCanvaSci = ROOT.TCanvas("cSci", "cSci", 1600, 600)
    pmtSweights.SetFillColor(ROOT.kRed); pmtSweights.SetLineWidth(0)
    pmtShist.SetLineWidth(2); pmtShist.SetLineColor(ROOT.kBlue+1)
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
    leg = ROOT.TLegend(0.80, 0.75, 0.95, 0.9)
    leg.AddEntry(pmtShist, "PMT energy before re-calibration", "l")
    leg.AddEntry(pmtSweights, "Weights associated to PMTs", "f")
    leg.Draw()

    cCanvaSci.SaveAs("testS.pdf")

    cCanvaCer = ROOT.TCanvas("cCer", "cCer", 1600, 600)
    pmtCweights.SetFillColor(ROOT.kCyan); pmtCweights.SetLineWidth(0)
    pmtChist.SetLineWidth(2); pmtChist.SetLineColor(ROOT.kBlue+1)
    pmtChist.SetLineWidth(2); pmtChist.SetLineColor(ROOT.kBlue+1)
    pmtChist.Scale(1/pmtCweights.GetMaximum())
    #pmtCweights.SetMaximum( max(pmtCweights.GetMaximum(), pmtChist.GetMaximum())*1.15 )
    pmtCweights.SetMaximum( 7 )
    line = ROOT.TLine(pmtChist.GetXaxis().GetXmin(), 1.0, pmtChist.GetXaxis().GetXmax(), 1.0)
    line.SetLineColor(ROOT.kGreen)  # Set line color to red
    line.SetLineWidth(2); line.SetLineStyle(2)          # Set line thickness    
    print("Max weight C: ", pmtCweights.GetMaximum())
    #pmtCweights.SetMaximum( pmtCweights.GetMaximum()*1.5 )
    pmtCweights.Draw("hist")
    pmtChist.Draw("sameHIST")
    line.Draw("same")
    leg = ROOT.TLegend(0.80, 0.75, 0.95, 0.9)
    leg.AddEntry(pmtChist, "PMT energy before re-calibration", "l")
    leg.AddEntry(pmtCweights, "Weights associated to PMTs", "f")
    leg.Draw()


    cCanvaCer.SaveAs("testC.pdf")

    ROOT.gStyle.SetOptStat(111)


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


    runs = ["0786", "0766", "0772", "0774", "0775", "0778", "0779", "0792"]
    #runs = ["0786", "0999", "0772", "0774", "0775", "0778", "0779", "0792"]

    #runs = ["0786", "0766", "0772", "0774", "0775", "0776", "0779", "0792"]
    #runs = ["0786", "1018", "0772", "0774", "0775", "0776", "0779", "0792"] 

    #energies = [10, 20, 30, 40, 60, 80, 100, 120]

    runs = ["0786", "0766", "0772", "0774", "0775", "0778", "0779"]
    energies = [10, 20, 30, 40, 60, 80, 100]


    #myCut = "(abs(XDWC2 - XDWC1) < 5) & (abs(YDWC2 - YDWC1)<5) & (MCounter<200) & (TailC<300) & (C2>160) & ( (totLeakage - L20)<5000 ) & (PShower>550)"
    myCut = "(abs(XDWC2 - XDWC1) < 5) & (abs(YDWC2 - YDWC1)<5) & (MCounter<200) & (TailC<300) & (C2>160) & ( (totLeakage - L20)<5000 )"

    varProf = "YDWC2"
    cut_x_min = [-19.83, -16.74, -16.22, -15.95, -15.60, -16.12, -16.07, -15.50]
    cut_x_max = [23.90, 22.19, 23.27, 23.44, 24.27, 23.79, 23.63, 24.12]
    cut_y_min = [-26.54, -25.80, -26.15, -26.15, -26.39, -25.63, -25.63, -26.03]
    cut_y_max = [13.38, 10.89, 9.72, 9.50, 9.86, 10.89, 10.54, 10.17]


    #cut_y_min = [-20, -20, -20, -20, -20, -20, -20, -20]
    #cut_y_max = [4, 4, 4, 4, 4, 4, 4, 4]


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

        data, funcS, funcC, eneSprof, eneCprof = GetDFparametrization(run, CurrentCut, filename, energy)
        # remove events on which the calibration was performed
        if(energy==calib_energy):
            data = data.iloc[cutCalib:]

        # Apply and test recalibration
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


        data["totPMTSene"] = np.dot(A_S.T, weights_Sci)
        data["totPMTCene"] = np.dot(A_C.T, weights_Cer)


        dfs.append(data)
        FitS.append(funcS)
        FitC.append(funcC)
        profS.append(eneSprof)
        profC.append(eneCprof)
        ROOT.gStyle.SetOptStat(0)
        ctest = ROOT.TCanvas("cEneProfAsymmetry{0}".format(energy), "cEneProfAsymmetry{0}".format(energy), 1400, 1200)
        ctest.SetLeftMargin(0.15)
        eneSprof.GetYaxis().SetTitle("PMT energy / Beam nominal energy")
        eneSprof.SetLineWidth(2); eneSprof.SetLineColor(ROOT.kRed); eneSprof.SetMarkerStyle(ROOT.kFullCircle); eneSprof.SetMarkerColor(ROOT.kRed)
        eneSprof.SetMinimum(0.8); eneSprof.SetMaximum(1.2); eneSprof.Draw()
        #ctestS.SaveAs("testS{0}GeV.png".format(energy))

        #ctestC = ROOT.TCanvas("c{0}".format(energy), "c{0}".format(energy), 1400, 1200)
        eneCprof.SetLineWidth(2); eneCprof.SetLineColor(ROOT.kBlue); eneCprof.SetMarkerStyle(ROOT.kFullCircle); eneCprof.SetMarkerColor(ROOT.kBlue)
        eneCprof.SetMinimum(0.8); eneCprof.SetMaximum(1.2); eneCprof.Draw("same")

        leg = ROOT.TLegend(0.65, 0.15, 0.85, 0.25)
        leg.AddEntry(eneSprof, "S channel", "L")
        leg.AddEntry(eneCprof, "C channel", "L")
        leg.Draw()
        ctest.SaveAs("EnergyProfOverAsymmetry_{0}GeV.png".format(energy))        


    fS10  = np.vectorize(FitS[0].Eval); fC10  = np.vectorize(FitC[0].Eval)  
    fS20  = np.vectorize(FitS[1].Eval); fC20  = np.vectorize(FitC[1].Eval)
    fS30  = np.vectorize(FitS[2].Eval); fC30  = np.vectorize(FitC[2].Eval)
    fS40  = np.vectorize(FitS[3].Eval); fC40  = np.vectorize(FitC[3].Eval)
    fS60  = np.vectorize(FitS[4].Eval); fC60  = np.vectorize(FitC[4].Eval)
    fS80  = np.vectorize(FitS[5].Eval); fC80  = np.vectorize(FitC[5].Eval)
    fS100  = np.vectorize(FitS[6].Eval); fC100  = np.vectorize(FitC[6].Eval)
    #fS120  = np.vectorize(FitS[7].Eval); fC120  = np.vectorize(FitC[7].Eval)

    binning_S = [0, 25, 53, 97]
    binning_C = [0, 25, 53, 97]


    


    DrawEnergyHist(dfs, energies, binning_S, binning_C, "totPMTSene", "totPMTCene")


    noiseS, noiseC = GetPMTnoise()


    

    ROOT.gStyle.SetOptStat(0)
        
    for index, energy in enumerate(energies):
        print("Index ", index, "\tUsing file with energy: ", energy)

        data = dfs[index]

        # Define conditions
        #binning_S = [0, 15, 25, 35, 50, 75, 100, 120]
        #binning_C = [0, 15, 25, 35, 50, 75, 100, 120]
        binning_S = [0, 25, 53, 97]
        binning_C = [0, 25, 53, 97]


        # Define binnings on PMT energy
        #conditions_S = [
        #    (data["totPMTSene"] >= binning_S[0]) & (data["totPMTSene"] < binning_S[1]),
        #    (data["totPMTSene"] >= binning_S[1]) & (data["totPMTSene"] < binning_S[2]),
        #    (data["totPMTSene"] >= binning_S[2]) & (data["totPMTSene"] < binning_S[3]),
        #    (data["totPMTSene"] >= binning_S[3])
        #]
        conditions_S = [data["totPMTSene"].any()>0]
        conditions_C = [data["totPMTCene"].any()>0]


        #conditions_C = [
        #    (data["totPMTCene"] >= binning_C[0]) & (data["totPMTCene"] < binning_C[1]),
        #    (data["totPMTCene"] >= binning_C[1]) & (data["totPMTCene"] < binning_C[2]),
        #    (data["totPMTCene"] >= binning_C[2]) & (data["totPMTCene"] < binning_C[3]),
        #    (data["totPMTCene"] >= binning_C[3])
        #]

        # Associate a parametrisation to each energy point
        # A parametrisation at each energy is extracted, it is possible to use only one for all
        # or one per bin
        #choices_S = [
        #    data["totPMTSene"] / fS40(data["AsymS"]),
        #    data["totPMTSene"] / fS40(data["AsymS"]),
        #    data["totPMTSene"] / fS40(data["AsymS"]),
        #    data["totPMTSene"] / fS40(data["AsymS"])
        #]
        #choices_C = [
        #    data["totPMTCene"] / fC40(data["AsymC"]),
        #    data["totPMTCene"] / fC40(data["AsymC"]),
        #    data["totPMTCene"] / fC40(data["AsymC"]),
        #    data["totPMTCene"] / fC40(data["AsymC"])
        #]
        choices_S = [data["totPMTSene"]/fS100(data["AsymS"])]
        choices_C = [data["totPMTCene"]/fC100(data["AsymC"])/0.95]
        #choices_C = [data["totPMTCene"]/get_profile_mean(profC[3], data["AsymC"])]

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
        AsymCut = 0.5
        #data = data[ (np.abs(data["AsymS"]<AsymCut) ) & (np.abs(data["AsymC"]<AsymCut) ) ]
        #data_filtered = data[ (np.abs(data["AsymS"])<AsymCut ) & (np.abs(data["AsymC"])<AsymCut ) & (np.abs(data["BaryS"])<4) & (np.abs(data["BaryC"])<4) ]
        data_filtered = data[ (np.abs(data["AsymS"])<AsymCut ) & (np.abs(data["AsymC"])<AsymCut ) ]

        for eneS, eneC in zip(data_filtered["energyS"].values, data_filtered["energyC"].values):
            HistScorrected_Asymcut.Fill(eneS)
            HistCcorrected_Asymcut.Fill(eneC)
            HistCombCorrected.Fill( (eneS+eneC)/2 )




        
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

        # Profile Barycenter position with asymmetry
        SpmtAsymBarHist = ROOT.TProfile("SpmtAsymBarHist_{0}".format(energy), "BarycenterY Vs Asymmetry (S channel) ({0} GeV); Asymmetry; Y Barycenter [mm]".format(energy), 50, -1, 1)
        CpmtAsymBarHist = ROOT.TProfile("CpmtAsymBarHist_{0}".format(energy), "BarycenterY Vs Asymmetry (C channel) ({0} GeV); Asymmetry; Y Barycenter [mm]".format(energy), 50, -1, 1)



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


        cAsymBarHist = ROOT.TCanvas("SpmtAsymBarHist_{0}".format(energy), "SpmtAsymBarHist_{0}".format(energy), 1400, 1200)
        SpmtAsymBarHist.SetLineColor(ROOT.kRed); SpmtAsymBarHist.SetMarkerStyle(23); SpmtAsymBarHist.SetMarkerColor(ROOT.kRed)
        SpmtAsymBarHist.Draw("")
        #cSpmtAsymBarHist.SaveAs("SpmtAsymBarHist_{0}GeV.png".format(energy))
        #cCpmtAsymBarHist = ROOT.TCanvas("CpmtAsymBarHist_{0}".format(energy), "CpmtAsymBarHist_{0}".format(energy), 1400, 1200)
        CpmtAsymBarHist.SetLineColor(ROOT.kBlue); CpmtAsymBarHist.SetMarkerStyle(23); CpmtAsymBarHist.SetMarkerColor(ROOT.kBlue)
        CpmtAsymBarHist.Draw("same")
        leg = ROOT.TLegend(0.73, 0.75, 0.88, 0.88)
        leg.AddEntry(SpmtAsymBarHist, "S Channel", "PL"); leg.AddEntry(CpmtAsymBarHist, "C Channel", "PL")
        leg.SetTextSize(0.019)
        leg.Draw()
        cAsymBarHist.SaveAs("AsymBarProfile_{0}GeV.png".format(energy))
        






        cS = ROOT.TCanvas("cS{0}".format(energy), "cS{0}".format(energy), 1400, 1200)
        cS.SetLeftMargin(0.15)
        cS.SetRightMargin(0.07)
        HistSrawNorm = HistSraw.Clone()
        HistScorrectedNorm = HistScorrected.Clone()
        HistScorrected_AsymcutNorm = HistScorrected_Asymcut.Clone()
        HistSrawNorm.Scale(1/HistSraw.Integral())
        HistScorrectedNorm.Scale(1/HistScorrectedNorm.Integral())
        HistScorrected_AsymcutNorm.Scale(1/HistScorrected_AsymcutNorm.Integral())

        HistSrawNorm.SetLineWidth(2); HistSrawNorm.SetLineColor(ROOT.kBlue+1); HistSrawNorm.SetFillColorAlpha(ROOT.kAzure+10, 0.08)
        HistScorrectedNorm.SetLineWidth(2); HistScorrectedNorm.SetLineColor(ROOT.kRed)
        HistScorrected_AsymcutNorm.SetLineWidth(2); HistScorrected_AsymcutNorm.SetLineColor(ROOT.kGreen+2)
        max_value = max(HistScorrectedNorm.GetMaximum(), HistSrawNorm.GetMaximum(), HistScorrected_AsymcutNorm.GetMaximum())*1.15
        HistSrawNorm.SetTitle("Reco S energy at {0}GeV".format(energy))
        HistSrawNorm.GetYaxis().SetTitle("Normalised Counts")
        HistSrawNorm.SetMaximum(max_value)
        HistSrawNorm.Draw("hist")
        HistScorrectedNorm.Draw("same hist")
        HistScorrected_AsymcutNorm.Draw("same hist")
        leg = ROOT.TLegend(0.62, 0.75, 0.92, 0.88)
        leg.SetTextSize(0.018)
        leg.AddEntry(HistSrawNorm, "totPMTSene (raw)", "L")
        leg.AddEntry(HistScorrectedNorm, "S energy (corrected)", "L")
        leg.AddEntry(HistScorrected_AsymcutNorm, "S energy (corrected, |asym|<{0})".format(AsymCut), "L")
        leg.Draw()
        cS.SaveAs("SciEnergy{0}GeV.png".format(energy))
        


        cC = ROOT.TCanvas("cC{0}".format(energy), "cC{0}".format(energy), 1400, 1200)
        cC.SetLeftMargin(0.15)
        cC.SetRightMargin(0.07)
        HistCrawNorm = HistCraw.Clone()
        HistCcorrectedNorm = HistCcorrected.Clone()
        HistCcorrected_AsymcutNorm = HistCcorrected_Asymcut.Clone()
        HistCrawNorm.Scale(1/HistCraw.Integral())
        HistCcorrectedNorm.Scale(1/HistCcorrectedNorm.Integral())
        HistCcorrected_AsymcutNorm.Scale(1/HistCcorrected_AsymcutNorm.Integral())

        HistCrawNorm.SetLineWidth(2); HistCrawNorm.SetLineColor(ROOT.kBlue+1); HistCrawNorm.SetFillColorAlpha(ROOT.kAzure+10, 0.08)
        HistCcorrectedNorm.SetLineWidth(2); HistCcorrectedNorm.SetLineColor(ROOT.kRed)
        HistCcorrected_AsymcutNorm.SetLineWidth(2); HistCcorrected_AsymcutNorm.SetLineColor(ROOT.kGreen+2)
        max_value = max(HistCcorrectedNorm.GetMaximum(), HistCrawNorm.GetMaximum(), HistCcorrected_AsymcutNorm.GetMaximum())*1.15
        HistCrawNorm.SetTitle("Reco C energy at {0}GeV".format(energy))
        HistCrawNorm.GetYaxis().SetTitle("Normalised Counts")
        HistCrawNorm.SetMaximum(max_value)
        HistCrawNorm.Draw("hist")
        HistCcorrectedNorm.Draw("same hist")
        HistCcorrected_AsymcutNorm.Draw("same hist")
        leg = ROOT.TLegend(0.62, 0.75, 0.92, 0.88)
        leg.SetTextSize(0.018)
        leg.AddEntry(HistCrawNorm, "totPMTCene (raw)", "L")
        leg.AddEntry(HistCcorrectedNorm, "C energy (corrected)", "L")
        leg.AddEntry(HistCcorrected_AsymcutNorm, "C energy (corrected, |asym|<{0})".format(AsymCut), "L")
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

        # With known RMSs for S and C channels, combine resolution
        # weighting for the relative resolutions

        sigmaS = BestFitS.GetParameter(2)
        sigmaC = BestFitC.GetParameter(2)

        weightedEneHist = ROOT.TH1D("WeightedEne_{0}".format(energy), "Combined Energy (weighted for RMS, {0}GeV);E [GeV];Counts".format(energy), 100, energy-0.4*energy, energy+0.4*energy)
        for eneS, eneC in zip(data["energyS"], data["energyC"]):
            s2 = sigmaS*sigmaS
            c2 = sigmaC*sigmaC
            num = (eneS/s2) + (eneC/c2)
            denum = (1/s2) + (1/c2)
            weightedE = num/denum
            weightedEneHist.Fill(weightedE)


        weightedEneHist.Fit("gaus","Q")
        fitWeightEne = weightedEneHist.GetFunction("gaus")
        weightedEneHist.Fit("gaus", "Q", "", fitWeightEne.GetParameter(1)-1.5*fitWeightEne.GetParameter(2), fitWeightEne.GetParameter(1)+3*fitWeightEne.GetParameter(2)) 
        BestFitWeighted = weightedEneHist.GetFunction("gaus")

        cWeight = ROOT.TCanvas("cWeightedEnergy{0}".format(energy), "cWeightedEnergy{0}".format(energy), 1400, 1200)
        weightedEneHist.Draw()
        cWeight.SaveAs("WeightedEne{0}.png".format(energy))


        MeanVec_S.append(BestFitS.GetParameter(1)); MeanErrVec_S.append(BestFitS.GetParError(1)); RmsVec_S.append(BestFitS.GetParameter(2)); RmsErrVec_S.append(BestFitS.GetParError(2))
        MeanVec_C.append(BestFitC.GetParameter(1)); MeanErrVec_C.append(BestFitC.GetParError(1)); RmsVec_C.append(BestFitC.GetParameter(2)); RmsErrVec_C.append(BestFitC.GetParError(2))
        MeanVec_Comb.append(BestFitComb.GetParameter(1)); MeanErrVec_Comb.append(BestFitComb.GetParError(1)); RmsVec_Comb.append(BestFitComb.GetParameter(2)); RmsErrVec_Comb.append(BestFitComb.GetParError(2))
        


        #MeanVec_Comb.append(BestFitWeighted.GetParameter(1)); MeanErrVec_Comb.append(BestFitWeighted.GetParError(1)); RmsVec_Comb.append(BestFitWeighted.GetParameter(2)); RmsErrVec_Comb.append(BestFitWeighted.GetParError(2))


    DrawEnergyHist(dfCorrected_array, energies, binning_S, binning_C, "energyS", "energyC")

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

                
    
    










if __name__=="__main__":
    main()