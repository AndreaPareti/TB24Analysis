import ROOT
import json
import numpy as np
import  pandas as pd 
import  time
#INPUTDIR="/afs/cern.ch/user/i/ideadr/scratch/TB2024_H8/physicsNtuples/"
INPUTDIR="/home/storage/data_apareti/TB24/EqualisationRuns/"


''' Show mean and rms of pedestals for all detectors
    as function of run number
'''
ROOT.gStyle.SetOptStat(0)


def interquatile(x):
    return (np.percentile(x, 75) - np.percentile(x, 25)) / 1.349

def smallestInterval(data, integral=0.6826894921370888):
    """smallestInterval(data, w = none, integral = 0.9) -->
      return the endpoints of the smallest interval
      that contains the given integral (90%) using weights w if given"""
    data = np.asarray(data)
    N = len(data)
    if not N:
       return None
    xsorted = np.sort(data)
    D = int(np.floor(integral * N))
    if D == 0:
        return None
    first_index = (xsorted[D:] - xsorted[:-D]).argmin()
    return xsorted[D + first_index] - xsorted[first_index]




def main():
    print("Hello")
    runs=[557, 558, 641, 681, 682, 694, 713, 767, 770, 771, 777, 777, 780, 781, 782, 783, 784, 796]
    runDictionary={62:727, 52:728, 42:729, 32:730, 22:731, 12:732, 13:734, 14:735, 23:736, 33:737, 43:738, 53:739, 54:740, 44:741, 34:742, 24:743, 15:744, 00:745, 11:747, 21:748, 31:749, 41:750, 51:752, 61:751, 60:753, 50:754, 40:755, 30:756, 20:757, 10:758, 17:759, 16:760, 25:761, 35:762, 45:763, 55:764}
    print(runDictionary[52])

    col1 = [55, 45, 35, 25, 16, 17, 10, 20, 30, 40, 50, 60]
    col2 = [54, 44, 34, 24, 15, 00, 11, 21, 31, 41, 51, 61]
    col3 = [53, 43, 33, 23, 14, 13, 12, 22, 32, 42, 52, 62]
    cols = [col1, col2, col3]

    enemean = []
    towname = []
    towadcpeak = []
    yMaxHitTow = []
    varProf = "YDWC2"
    cut_x_min = [-19.83, -16.74, -16.22, -15.95, -15.60, -16.12, -16.07, -15.50]
    cut_x_max = [23.90, 22.19, 23.27, 23.44, 24.27, 23.79, 23.63, 24.12]
    cut_y_min = [-26.54, -25.80, -26.15, -26.15, -26.39, -25.63, -25.63, -26.03]
    cut_y_max = [13.38, 10.89, 9.72, 9.50, 9.86, 10.89, 10.54, 10.17]

        
    for col in cols:    
        for tow in range(1, len(col2)-1):
            #print("Column: ", col, "\tTow: ", tow)
            print("Tower: ", col[tow], "run: ", runDictionary[col[tow]])
            run = runDictionary[col[tow]]
            rootfile = INPUTDIR+"physics_sps2024_run0"+str(run)+".root"
            print(rootfile)

            #cut = "abs(YDWC2-YDWC1)<5 && abs(XDWC2-XDWC1)<5 && PShower>550 && C2>200 && XDWC2>-20 && XDWC2<20"
            
            
            #cut = "abs(YDWC2-YDWC1)<5 && abs(XDWC2-XDWC1)<5 && PShower>550 && C1>150 && C2>300 && C3>100 && XDWC2>-20 && XDWC2<20"


            myCut = "(abs(XDWC2 - XDWC1) < 5) & (abs(YDWC2 - YDWC1)<5) & (MCounter<200) & (TailC<300) & (C2>160) & ( (totLeakage - L20)<5000 ) & (PShower>550)"


            energy = 20
            index = 1
            cut = myCut + " & (XDWC2 > {0}) & (XDWC2 < {1}) & (YDWC2 > {2}) & (YDWC2 < {3})".format(cut_x_min[index], cut_x_max[index], cut_y_min[index], cut_y_max[index]) 


            #cut = "totPMTSene>40"
            # cut = "abs(YDWC2-YDWC1)<5 && abs(XDWC2-XDWC1)<5 && PShower>550 && C2>200 && YDWC2>-5"

            var0 = "TC"+str(col[tow])
            varUp = "TC"+str(col[tow-1])
            varDown = "TC"+str(col[tow+1])
            if(var0 == "TC0"): var0 = "TC00"
            if(varUp == "TC0"): varUp = "TC00"
            if(varDown == "TC0"): varDown = "TC00"
            varTotS = "totPMTCene"


            print(var0, varUp, varDown)
            #print("\nCreating Profiles...\n")
            

            # Tower where the beam is shot 
            towCenterProf = ROOT.TProfile("towCenterProf{0}".format(tow), "towCenterProf{0}".format(tow), 20, -20, 10)
            towUpProf = ROOT.TProfile("towUpProf{0}".format(tow), "towUpProf{0}".format(tow), 20, -20, 10)
            towDownProf = ROOT.TProfile("towDownProf{0}".format(tow), "towDownProf{0}".format(tow), 20, -20, 10)
            towSumProf = ROOT.TProfile("towSumProf{0}".format(tow), "towSumProf{0}".format(tow), 20, -20, 10)
            totSProf = ROOT.TProfile("towCProf{0}".format(tow), "towCProf{0}".format(tow), 20, -20, 10)

            # Central tower profile, with low-energy events removed
            HitTowProfCut = ROOT.TProfile("HitTowProfCut{0}".format(tow), "HitTowProfCut{0}".format(tow), 20, -20, 10)


            file = ROOT.TFile.Open(rootfile)
            tree = file.Get("Ftree")
            file.cd()
            mean = 0

            # array with energy content of each event
            data = []

            for i in range(tree.GetEntries()):
                tree.GetEntry(i)
                cut_formula = ROOT.TTreeFormula("cut", cut, tree)

                tow0val = tree.GetLeaf(var0).GetValue()
                towUpval = tree.GetLeaf(varUp).GetValue()
                towDownval = tree.GetLeaf(varDown).GetValue()
                totS = tree.GetLeaf(varTotS).GetValue()


                ydwc2 = tree.GetLeaf("YDWC2").GetValue()
                #print("\n value: ",value)
                #value = getattr(tree, key)

                if (cut_formula.EvalInstance()):
                    towCenterProf.Fill(ydwc2, tow0val, 1) 
                    towUpProf.Fill(ydwc2, towUpval, 1)
                    towDownProf.Fill(ydwc2, towDownval, 1)
                    towSumProf.Fill(ydwc2, tow0val+towUpval+towDownval, 1)
                    totSProf.Fill(ydwc2, totS, 1)
                    data.append(tow0val)

                    #if(tow0val>8.):
                    #    HitTowProfCut.Fill(ydwc2, tow0val)
                    #if(ydwc2<0.1 and ydwc2>-1): print("YDWC2==0! TotS: ", totS, "\tYDWC2 value: ", ydwc2)

            #x = smallestInterval(data)
            #interquantile = interquatile(data)

            x0 = np.percentile(data, 80)
            x1 = np.percentile(data, 20)

            x2 = np.percentile(data, 20)+(np.percentile(data, 80) - np.percentile(data, 20))/2 

            print("\nData max: ", np.array(data).max(), "\tInterquantile: ",  x0, x1, x2)

            for i in range(tree.GetEntries()):
                tree.GetEntry(i)
                cut_formula = ROOT.TTreeFormula("cut", cut, tree)
                ydwc2 = tree.GetLeaf("YDWC2").GetValue()
                tow0val = tree.GetLeaf(var0).GetValue()
                if (cut_formula.EvalInstance() and tow0val>x1):
                    HitTowProfCut.Fill(ydwc2, tow0val)
                    #if(ydwc2<0.1 and ydwc2>-1): print("YDWC2==0! TotS: ", totS, "\tYDWC2 value: ", ydwc2)            


            #enemin = np.percentile(x, 25)+(np.percentile(x, 75) - np.percentile(x, 25))/2 


            towCenterProf.SetLineColor(ROOT.kRed+1)
            towUpProf.SetLineColor(ROOT.kBlue)
            towDownProf.SetLineColor(ROOT.kGreen+2)
            towSumProf.SetLineColor(ROOT.kBlack)
            towSumProf.SetLineColor(ROOT.kBlack)
            totSProf.SetLineColor(ROOT.kAzure+1)
            towCenterProf.SetMinimum(0.)
            towCenterProf.SetMaximum(23.)
            #HitTowProfCut.SetLineColor(ROOT.kOrange+1)
            #HitTowProfCut.SetLineWidth(2)




            maxbin = towCenterProf.GetMaximumBin()
            xmax = towCenterProf.GetXaxis().GetBinCenter(maxbin)
            xlow = xmax - 8
            xhigh = xmax + 8
 
            adcvar = var0 + "_adc"
            parabola = ROOT.TF1("parabola", "pol2", xlow, xhigh)
            HitTowProfCut.Fit(parabola, "Q", "", xlow, xhigh)

            # Get the fit parameters
            a = parabola.GetParameter(0)
            b = parabola.GetParameter(1)
            c = parabola.GetParameter(2)

            # Find the x-value of the maximum point of the parabola
            if c < 0:  # Check if the parabola opens downwards
                # update xmax with parabola maximum
                xmax = -b / (2 * c)
                # update Y axis range for adc value 
                xlow = xmax - 3
                xhigh = xmax + 3

                y_max = parabola.Eval(xmax)
                print("Maximum at x =", xmax, "with y =", y_max)
            else:
                print("The parabola opens upwards; no maximum within this range.")

            yMaxHitTow.append(xmax)

            c = ROOT.TCanvas("c", "c", 700, 600)
            towCenterProf.SetTitle("Beam in tower {}".format(col[tow]))
            towCenterProf.GetXaxis().SetTitle("YDWC2")
            towCenterProf.GetYaxis().SetTitle(var0)
            towCenterProf.Draw() 
            towUpProf.Draw("same") 
            towDownProf.Draw("same") 
            towSumProf.Draw("same") 
            totSProf.Draw("same") 
            #HitTowProfCut.Draw("same")

            parabola.SetLineWidth(2)
            parabola.SetLineStyle(2)
            parabola.SetLineColor(ROOT.kMagenta)
            parabola.Draw("same")


            sumString = var0+" + "+varUp+" + "+varDown
            leg = ROOT.TLegend(0.75, 0.82, 0.96, 0.95)
            leg.AddEntry(towCenterProf, var0, "l")
            leg.AddEntry(towUpProf, varUp, "l")
            leg.AddEntry(towDownProf, varDown, "l")
            leg.AddEntry(towSumProf, sumString, "l")
            leg.AddEntry(totSProf, "totC", "l")

            enemean.append(totSProf.GetMean(2))
            towname.append(var0)

            leg.Draw()
            c.SaveAs("CerTowProf{0}.png".format(col[tow]))  
    


            adchist = ROOT.TH1F("adchist{0}".format(col[tow]), "adchist{0}".format(col[tow]), 512, 0., 4096)

            for i in range(tree.GetEntries()):
                tree.GetEntry(i)
                cut_formula = ROOT.TTreeFormula("cut", cut, tree)

                adc_val = tree.GetLeaf(adcvar).GetValue()
                ydwc2 = tree.GetLeaf("YDWC2").GetValue()
                #print("\n value: ",value)
                #value = getattr(tree, key)

                if (cut_formula.EvalInstance() and ydwc2>xlow and ydwc2<xhigh):
                    adchist.Fill(adc_val)

            adchist.Fit("gaus", "Q")
            func = adchist.GetFunction("gaus")
            mean = func.GetParameter(1)
            rms = func.GetParameter(2)
            adchist.Fit("gaus", "Q", "", mean-rms, mean+rms)
            gaus2 = adchist.GetFunction("gaus")
            adcmax = gaus2.GetParameter(1)




            c1 = ROOT.TCanvas("c1", "c1", 700, 600)
            adchist.Draw()
            c1.SaveAs("CerAdcPeak{0}.png".format(col[tow]))
            #adcmaxbin = adchist.GetMaximumBin()
            #adcmax = adchist.GetXaxis().GetBinCenter(adcmaxbin)
            print("Tower ", var0, "\tmean: ", totSProf.GetMean(2), "\tRange: ", xlow, ", ", xhigh, "\tmax position: Y=", xmax,  "\tADC value for max: ", adcmax)
            towadcpeak.append(adcmax)




            file.Close()
            #time.sleep(10)

        
    print(towname)
    print(enemean)

    towname = np.array(towname)
    enemean = np.array(enemean)
    yMaxHitTow = np.array(yMaxHitTow)

    df = pd.DataFrame({"TowerName" : towname, "MeanEnergy" : enemean, "AdcPeak" : towadcpeak, "YDWCMaxSignal" : yMaxHitTow})
    print(df)
    df.to_csv('TowerMeansCer.csv', index=False)






if __name__=="__main__":
    main()