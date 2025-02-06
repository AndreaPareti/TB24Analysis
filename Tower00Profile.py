import pandas as pd
import numpy as np
import ROOT
import uproot


ROOT.gStyle.SetOptStat(0)



def main():
    print("Hello there")

    LeftColumn = ["1009", "1007", "1008"]
    MidLeftColumn = ["1013", "1006", "1014"]
    CenterColumn = ["1010", "1018", "1011"]
    RightColumn = ["1004", "1002", "1003"]

    column_names = ["Left", "MidLeft", "Middle", "Right"]
    column_runs = [LeftColumn, MidLeftColumn, CenterColumn, RightColumn] 
    myCut = "(abs(XDWC2 - XDWC1) < 5) & (abs(YDWC2 - YDWC1)<5) & (MCounter<200) & (TailC<300) & (C2>160) & ( (totLeakage - L20)<5000 ) & (PShower>550)"

    Rows = ["Up", "Center", "Down"]


    # variable used for profile plots
    varProf = "XDWC2"
    infolder = "/home/storage/data_apareti/TB24/T00scan/"
    treename = "Ftree" 

    t00 = ""; t10=""; t11=""; t12=""; t13=""; t14=""; t15=""; t16=""; t17=""  

    #t00 = "TC00"; t00="TC00"; t11="TC11"; t12="TC12"; t13="TC13"; t14="TC14"; t15="TC15"; t16="TC16"; t17="TC17" 
    channel = "C" # choose "S" or "C"

    totE = ""
    if (channel=="S"):
        totE = "totPMTSene"
        t00 = "TS00"; t10="TS10"; t11="TS11"; t12="TS12"; t13="TS13"; t14="TS14"; t15="TS15"; t16="TS16"; t17="TS17"  
    elif(channel=="C"): 
        totE = "totPMTCene"  
        t00 = "TC00"; t10="TC10"; t11="TC11"; t12="TC12"; t13="TC13"; t14="TC14"; t15="TC15"; t16="TC16"; t17="TC17" 


    branches = [varProf, t00, t11, t12, t13, t14, t15, t16, t17, totE]


    min_bin = 0.; max_bin=0
    if(varProf=="YDWC2"): min_bin=-30; max_bin=20
    elif(varProf=="XDWC2"): min_bin=-30; max_bin=30


    for colname, colrun in zip(column_names, column_runs):
        run_up = colrun[0]; run_center=colrun[1]; run_down=colrun[2]
        print(colname, run_up, run_center, run_down)

        filename_up = "physics_sps2024_run" + run_up + ".root"
        filename_center = "physics_sps2024_run" + run_center + ".root"
        filename_down = "physics_sps2024_run" + run_down + ".root"      

        for index, row in enumerate(Rows):
            #print(index, row, colrun[index]) 
            filename = "physics_sps2024_run" + colrun[index] + ".root"
            #filename = infolder+colrun[index]
            print(index, colname, row, filename)
            file = uproot.open(infolder+filename)
            tree = file[treename]
            data = tree.arrays(branches, cut=myCut, library="pd")
            print(branches)


            t00Prof = ROOT.TProfile("t00ProfRun{0}".format(colrun[index]), "T{1}00ProfRun{0}{2}; {3} [mm]; E [GeV]".format(colrun[index], channel, colname, varProf), 30, min_bin, max_bin)
            t11Prof = ROOT.TProfile("t11ProfRun{0}".format(colrun[index]), "T{1}11ProfRun{0}{2}; {3} [mm]; E [GeV]".format(colrun[index], channel, colname, varProf), 30, min_bin, max_bin)
            t15Prof = ROOT.TProfile("t15ProfRun{0}".format(colrun[index]), "T{1}15ProfRun{0}{2}; {3} [mm]; E [GeV]".format(colrun[index], channel, colname, varProf), 30, min_bin, max_bin)
            t16Prof = ROOT.TProfile("t16ProfRun{0}".format(colrun[index]), "T{1}16ProfRun{0}{2}; {3} [mm]; E [GeV]".format(colrun[index], channel, colname, varProf), 30, min_bin, max_bin)
            t17Prof = ROOT.TProfile("t17ProfRun{0}".format(colrun[index]), "T{1}17ProfRun{0}{2}; {3} [mm]; E [GeV]".format(colrun[index], channel, colname, varProf), 30, min_bin, max_bin)
            t10Prof = ROOT.TProfile("t10ProfRun{0}".format(colrun[index]), "T{1}10ProfRun{0}{2}; {3} [mm]; E [GeV]".format(colrun[index], channel, colname, varProf), 30, min_bin, max_bin)
            t14Prof = ROOT.TProfile("t14ProfRun{0}".format(colrun[index]), "T{1}14ProfRun{0}{2}; {3} [mm]; E [GeV]".format(colrun[index], channel, colname, varProf), 30, min_bin, max_bin)
            t13Prof = ROOT.TProfile("t13ProfRun{0}".format(colrun[index]), "T{1}13ProfRun{0}{2}; {3} [mm]; E [GeV]".format(colrun[index], channel, colname, varProf), 30, min_bin, max_bin)
            t12Prof = ROOT.TProfile("t12ProfRun{0}".format(colrun[index]), "T{1}12ProfRun{0}{2}; {3} [mm]; E [GeV]".format(colrun[index], channel, colname, varProf), 30, min_bin, max_bin)
            totEProf = ROOT.TProfile("totEProfRun{0}".format(colrun[index]), "T{1}ProfRun{0}{2}; {3} [mm]; E [GeV]".format(colrun[index], totE, colname, varProf), 30, min_bin, max_bin)


            for datarow in data.itertuples(index=False):
                t00Prof.Fill(getattr(datarow, varProf), getattr(datarow, t00))
                t11Prof.Fill(getattr(datarow, varProf), getattr(datarow, t11))
                t12Prof.Fill(getattr(datarow, varProf), getattr(datarow, t12))
                t13Prof.Fill(getattr(datarow, varProf), getattr(datarow, t13))
                t14Prof.Fill(getattr(datarow, varProf), getattr(datarow, t14))
                t15Prof.Fill(getattr(datarow, varProf), getattr(datarow, t15))
                t16Prof.Fill(getattr(datarow, varProf), getattr(datarow, t16))
                t17Prof.Fill(getattr(datarow, varProf), getattr(datarow, t17))
                totEProf.Fill(getattr(datarow, varProf), getattr(datarow, totE))


            totEProf.SetLineWidth(2); totEProf.SetLineColor(ROOT.kAzure+1); totEProf.SetMarkerStyle(53);  totEProf.SetMarkerColor(ROOT.kAzure+1);  totEProf.SetMarkerSize(2) 
            t15Prof.SetLineWidth(2); t15Prof.SetLineColor(ROOT.kBlue+1); t15Prof.SetMarkerStyle(53);  t15Prof.SetMarkerColor(ROOT.kBlue+1);  t15Prof.SetMarkerSize(2)
            t11Prof.SetLineWidth(2); t11Prof.SetLineColor(ROOT.kGreen+2); t11Prof.SetMarkerStyle(53); t11Prof.SetMarkerColor(ROOT.kGreen+2);  t11Prof.SetMarkerSize(2)
            t00Prof.SetLineWidth(2); t00Prof.SetLineColor(ROOT.kRed+1); t00Prof.SetMarkerStyle(53);   t00Prof.SetMarkerColor(ROOT.kRed+1);  t00Prof.SetMarkerSize(2)

            t16Prof.SetLineWidth(2); t16Prof.SetLineColor(ROOT.kBlue-7); t16Prof.SetMarkerStyle(54);  t16Prof.SetMarkerColor(ROOT.kBlue-7);  t16Prof.SetMarkerSize(2)
            t17Prof.SetLineWidth(2); t17Prof.SetLineColor(ROOT.kGreen-7); t17Prof.SetMarkerStyle(54); t17Prof.SetMarkerColor(ROOT.kGreen-7);  t17Prof.SetMarkerSize(2)
            t10Prof.SetLineWidth(2); t10Prof.SetLineColor(ROOT.kRed-7); t10Prof.SetMarkerStyle(54);   t10Prof.SetMarkerColor(ROOT.kRed-7);  t10Prof.SetMarkerSize(2)

            t14Prof.SetLineWidth(2); t14Prof.SetLineColor(ROOT.kBlue-7); t14Prof.SetMarkerStyle(55);  t14Prof.SetMarkerColor(ROOT.kBlue-7);  t14Prof.SetMarkerSize(2)
            t13Prof.SetLineWidth(2); t13Prof.SetLineColor(ROOT.kGreen-7); t13Prof.SetMarkerStyle(55); t13Prof.SetMarkerColor(ROOT.kGreen-7);  t13Prof.SetMarkerSize(2)
            t12Prof.SetLineWidth(2); t12Prof.SetLineColor(ROOT.kRed-7); t12Prof.SetMarkerStyle(55);   t12Prof.SetMarkerColor(ROOT.kRed-7);  t12Prof.SetMarkerSize(2)



            ctest = ROOT.TCanvas("c{0}{1}".format(row, colname), "c{0}{1}".format(row, colname), 1400, 1200)
            totEProf.SetMinimum(0.); totEProf.SetMaximum(24.)
            totEProf.SetTitle("{1} Energy Profile Run {0}".format(colrun[index], channel))
            totEProf.Draw()
            t00Prof.Draw("same")
            t15Prof.Draw("same")
            t11Prof.Draw("same")
            t10Prof.Draw("same")
            t16Prof.Draw("same")
            t17Prof.Draw("same")
            t14Prof.Draw("same")
            t13Prof.Draw("same")
            t12Prof.Draw("same")

            leg = ROOT.TLegend(0.85, 0.73, 0.995, 0.995)
            leg.SetHeader("Position: {0}, {1}".format(row, colname), "c")
            leg.SetTextSize(0.0157)
            leg.AddEntry(totEProf, totE, "PL")
            leg.AddEntry(t00Prof, t00, "PL")
            leg.AddEntry(t11Prof, t11, "PL")
            leg.AddEntry(t15Prof, t15, "PL")
            leg.AddEntry(t16Prof, t16, "PL")
            leg.AddEntry(t17Prof, t17, "PL")
            leg.AddEntry(t10Prof, t10, "PL")
            leg.AddEntry(t12Prof, t12, "PL")
            leg.AddEntry(t13Prof, t13, "PL")
            leg.AddEntry(t14Prof, t14, "PL")

            leg.Draw()

            ctest.SaveAs("T{3}00scan_{0}{1}_{2}.png".format(row, colname, varProf, channel))






if __name__=="__main__":
    main()