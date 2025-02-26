import sys
sys.path.append('../analysis/')
import SiemensQuadraProperties as sqp
from ActivityTools import *
from SimulationDataset import *
from SinogramTools import *

import multiprocessing as mp
mp.set_start_method('fork')

from multiprocessing import Pool
import random
import numpy as np
import matplotlib.pyplot as mpl
from array import array

from ROOT import TH2F, TCanvas, TMath, TGraph, TProfile
# myColours=[royalblue, darkorange, yellowgreen, hotpink, orangered]
myColours = []
for colour in mpl.cm.viridis( np.linspace( 0.05, 0.8, 4 ) ):
    myColours.append( colour )
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=myColours)

params = {'legend.fontsize': 15,
          'legend.title_fontsize': 15,
          'legend.loc': "lower right",
          'axes.labelsize': 15,
          'xtick.labelsize': 15,
          'ytick.labelsize': 15}
mpl.rcParams.update(params)

def FindNearestBins(projectionShifted, windowEdge) :
    bin1, bin2 = 0, 0
    binContainingVal = projectionShifted.FindBin(windowEdge)
    binContainingValLowEdge = projectionShifted.GetBinLowEdge(binContainingVal)
    binWidth = projectionShifted.GetBinWidth(binContainingVal)
    distanceBinEdgeLow = np.absolute( windowEdge - binContainingValLowEdge )
    distanceBinEdgeHigh = np.absolute( windowEdge - (binContainingValLowEdge + binWidth) )
    
    if distanceBinEdgeLow < distanceBinEdgeHigh : 
        bin1 = binContainingVal-1
        bin2 = binContainingVal
    else :
        bin1 = binContainingVal
        bin2 = binContainingVal+1
    return bin1, bin2

def ConstructLinearFunction(projectionShifted, windowEdge) :
    
    nearestBin1, nearestBin2 = FindNearestBins(projectionShifted, windowEdge)
    
    x1 = projectionShifted.GetBinCenter(nearestBin1)
    y1 = projectionShifted.GetBinContent(nearestBin1)
    x2 = projectionShifted.GetBinCenter(nearestBin2)
    y2 = projectionShifted.GetBinContent(nearestBin2)
    
    slope = (y2 - y1)/(x2 - x1)
    offset = y1 - (y2-y1)*(x1)/(x2-x1)
    # print("The linear function formula is: y = ", slope, "*x + ", offset)
    return slope, offset

def GetCinter(slope, offset, windowEdge) :

    return slope*windowEdge + offset

def CalcEntriesInPartialPixel(projectionShifted, windowEdge, slope, offset):
  
    binWidth = projectionShifted.GetBinWidth(projectionShifted.FindBin(windowEdge))

    p1x = 0.0
    p2x = 0.0
    if windowEdge < 0 :
        p1x = projectionShifted.GetBinLowEdge(projectionShifted.FindBin(windowEdge))
        p2x = windowEdge
    
    else :
        p1x = windowEdge
        p2x = projectionShifted.GetBinLowEdge(projectionShifted.FindBin(windowEdge)+1)
    
    p1y = p1x*slope + offset
    p2y = p2x*slope + offset
    
    nentries = ((p1y+p2y)/2.0) * np.absolute(p2x - p1x)/binWidth

    return nentries

def CalcCSR(projectionShifted): 
    #get linear functions needed for interpolation points
    slopeLowEdge, offsetLowEdge = ConstructLinearFunction(projectionShifted, -20.0)
    slopeHighEdge, offsetHighEdge = ConstructLinearFunction(projectionShifted, 20.0)

    #R+S counts outside the (-20mm, +20mm) range
    #left side
    entriesPartialPixelLowEdge = CalcEntriesInPartialPixel(projectionShifted, -20.0, slopeLowEdge, offsetLowEdge)
    entriesOutsideStripLeft = projectionShifted.Integral(0, projectionShifted.FindBin(-20.0)-1)
    #right side
    entriesPartialPixelHighEdge = CalcEntriesInPartialPixel(projectionShifted, 20.0, slopeHighEdge, offsetHighEdge)
    entriesOutsideStripRight = projectionShifted.Integral(projectionShifted.FindBin(20.0)+1, -1)
    totCSRoutsideStrip = entriesOutsideStripLeft + entriesPartialPixelLowEdge + entriesPartialPixelHighEdge + entriesOutsideStripRight

    #R+S inside the (-20mm, +20mm) range as an average of CL and CR
    binWidth = projectionShifted.GetBinWidth(projectionShifted.FindFixBin(-20.0))
    fracOfPixInsideLeft = (projectionShifted.GetBinLowEdge(projectionShifted.FindBin(-20.0)+1) + 20.0)/binWidth
    fracOfPixInsideRight = (20.0 - projectionShifted.GetBinLowEdge(projectionShifted.FindBin(20.0)))/binWidth

    pixelsInStrip = projectionShifted.FindBin(20.0) - projectionShifted.FindBin(-20.0) -1

    CLinter = GetCinter(slopeLowEdge, offsetLowEdge, -20.0)
    CRinter = GetCinter(slopeHighEdge, offsetHighEdge, 20.0)
    totCSRinsideStrip = ((CLinter+CRinter)/2.0)*(fracOfPixInsideLeft + pixelsInStrip + fracOfPixInsideRight)

    #return total S+R counts
    return totCSRoutsideStrip + totCSRinsideStrip

def CountRatePerformance(detectorMaterial, nevents, Emin, Emax, simulationWindow, coincidenceWindow, activity, detectorLength, phantomLength):

    tracerData = CreateDataset( detectorLength, "Siemens", phantomLength, "LinearF18", nevents, Emin, Emax, detectorMaterial, SourceOffset=45)
    #ClusterLimitMM=16
    crystalData = None
    crystalActivity = None
    activityList = []
    dataList = []

    # calculate crystalActivity and add it to the activityList only if the crystal material is radioactive
    if detectorMaterial == "LSO" or detectorMaterial == "LYSO" :
        crystalActivity= sqp.Lu176decaysInMass( sqp.DetectorMassLength( detectorLength, detectorMaterial ) )
        crystalData = CreateDataset( detectorLength, "Siemens", phantomLength, "Siemens", nevents, Emin, Emax, detectorMaterial )
        activityList = [activity, crystalActivity]
        dataList = [tracerData, crystalData]
    else :
        activityList = [activity]
        dataList = [tracerData]

    nbinsx = 250
    nbinsy = 380
    xmin = -410.
    xmax = 410.
    dist = (xmax - xmin)/nbinsx

    binEdges = np.arange(-410, 410-dist, dist)
    binEdges = array('d', binEdges)

    #original unshifted sinogram 
    sinogram = TH2F("sinogram", "; Projection displacement [mm]; Projection angle [rad]; Events", nbinsx, -410, 410, nbinsy, 0, 3.14)
    #shifted sinogram needed for NECR calculation
    sinogramShifted = TH2F("sinogramShifted", "; Projection displacement [mm]; Projection angle [rad]; Events", nbinsx, -410, 410, nbinsy, 0, 3.14)
    #sinogram with random coincidences only 
    sinogramRandoms = TH2F("sinogramRandoms", "; Projection displacement [mm]; Projection angle [rad]; Events", nbinsx, -410, 410, nbinsy, 0, 3.14)
    #profiles needed for cut on SinogramS 
    profile = TProfile("profile", "profile", len(binEdges)-1, binEdges)
    profileRandoms = TProfile("profileRandoms", "profileRandoms", len(binEdges)-1, binEdges)

    for event in GenerateCoincidences_new( activityList, dataList, simulationWindow, coincidenceWindow, UsePhotonTime=False, EnergyResolution=0.0, TimeResolution=0.0 ):

        if len( event ) > 1:
            #only slices within the central 650 mm are used 
            zmin = -325 #mm
            zmax = 325 #mm
            z1 = event[0][6]
            z2 = event[1][6]
            zmean = (z1 + z2)/2.
            if zmean > zmax or zmean < zmin :
                continue

            sinogramS, sinogramTheta = CalcSinogramCoords(event)

            #fill in sinogram
            sinogram.Fill(sinogramS, sinogramTheta)
            profile.Fill(sinogramS, sinogramS)
            
            #fill in sinogram with random coincidences 
            if event[0][0] != event[1][0]:
                sinogramRandoms.Fill(sinogramS, sinogramTheta)
                profileRandoms.Fill(sinogramS, sinogramS)


    # canv = TCanvas("canv", "canv", 800, 600)
    # sinogram.Draw("colz")
    # pdfName = "sinogram_unshifted_" + str(activity) + ".pdf"
    # canv.SaveAs(pdfName)

    #apply cut on sinogramS
    for b in range (1, nbinsx+1) :
        if profile.GetBinContent(b) > 120 or profile.GetBinContent(b) < -120 :
            for by in range(0, nbinsy+1):
                sinogram.SetBinContent(b, by, 0.)

    for b in range (1, nbinsx+1) :    
        if profileRandoms.GetBinContent(b) > 120 or profileRandoms.GetBinContent(b) < -120 :
            for by in range(0, nbinsy+1):
                sinogramRandoms.SetBinContent(b, by, 0.)

    # canv.Clear()
    # sinogramRandoms.Draw("colz")
    # pdfName = "sinogram_randoms_" + str(activity) + ".pdf"
    # canv.SaveAs(pdfName)

    #zero the overflow bin
    for ybin in range(1, sinogram.GetNbinsY()+1):
        sinogram.SetBinContent(sinogram.GetNbinsX()+1, ybin, 0)

    #find the max pixel in radial distance for each projection angle bin and shift the sinogram 
    for ybin in range(1, sinogram.GetNbinsY()+1):
        shift = sinogram.ProjectionX("proj", ybin, ybin+1, "").GetMaximumBin() - (sinogram.GetNbinsX()/2)
        shift = int(shift)
        # print("shift = ", shift)
        for xbin in range(1, sinogram.GetNbinsX()+1):
            sinogramShifted.SetBinContent(xbin, ybin, sinogram.GetBinContent(xbin+shift, ybin))
            
    # canv.Clear()
    # sinogramShifted.Draw("colz")
    # pdfName = "sinogram_shifted_" + str(activity) + ".pdf"
    # canv.SaveAs(pdfName)

    # canv.Clear()
    projectionShifted = sinogramShifted.ProjectionX()
    # projectionShifted.Draw("hist")
    # pdfName = "projection_shifted_" + str(activity) + ".pdf"
    # canv.SaveAs(pdfName)

    #total counts and rate
    CTOT = projectionShifted.Integral()
    RTOT = CTOT/simulationWindow

    #scatter+random counts and rate
    Csr = CalcCSR(projectionShifted)
    Rsr = Csr/simulationWindow

    #true counts and rate 
    Ct = CTOT - Csr
    Rt = Ct/simulationWindow

    #random counts and rate
    Cr = sinogramRandoms.ProjectionX().Integral()
    Rr = Cr/simulationWindow

    #scatter counts and rate
    Cs = Csr - Cr 
    Rs = Cs/simulationWindow

    #NECR 
    NECR = Rt*Rt/RTOT

    return RTOT, Rsr, Rt, Rr, Rs, NECR

detectorMaterial = "LSO"
nevents = 1000000
Emin = 435.0
Emax = 585.0
simulationWindow = 0.001
coincidenceWindow = 4.7E-9
detectorLength = 1024
phantomLength = 700
startingActivity = 1100E6 #Bq

NECRs = []
RTOTs = []
Rss = []
Rsrs = []
Rts = []
Rrs = []
activities = []
phantomRadius = 20.3 / 2.0
phantomVolume = phantomRadius * phantomRadius * math.pi * phantomLength / 10.0

for t in range(0, 700, 20):
    tsec = 60*t
    activity = F18ActivityAtTime( startingActivity, tsec )

    RTOTatTime, RsrAtTime, RtAtTime, RrAtTime, RsAtTime, NECRAtTime = CountRatePerformance(detectorMaterial, nevents, Emin, Emax, simulationWindow, coincidenceWindow, activity, detectorLength, phantomLength)
    NECRs.append(NECRAtTime*0.000001)
    RTOTs.append(RTOTatTime*0.000001)
    Rss.append(RsAtTime*0.000001)
    Rsrs.append(RsrAtTime*0.000001)
    Rts.append(RtAtTime*0.000001)
    Rrs.append(RrAtTime*0.000001)
    activities.append(activity*0.001/phantomVolume)

#Rates plot
mpl.plot(activities, NECRs, marker='.', linestyle='-', color='red', label='NECR')
mpl.plot(activities, RTOTs, marker='.', linestyle=':', color='black', label='Prompts')
mpl.plot(activities, Rss, marker='.', linestyle='-', color='blue', label='Scatter')
mpl.plot(activities, Rrs, marker='.', linestyle='--', color='black', label='Randoms')
mpl.plot(activities, Rts, marker='.', linestyle='-', color='black', label='Trues')
mpl.xlabel( "Activity concentration [kBq/ml]")
mpl.ylabel( "Rate [Mcps]")
mpl.legend(loc='upper left')
mpl.savefig("Rates-Siemens.pdf")

mpl.clf()
#Zoomed in rates plot for comparison with siemens paper results
mpl.plot(activities, NECRs, marker='.', linestyle='-', color='red', label='NECR')
mpl.plot(activities, RTOTs, marker='.', linestyle=':', color='black', label='Prompts')
mpl.plot(activities, Rss, marker='.', linestyle='-', color='blue', label='Scatter')
mpl.plot(activities, Rrs, marker='.', linestyle='--', color='black', label='Randoms')
mpl.plot(activities, Rts, marker='.', linestyle='-', color='black', label='Trues')
mpl.xlabel( "Activity concentration [kBq/ml]")
mpl.ylabel( "Rate [Mcps]")
mpl.legend(loc='upper right')
mpl.axis([0, 40, 0, 10])
mpl.savefig("Rates-Siemens-zoom.pdf")

mpl.clf()

#Calculate the scatter fraction for different activities and plot results
activities = np.array(activities)
Rsrs = np.array(Rsrs)
Rrs = np.array(Rrs)
RTOTs = np.array(RTOTs)
num = np.subtract(Rsrs, Rrs)
denom = np.subtract(RTOTs, Rrs)
SF = np.divide(num, denom)*100.

mpl.plot(activities, SF, marker='.', linestyle='-', color='blue')
mpl.axis([0, 50, 0, 100])
mpl.xlabel( "Activity concentration [kBq/ml]")
mpl.ylabel( "Scatter Fraction [%]")
mpl.savefig("Scatter-fraction-Siemens.pdf")

#Scatter fraction measured as a mean of the scatter fractions correspodning to the last 3 measurements
SFval = np.mean(SF[-3:])
print("SFval = ", SFval)




