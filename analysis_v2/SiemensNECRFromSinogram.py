# NEMA standard NECR measurement 
# Only high granularity mode for Siemens Quadra allows for the minSectorDifference cut 
# If running for Explorer or with different granularity the minSectorDifference cut should 
# be replaced by a simple DeltaPhi cut

import SiemensQuadraProperties as sqp
import CoincidenceGeneration as cg
from SimulationDataset import *
import PhysicsConstants as pc
import math
from array import array
import sys
sys.path.append('../analysis/')
from SinogramTools import *

import matplotlib.pyplot as mpl
params = {'legend.fontsize': 15,
          'legend.title_fontsize': 15,
          'legend.loc': "upper left",
          'axes.labelsize': 15,
          'xtick.labelsize': 15,
          'ytick.labelsize': 15}
mpl.rcParams.update(params)

# Fix random seed for reproducibility, or leave blank to allow the results to vary
import numpy as np
RNG = np.random.default_rng(1234)

# Set batch size for coincidence generation
# Larger batches will use more RAM, and will eventually become inefficient due to wasted events
# Smaller batches are inefficient due to overheads
# 256-1024 seems to be about right
BATCH_SIZE=1024

from ROOT import TH2F, TCanvas, TMath, TGraph, TProfile, TH1F

def IsTwoHitEvent(event) :
    if len(event) != 2 :
        return False
    return True

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

def CountRatePerformanceData(detectorMaterial, nevents, Emin, Emax, detectorLength, phantomLength) :

    tracerData = CreateDataset( detectorLength, "SiemensCrystal", phantomLength, "LinearF18", nevents, Emin, Emax, detectorMaterial, SourceOffset=45, UseNumpy=True )
    crystalData = None
    crystalActivity = None
    activityList = []
    dataList = []
    moduleIDs = tracerData.GetModuleIDs()

    # calculate crystalActivity and add it to the activityList only if the crystal material is radioactive
    if detectorMaterial == "LSO" or detectorMaterial == "LYSO" :
        crystalActivity= sqp.Lu176decaysInMass( sqp.DetectorMassLength( detectorLength, detectorMaterial ) )
        crystalData = CreateDataset( detectorLength, "SiemensCrystal", phantomLength, "Siemens", nevents, Emin, Emax, detectorMaterial, UseNumpy=True )
        activityList = [0.0, crystalActivity]
        dataList = [tracerData, crystalData]
        moduleIDs.update( crystalData.GetModuleIDs() )
    else :
        activityList = [0.0]
        dataList = [tracerData]

    return activityList, dataList, moduleIDs

def CountRatePerformance(generator, delayedGenerator, simulationWindow, PairMode, moduleIDs, Nsectors):

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
    sinogramDelayed = TH2F("sinogramDelayed", "; Projection displacement [mm]; Projection angle [rad]; Events", nbinsx, -410, 410, nbinsy, 0, 3.14)
    profileDelayed = TProfile("profileDelayed", "profileDelayed", len(binEdges)-1, binEdges)
    # deltaPhiHist = TH1F("deltaPhiHist", "; #Delta#Phi; Entries", 100, 0, np.pi)
    # deltaPhiHistTrue = TH1F("deltaPhiHistTrue", "; #Delta#Phi; Entries", 100, 0, np.pi)

    energyHist1 = TH1F("energyHist1", "First photon energy; Energy [keV]; Entries", 100, 0, 1000)
    energyHist2 = TH1F("energyHist2", "Second photon energy; Energy [keV]; Entries", 100, 0, 1000)

    if PairMode == "Exclusive":
        for event in generator: 
            if IsTwoHitEvent(event) == True:
                #only slices within the central 650 mm are used 
                zmin = -325 #mm
                zmax = 325 #mm
                z1 = event[0][DATASET_Z]
                z2 = event[1][DATASET_Z]
                zmean = (z1 + z2)/2.
                if zmean > zmax or zmean < zmin :
                    continue

                photon1GP = (event[0][DATASET_R], event[0][DATASET_PHI], event[0][DATASET_Z])
                module1ID = moduleIDs[photon1GP]
                photon2GP = (event[1][DATASET_R], event[1][DATASET_PHI], event[1][DATASET_Z])
                module2ID = moduleIDs[photon2GP]
                sectorDiff1 = module1ID[1] - module2ID[1]
                if sectorDiff1<0:
                    sectorDiff1 += Nsectors
                sectorDiff2 = module2ID[1] - module1ID[1]
                if sectorDiff2<0:
                    sectorDiff2 += Nsectors
                sectorDiff = min(sectorDiff1,sectorDiff2)
                if sectorDiff < 4:
                    continue

                # if sector information unavailable (running granularity != crystal) then do deltaPhi cut
                # deltaPhi = np.absolute(event[1][DATASET_PHI] - event[0][DATASET_PHI])
                    # if deltaPhi < 0.66 :
                    #     continue
                        
                # deltaPhiHist.Fill(deltaPhi)

                # if event[0][DATASET_EVENT] == event[1][DATASET_EVENT] :
                #     deltaPhiHistTrue.Fill(deltaPhi)

                sinogramS, sinogramTheta = CalcSinogramCoords(event)

                #fill in sinogram
                sinogram.Fill(sinogramS, sinogramTheta)
                profile.Fill(sinogramS, sinogramS)
                
                #fill in sinogram with random coincidences 
                if event[0][DATASET_EVENT] != event[1][DATASET_EVENT]:
                    sinogramRandoms.Fill(sinogramS, sinogramTheta)
                    profileRandoms.Fill(sinogramS, sinogramS)

    elif PairMode == "TakeAllGoods":
        for event in generator:
            if len(event) > 1:
                firstPhoton = event[0]
                for secondPhotonIndex in range( 1, len( event ) ):
                    #pair = np.take_along_axis( event, np.array([ 0, secondPhotonIndex ]), axis=0 )
                    pair = [ firstPhoton, event[secondPhotonIndex] ]
                    
                    # Now just repeat the "Exclusive" calculation
                    if IsTwoHitEvent(pair) == True:
                        #only slices within the central 650 mm are used 
                        zmin = -325 #mm
                        zmax = 325 #mm
                        z1 = pair[0][DATASET_Z]
                        z2 = pair[1][DATASET_Z]
                        zmean = (z1 + z2)/2.
                        if zmean > zmax or zmean < zmin :
                            continue
                        
                        photon1GP = (pair[0][DATASET_R], pair[0][DATASET_PHI], pair[0][DATASET_Z])
                        module1ID = moduleIDs[photon1GP]
                        photon2GP = (pair[1][DATASET_R], pair[1][DATASET_PHI], pair[1][DATASET_Z])
                        module2ID = moduleIDs[photon2GP]
                        sectorDiff1 = module1ID[1] - module2ID[1]
                        if sectorDiff1<0:
                            sectorDiff1 += Nsectors
                        sectorDiff2 = module2ID[1] - module1ID[1]
                        if sectorDiff2<0:
                            sectorDiff2 += Nsectors
                        sectorDiff = min(sectorDiff1,sectorDiff2)
                        if sectorDiff < 4:
                            continue



                        # if sector information unavailable (running granularity != crystal) then do deltaPhi cut
                        # deltaPhi = np.absolute(pair[1][DATASET_PHI] - pair[0][DATASET_PHI])
                        # if deltaPhi < 0.66 :
                        #     continue
                        
                        # deltaPhiHist.Fill(deltaPhi)

                        # if pair[0][DATASET_EVENT] == pair[1][DATASET_EVENT] :
                            # deltaPhiHistTrue.Fill(deltaPhi)

                        sinogramS, sinogramTheta = CalcSinogramCoords(pair)

                        #fill in sinogram
                        sinogram.Fill(sinogramS, sinogramTheta)
                        profile.Fill(sinogramS, sinogramS)

                        energyHist1.Fill(pair[0][DATASET_ENERGY])
                        energyHist2.Fill(pair[1][DATASET_ENERGY])
                        if pair[0][DATASET_ENERGY] < 435.0 or pair[0][DATASET_ENERGY] > 585.0:
                            print("First photon not falling into the energy window")
                        if pair[1][DATASET_ENERGY] < 435.0 or pair[1][DATASET_ENERGY] > 585.0:
                            print("Second photon not falling into the energy window")

                        #fill in sinogram with random coincidences 
                        if pair[0][DATASET_EVENT] != pair[1][DATASET_EVENT]:
                            sinogramRandoms.Fill(sinogramS, sinogramTheta)
                            profileRandoms.Fill(sinogramS, sinogramS)

        # Fill the sinograms using delayed coincidences
        for event in delayedGenerator:
            if len(event) > 1:
                firstPhoton = event[0]
                for secondPhotonIndex in range( 1, len( event ) ):
                    #pair = np.take_along_axis( event, np.array([ 0, secondPhotonIndex ]), axis=0 )
                    pair = [ firstPhoton, event[secondPhotonIndex] ]
                    
                    # Now just repeat the "Exclusive" calculation
                    if IsTwoHitEvent(pair) == True:
                        
                        photon1GP = (pair[0][DATASET_R], pair[0][DATASET_PHI], pair[0][DATASET_Z])
                        module1ID = moduleIDs[photon1GP]
                        photon2GP = (pair[1][DATASET_R], pair[1][DATASET_PHI], pair[1][DATASET_Z])
                        module2ID = moduleIDs[photon2GP]
                        sectorDiff1 = module1ID[1] - module2ID[1]
                        if sectorDiff1<0:
                            sectorDiff1 += Nsectors
                        sectorDiff2 = module2ID[1] - module1ID[1]
                        if sectorDiff2<0:
                            sectorDiff2 += Nsectors
                        sectorDiff = min(sectorDiff1,sectorDiff2)
                        if sectorDiff < 4:
                            continue

                        #only slices within the central 650 mm are used  
                        zmin = -325 #mm
                        zmax = 325 #mm
                        z1 = pair[0][DATASET_Z]
                        z2 = pair[1][DATASET_Z]
                        zmean = (z1 + z2)/2.
                        if zmean > zmax or zmean < zmin :
                            continue

                        sinogramS, sinogramTheta = CalcSinogramCoords(pair)

                        #fill in sinogram
                        sinogramDelayed.Fill(sinogramS, sinogramTheta)
                        profileDelayed.Fill(sinogramS, sinogramS)
                        

    else :
        print("Unrecognised coincidence pairing mode: ", PairMode)
        return
    canv = TCanvas("canv", "canv", 800, 600)

    canv.Clear()
    energyHist1.Draw("hist")
    pdfName = "Photon1Energy_" + str(activity) + ".pdf"
    canv.SaveAs(pdfName)

    canv.Clear()
    energyHist2.Draw("hist")
    pdfName = "Photon2Energy_" + str(activity) + ".pdf"
    canv.SaveAs(pdfName)


    # canv.Clear()
    # deltaPhiHist.Draw("hist")
    # pdfName = "DeltaPhi_" + str(activity) + ".pdf"
    # canv.SaveAs(pdfName)
    # canv.Clear()
    # deltaPhiHistTrue.Draw("hist")
    # pdfName = "DeltaPhiTrue_" + str(activity) + ".pdf"
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

    for b in range (1, nbinsx+1) :    
        if profileDelayed.GetBinContent(b) > 120 or profileDelayed.GetBinContent(b) < -120 :
            for by in range(0, nbinsy+1):
                sinogramDelayed.SetBinContent(b, by, 0.)            

        
    canv.Clear()
    sinogram.Draw("colz")
    print("Sinogram entries = ", sinogram.GetEntries())
    pdfName = "sinogram_unshifted_" + str(activity) + ".pdf"
    canv.SaveAs(pdfName)

    canv.Clear()
    sinogramRandoms.Draw("colz")
    pdfName = "sinogram_randoms_" + str(activity) + ".pdf"
    canv.SaveAs(pdfName)

    canv.Clear()
    sinogramDelayed.Draw("colz")
    print("Delayed sinogram entries = ", sinogramDelayed.GetEntries())
    pdfName = "sinogram_delayed_" + str(activity) + ".pdf"
    canv.SaveAs(pdfName)

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

    simulationWindowS = simulationWindow*1E-9
    #total counts and rate
    CTOT = projectionShifted.Integral()
    RTOT = CTOT/simulationWindowS

    #scatter+random counts and rate
    Csr = CalcCSR(projectionShifted)
    Rsr = Csr/simulationWindowS

    #true counts and rate 
    Ct = CTOT - Csr
    Rt = Ct/simulationWindowS

    #random counts and rate
    # Cr = sinogramRandoms.ProjectionX().Integral()
    Cr = sinogramDelayed.ProjectionX().Integral()
    Rr = Cr/simulationWindowS

    #scatter counts and rate
    Cs = Csr - Cr 
    Rs = Cs/simulationWindowS

    #NECR 
    NECR = Rt*Rt/RTOT

    return RTOT, Rsr, Rt, Rr, Rs, NECR

detectorMaterial = "LSO"
nevents = 1000000
Emin = 435.0
Emax = 585.0
detectorLength = 1024
phantomLength = 700
simulationWindow = 1E7
coincidenceWindow = 4.7
startingActivity = 1100E6 #Bq
PairMode = "TakeAllGoods"

NECRs = []
RTOTs = []
Rss = []
Rsrs = []
Rts = []
Rrs = []
activities = []
phantomRadius = 20.3 / 2.0
phantomVolume = phantomRadius * phantomRadius * math.pi * phantomLength / 10.0

activityList, dataList, moduleIDsTable = CountRatePerformanceData(detectorMaterial, nevents, Emin, Emax, detectorLength, phantomLength)

modIDs = moduleIDsTable.values()
#needed for minSectorDifference calculation
sectorMax = max({sector[1] for sector in modIDs})
nsectors = sectorMax + 1 

for t in range(0, 700, 20):
    tsec = 60*t
    activity = pc.TracerActivityAtTime( startingActivity, tsec, "F18" )
    activityList[0] = activity

    generator = cg.GenerateCoincidences( BATCH_SIZE, activityList, dataList, RNG, coincidenceWindow, simulationWindow, MultiWindow=False, EnergyResolution=0.0, EnergyMin=Emin, EnergyMax=Emax, TimeResolution=0.0 )

    delayedGenerator = cg.GenerateDelayedCoincidences( BATCH_SIZE, activityList, dataList, RNG, coincidenceWindow, simulationWindow, MultiWindow=False, EnergyResolution=0.0, EnergyMin=Emin, EnergyMax=Emax, TimeResolution=0.0, Delay=50 )

    RTOTatTime, RsrAtTime, RtAtTime, RrAtTime, RsAtTime, NECRAtTime = CountRatePerformance( generator, delayedGenerator, simulationWindow, PairMode, moduleIDsTable, nsectors)

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

# Save the numerical results to a csv files
np.set_printoptions(precision = 3, suppress = True)
results=np.vstack((NECRs, Rts, Rss, Rrs, activities))
print("results = ")
print(results)
np.savetxt("results.csv", results.T, '%1.3f', delimiter=",")

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