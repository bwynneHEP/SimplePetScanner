import sys
sys.path.append('../analysis/')
import SiemensQuadraProperties as sqp
from ActivityTools import *
from SimulationDataset import *
import PhysicsConstants as PC

import multiprocessing as mp
mp.set_start_method('fork')

from multiprocessing import Pool
import numpy as np
import matplotlib.pyplot as mpl
import bisect
from scipy.optimize import curve_fit
mpl.rcParams['text.usetex'] = True

params = {'legend.fontsize': 15,
          'legend.title_fontsize': 15,
          'legend.loc': "lower right",
          'axes.labelsize': 15,
          'xtick.labelsize': 15,
          'ytick.labelsize': 15}
mpl.rcParams.update(params)

def GetSliceNumber(detectorLength, nSlices, zmean):
    halfLength = detectorLength/2.
    sliceWidth = detectorLength/nSlices
    slices = np.arange(-1.*halfLength, halfLength, sliceWidth).tolist()
    return bisect.bisect(slices, zmean)-1 #returns the slice number where slices start at 0

def Fitfunction(x, a, b):
        return a*np.exp(-2*b*x)

def CalcStot(Rates, startingActivity):

    x = np.array([2.5, 5, 7.5, 10, 12.5])
    y = np.array(Rates)
    fig, ax = mpl.subplots()
    ax.scatter(x, y, label='Data', color='black')

    param_bounds = ([0.1e+06, 0.019], [1.3e+06, 0.024])
    params, cov = curve_fit(Fitfunction, x, y, bounds=param_bounds)
    a_fit, b_fit = params
    x_fit = np.linspace(2.5, 12.5, 75)
    y_fit = Fitfunction(x_fit, a_fit, b_fit)

    print("Fitted R0 = ", a_fit)
    print("Fitted mu = ", b_fit)

    startingActivitykBq = startingActivity/1000.

    Stot = a_fit/startingActivitykBq
    print("Stot = R0/Acal = ", Stot, " cps/kBq")

    mpl.plot(x_fit, y_fit, label='Fit', color='red')
    mpl.subplots_adjust(left=0.16)
    mpl.subplots_adjust(bottom=0.12)
    mpl.legend(loc='upper right')
    mpl.xlabel('$X_j$ [mm]')
    mpl.ylabel('$R_{CORR,j}$ [cps]')
    mpl.savefig("StotFit.pdf")
    mpl.clf()
    return Stot

def PlotAxialSensitivityProfile(R1i, R1, STOT, nSlices):
    R1i = np.array(R1i)
    Si = R1i*STOT*1000./R1
    sliceno = range(0, nSlices)
    mpl.plot(sliceno, Si, '.', color='black')
    mpl.xlabel("Slice number")
    mpl.ylabel("Count Rate [cps/MBq]")
    mpl.savefig("AxialSensitivity.pdf")

def OneSim(detectorLength, phantomLength, nEvents, Emin, Emax, detectorMaterial, simulationWindow, nSleeves, nSlices) :
    
    tracerData = CreateDataset( detectorLength, "SiemensCrystal", phantomLength, "LinearF18", nEvents, Emin, Emax, detectorMaterial, SourceOffset=0, NAluminiumSleeves=nSleeves, ClusterLimitMM=6.4)
    crystalData = None
    crystalActivity = None
    activityList = []
    dataList = []

    # calculate crystalActivity and add it to the activityList only if the crystal material is radioactive
    if detectorMaterial == "LSO" or detectorMaterial == "LYSO" :
        crystalActivity= sqp.Lu176decaysInMass( sqp.DetectorMassLength( detectorLength, detectorMaterial ) )
        crystalData = CreateDataset( detectorLength, "SiemensCrystal", phantomLength, "Siemens", nEvents, Emin, Emax, detectorMaterial, ClusterLimitMM=6.4 )
        activityList = [startingActivity, crystalActivity]
        dataList = [tracerData, crystalData]
    else :
        activityList = [startingActivity]
        dataList = [tracerData]

    countsPerSlice = np.zeros(nSlices)

    for event in GenerateCoincidences_new( activityList, dataList, simulationWindow, coincidenceWindow, UsePhotonTime=False, EnergyResolution=0.0, TimeResolution=0.0 ): 
        if len( event ) > 1:
            if event[0][DATASET_EVENT] != event[1][DATASET_EVENT]:
                continue
            zmean = (event[0][DATASET_Z] + event[1][DATASET_Z])/2.
            sliceNo = GetSliceNumber(detectorLength, nSlices, zmean)
            countsPerSlice[sliceNo] = countsPerSlice[sliceNo] + 1

    ratesPerSlice = countsPerSlice/simulationWindow

    RjTOT = np.sum(ratesPerSlice)

    Rcorr = simulationWindow*np.log(2)/( PC.halflives["F18"] * (1 - np.exp(-1.*simulationWindow*np.log(2)/PC.halflives["F18"])) )

    Rcorrj = Rcorr*RjTOT

    return ratesPerSlice, Rcorrj

CorrectedRates = []

#adapt as needed
nSleeves = [1, 2, 3, 4, 5]
detectorMaterial = "LSO"
nEvents = 2000000
Emin = 435.0
Emax = 585.0
simulationWindow = 0.219
coincidenceWindow = 4.7E-9
detectorLength = 1024
phantomLength = 700
startingActivity = 4.56E6 #Bq
nSlices = 640

simArgs = []
for sleeves in nSleeves:
    simArgs.append( (detectorLength, phantomLength, nEvents, Emin, Emax, detectorMaterial, simulationWindow, sleeves, nSlices) )

result = None

with Pool( processes=5 ) as p:
    result = p.starmap( OneSim, simArgs )

RatesPerSlice = []
CorrectedRates = []
for entry in result:
    RatesPerSlice.append( entry[0] )
    CorrectedRates.append( entry[1] )

#Total sensitivity
STOT = CalcStot(CorrectedRates, startingActivity)

R1 = np.sum(np.array(RatesPerSlice[0]))
PlotAxialSensitivityProfile(RatesPerSlice[0], R1, STOT, nSlices)


