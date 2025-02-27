import numpy as np
from SimulationDataset import *

def CalcSinogramCoords(event):
    x1 = event[0][DATASET_R]*np.cos(event[0][DATASET_PHI])
    y1 = event[0][DATASET_R]*np.sin(event[0][DATASET_PHI])
    
    x2 = event[1][DATASET_R]*np.cos(event[1][DATASET_PHI])
    y2 = event[1][DATASET_R]*np.sin(event[1][DATASET_PHI])

    sinogramTheta = np.atan2(x1 - x2, y1 - y2)

    denom = (y1 - y2) * (y1 - y2) + (x2 - x1) * (x2 - x1)
    sinogramS = 0.

    if denom != 0. :
        denom = np.sqrt(denom)
        sinogramS = (x1 * (y1 - y2) + y1 * (x2 - x1))/ denom

    if sinogramTheta < 0.0 :
        sinogramTheta = sinogramTheta + np.pi
        sinogramS = -sinogramS
    
    return sinogramS, sinogramTheta
