# Scintillator densities [g/cm^3]
densities = {
    "BGO" : 7.13,
    "NaI" : 3.667,
    "CsF" : 4.115,
    "LYSO" : 7.1,
    "LSO" : 7.4,
    "BaF2" : 4.89,
    "CsI" : 4.51,
    "CdWO4" : 7.9,
    "CaF2" : 3.18,
    "CZT" : 5.8,
    "FAPbI3" : 4.0,
    "CsPbBr3" : 5.8,
    "MAPbI3" : 4.0,
    "MAPbBr3" : 3.83,
    "Cs2AgBiBr6" : 4.65,
    "CaTiO3" : 0.986,
    "eCsI" : 4.51,
}

# Tracer half-lives [s]
halflives = {
    "F18" : 109.77*60.0,
    "Zr89" : 78.41*60.0*60.0,
    "C11" : 20.4*60.0,
    "O15" : 2*60.0,
    "N13" : 10*60.0,
    "Rb82" : 75.0,
    "Ga68" : 68*60.0,
}

def ActivityAtTime( StartingActivity, TimeElapsed, HalfLife ):
    return StartingActivity * ( 2.0 ** ( -TimeElapsed / HalfLife ) )

def TracerActivityAtTime( StartingActivity, TimeElapsed, Isotope ):
    return ActivityAtTime( StartingActivity, TimeElapsed, halflives[Isotope] )
