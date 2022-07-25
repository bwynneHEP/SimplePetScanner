import math

crystalCM3 = 0.32*0.32*2.0
crystalGram = crystalCM3*7.4
LSOmoleGram = 2.0*175.0 + 28.1 + 5.0*16.0 #Lu 2 si 1 O 5
NAvogadro = 6.022E23
LSOunits = NAvogadro * crystalGram / LSOmoleGram
Lu176atoms = LSOunits * 2.0 * 0.026 #2.6% Lu176
Lu176halfYears = 3.76E10
Lu176halfSec = Lu176halfYears * 3.154E7

#https://scipython.com/book2/chapter-6-numpy/examples/simulating-radioactive-decay/
#The probability that a given nucleus will decay in time dt is dt/tau
#tau = tHalf/ln2

Lu176crystalDecaysPerSec = Lu176atoms * math.log(2.0) / Lu176halfSec
print( Lu176crystalDecaysPerSec )
