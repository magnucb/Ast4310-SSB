from ssbfunctions import falc_tools, earth_tools
import pylab as pl
from matplotlib import rc
rc('font',**{'family':'serif'})

falc = falc_tools()
earth = earth_tools()

hearth    = earth.h
hfalc     = falc.h

tearth    = earth.temp
tfalc     = falc.temp

pearth    = 10**earth.logP
pfalc     = falc.ptot

densearth = 10**earth.logrho
densfalc  = falc.dens

Nearth    = 10**earth.logN
Nfalc     = falc.nhyd + falc.nprot + falc.nel + falc.nhyd*0.1 # last is He

mcearth   = pearth/980.665 
mfalc     = falc.colm


print "Variable ratios of Earth / Sun on their respective surfaces:"
print
print "Temperature      :", tearth[hearth==min(hearth)][0]   / tfalc[pl.where(hfalc==0)][0]
print
print "Total pressure   :", pearth[hearth==min(hearth)][0]   / pfalc[pl.where(hfalc==0)][0]
print
print "Total density    :", densearth[hearth==min(hearth)][0]/ densfalc[pl.where(hfalc==0)][0]
print
print "Particle density :", Nearth[hearth==min(hearth)][0]   / Nfalc[pl.where(hfalc==0)][0]
print
print "Column density   :", mcearth[hearth==min(hearth)][0]  / mfalc[pl.where(hfalc==0)][0]
print

Rsun = 695700e5 # in cm
Dist = 1.496e+13 # 1 AU in cm

nphot_top = falc.nphot_high(5770)
Nphot_atEarth = pl.pi*nphot_top*(Rsun/Dist)**2
print "Nphot_atEarth                 =", Nphot_atEarth, "[1/cm**3]"
print 
print "Nphot_atEarth/Npart_Earth     =", Nphot_atEarth/Nearth[hearth==min(hearth)][0]
print
print "Nphot_atEarth/Npart_fromEarth =", Nphot_atEarth/falc.nphot_lower(tearth[hearth==min(hearth)][0])
print
print "Nphot_atEarth/Nphot_atSun     =", Nphot_atEarth/nphot_top 

