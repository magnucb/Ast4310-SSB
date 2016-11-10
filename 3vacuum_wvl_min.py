from ssbfunctions import specline, falc_tools
import pylab as pl
from scipy.signal import argrelextrema
from matplotlib import rc
rc('font',**{'family':'serif'})

sl = specline()

wvn = sl.wvn # cm**-1
wvl = 1e8/wvn # cm to Angstrom
wvl_air = sl.wvl_vac2air(wvl)

IntUncor = sl.IntScaleSunUncorEarth

extremalargs = argrelextrema(IntUncor, pl.less)[0]
extremalargs = extremalargs[IntUncor[extremalargs]<0.2]
minima = IntUncor[extremalargs]

pl.figure()
pl.plot(wvn, IntUncor, label=r"Solar intensity uncorrected")
for i in range(len(minima)):
    pl.plot([wvn[extremalargs[i]]], [minima[i]], marker='o', label=r"Na I D$_%g$, wvn=%g [cm$^-1$], $\lambda_{air}$=%g $\AA$" % (i+1, wvn[extremalargs[i]], wvl_air[extremalargs[i]]))

pl.legend(loc='upper left')
pl.grid('on')
pl.xlim([wvn[0], wvn[-1]])
pl.ylim([0.,1.35])
pl.title("Solar disk intensity atlas")
pl.xlabel(r"Wavenumber [cm$^-1$]")
pl.ylabel(r"Relative intensity")
pl.savefig("figs/3vac_wvl_min.png")
