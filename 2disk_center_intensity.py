from ssbfunctions import solcont_tools
import pylab as pl
import sys
from matplotlib import rc
rc('font',**{'family':'serif'})

sc = solcont_tools()
Icont   = sc.Icont
wav     = sc.wav
h       = sc.h

calcfalcInt = pl.zeros(len(wav))

for i in pl.arange(len(wav)):

    hmean, contfunc, hint, intt = sc.hmean_contfunc_hint_intt(wav[i])
    calcfalcInt[i] = intt/1e14

pl.plot(wav, calcfalcInt, label="FALC calculation")
pl.plot(wav, Icont, label="Observed data")
pl.legend()
pl.grid('on')
pl.title("Observed and computed continuum intensity")
pl.xlabel(r"Wavelength [$\mu$m]")
pl.ylabel(r"Intensity [10$^{14}$ erg / (s cm$^2$ ster cm)]")
pl.savefig("figs/2disk_center_int.png")