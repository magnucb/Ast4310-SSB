from ssbfunctions import falc_tools
import pylab as pl
from matplotlib import rc
rc('font',**{'family':'serif'})

falc = falc_tools()

colm = falc.colm
h = falc.h

pl.figure()
pl.plot(h, colm)
pl.title("Column mass vs. height")
pl.xlabel("Height [km]")
pl.ylabel(r"m$_c$ [g/cm$^2$]")
pl.grid('on')
pl.savefig("figs/colm_vs_h.png")

pl.figure()
pl.plot(h, colm)
pl.title("Column mass vs. height")
pl.yscale('log')
pl.xlabel("Height [km]")
pl.ylabel(r"m$_c$ [g/cm$^2$]")
pl.grid('on')
pl.savefig("figs/colm_vs_h_log.png")