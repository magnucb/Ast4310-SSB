from ssbfunctions import falc_tools
import pylab as pl
import sys
from matplotlib import rc
rc('font',**{'family':'serif'})

falc = falc_tools()
ptot = falc.ptot
colm = falc.colm

pl.figure()
pl.plot(colm, ptot)
pl.xscale('log')
pl.yscale('log')
pl.title("Total pressure against column mass, log scaled")
pl.xlabel(r"Column mass [g cm$^{-2}$]")
pl.ylabel(r"Total pressure [dyn cm$^{-2}$]")
pl.grid('on')
pl.savefig("figs/totpress_vs_m_log.png")

pl.figure()
pl.plot(colm, ptot/1e4)
pl.title("Total pressure against column mass")
pl.xlabel(r"Column mass [g cm$^{-2}$]")
pl.ylabel(r"Total pressure [10$^4$ dyn cm$^{-2}$]")
pl.grid('on')
pl.savefig("figs/totpress_vs_m.png")