from ssbfunctions import earth_tools
import pylab as pl
from matplotlib import rc
rc('font',**{'family':'serif'})

earth = earth_tools()

h, P, rho, N = earth.h, 10**earth.logP, 10**earth.logrho, 10**earth.logN

maxP = max(P)
maxrho = max(rho)
maxN = max(N)

Pnf = P/maxP
rhonf = rho/maxrho
Nnf = N/maxN

pl.figure()
pl.plot(h, Pnf, label=r'P/P$_{max}$')
pl.plot(h, rhonf, label=r'$\rho$/$\rho_{max}$')
pl.plot(h, Nnf, label=r'N/N$_{max}$')
pl.grid('on')
pl.yscale('log')
pl.legend(loc='best')
pl.title('Normalised pressure, gas- and pressure densities')
pl.xlabel('Height [km]')
pl.ylabel('Ratios')
pl.savefig('figs/earth_pres_dens_vs_h.png')