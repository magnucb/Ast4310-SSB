from ssbfunctions import earth_tools
import pylab as pl
from matplotlib import rc
rc('font',**{'family':'serif'})

earth = earth_tools()
h, rho, N = earth.h, 10**earth.logrho, 10**earth.logN

mH  = 1.67352e-24 # g
mmw = rho/(N*mH)

pl.figure()
pl.plot(h, mmw)
pl.title('Mean molecular weight over height in Earth atmosphere')
pl.grid('on')
pl.xlabel('Height [km]')
pl.ylabel(r'$\mu_E$')
pl.savefig('figs/earth_meanmolweight.png')