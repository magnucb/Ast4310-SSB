from ssbfunctions import falc_tools
import pylab as pl
import sys
from matplotlib import rc
rc('font',**{'family':'serif'})

falc = falc_tools()

h = falc.h
temp = falc.temp

fig = pl.figure()
pl.plot(h, temp)
pl.title('Temperature stratification')
pl.xlabel('Height [km]')
pl.ylabel('Temperature [K]')
pl.ylim([2e3,1e4])
pl.grid('on')
pl.savefig('figs/temp_vs_h.png')