from ssbfunctions import falc_tools
import pylab as pl
from matplotlib import rc
rc('font',**{'family':'serif'})

falc = falc_tools()
h = falc.h
nhyd = falc.nhyd
nprot = falc.nprot
ionfracH = nprot/(nhyd+nprot)

pl.figure()
pl.plot(h, ionfracH)
pl.grid('on')
pl.yscale('log')
pl.title('Ionization fraction over height')
pl.xlabel('Height [km]')
pl.ylabel(r'n$_{p}$/(n$_{p}$ + n$_{H}$)')
pl.savefig('figs/ion_frac.png')