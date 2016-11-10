from ssbfunctions import falc_tools
import pylab as pl
from matplotlib import rc
rc('font',**{'family':'serif'})

falc = falc_tools()
h = falc.h
nhyd = falc.nhyd
nel = falc.nel
nprot = falc.nprot
nel_anion = nel - nprot

pl.figure()
pl.plot(h,nhyd, label=r'n$_H$')
pl.plot(h,nel, label=r'n$_e$')
pl.plot(h,nprot, label=r'n$_p$')
pl.plot(h,nel_anion, label=r'n$_{e,free}$')
pl.grid('on')
pl.legend(loc='best')
pl.yscale('log')
pl.title('Density comparisons')
pl.xlabel('Height [km]')
pl.ylabel(r'No. densities [1/cm$^3$]')
pl.savefig('figs/dens_comp.png')
