from ssbfunctions import falc_tools
import pylab as pl
from matplotlib import rc
rc('font',**{'family':'serif'})

falc = falc_tools()

h = falc.h
dens = falc.dens

h1 = 2017 # 2 Mm, but h is defined in km
rho0 = dens[pl.where(h == 0.)]
Hrho = h1/pl.log(rho0/dens[pl.where(h==h1)])
rho = rho0*pl.exp(-h/Hrho)

pl.figure()
pl.plot(h, dens, label='Data')
pl.plot(h, rho, label='Exp. Aprx.')
pl.grid('on')
pl.yscale('log')
pl.xlabel('Height [km]')
pl.ylabel(r'Density [g/cm$^3$]')
pl.title('Gas density against height')
pl.legend(loc='best')
pl.savefig('figs/gasdens_vs_h.png')