from ssbfunctions import earth_tools
import pylab as pl
from matplotlib import rc
rc('font',**{'family':'serif'})

earth = earth_tools()

h, rho = earth.h, 10**earth.logrho

rho0 = rho[0]
h1 = 120 # km
Hrho = h1/(pl.log(rho0/rho[pl.where(h == h1)]))
rho_h = rho0*pl.exp(-h/(Hrho))

pl.figure()
pl.plot(h, rho, label='Data')
pl.plot(h, rho_h, label='Density scale height approx.')
pl.legend()
pl.title('Observed density and lower atmosphere density scale height approximation')
pl.ylim([1e-14, 1e-3])
pl.yscale('log')
pl.xlabel('Height [km]')
pl.ylabel(r'$\rho$ [g/cm$^3$]')
pl.grid('on')
pl.savefig('figs/earth_densityscales.png')