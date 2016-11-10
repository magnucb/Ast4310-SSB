from ssbfunctions import earth_tools
import pylab as pl
from matplotlib import rc
rc('font',**{'family':'serif'})

earth = earth_tools()
h, logP, temp, logrho, logN = earth.h, earth.logP, earth.temp, earth.logrho, earth.logN

pl.figure()
pl.plot(h, temp)
pl.grid('on')
pl.title('Temperatures in Earth atmosphere')
pl.xlabel('Height [km]')
pl.ylabel('Temperature [K]')
pl.savefig('figs/earth_temp.png')


pl.figure()
pl.plot(h, 10**logP)
pl.grid('on')
pl.title('Pressure in Earth atmosphere')
pl.xlabel('Height [km]')
pl.ylabel(r'Pressure [dyn/cm$^2$]')
pl.yscale('log')
pl.savefig('figs/earth_pressure.png')

pl.figure()
pl.plot(h, 10**logN)
pl.grid('on')
pl.title('Particle. density in Earth atmosphere')
pl.xlabel('Height [km]')
pl.ylabel(r'Particle density [1/cm$^3$]')
pl.yscale('log')
pl.savefig('figs/earth_pardens.png')

pl.figure()
pl.plot(h, 10**logrho)
pl.grid('on')
pl.title('Gas density in Earth atmosphere')
pl.xlabel('Height [km]')
pl.ylabel(r'Gas density [g/cm$^3$]')
pl.yscale('log')
pl.savefig('figs/earth_gasdens.png')
