from ssbfunctions import falc_tools
import pylab as pl
from matplotlib import rc
rc('font',**{'family':'serif'})

falc = falc_tools()

kerg = falc.kerg
temp = falc.temp
h = falc.h
nhyd = falc.nhyd
nel = falc.nel
nHe = 0.1*nhyd

ptot = falc.ptot
pgasptot = falc.pgasptot
pgas = pgasptot*ptot

label1 = r'P$_{gas}$'
label2 = r'(n$_{H}$+n$_{e}$)k$_{B}$T'
pl.figure()
pl.plot(h, pgas, label=label1)#'Pgas')
pl.plot(h, (nhyd+nel)*kerg*temp, label=label2)#'(nH + nel)*kT')
pl.legend(loc='best')
pl.grid('on')
pl.yscale('log')
pl.title('Gas densities')
pl.xlabel('Height [km]')
pl.ylabel(r'Pressure [dyne/cm$^2$]')
pl.savefig('figs/gasP_vs_h.png')


label3 = r'P$_{gas}$/(n$_{H}$+n$_{e}$)k$_{B}$T'
label4 = r'P$_{gas}$/(n$_{H}$+n$_{e}$+n$_{He}$)k$_{B}$T'
pl.figure()
pl.plot(h, pgas/((nhyd+nel)*kerg*temp), label=label3)#'Pgas/(nhyd + nel)*kT')
pl.plot(h, pgas/((nhyd+nel+nHe)*kerg*temp), label=label4)#'Pgas/(nhyd + nel + nHe)*kT')
pl.grid('on')
pl.legend(loc='best')
pl.title('Gas density proportions')
pl.xlabel('Height [km]')
pl.ylabel('Pressure / Ideal gas law pressure')
pl.savefig('figs/gasP_vs_h_ideal.png')