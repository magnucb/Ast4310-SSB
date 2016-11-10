from ssbfunctions import solcont_tools
import pylab as pl
import sys
from matplotlib import rc
rc('font',**{'family':'serif'})

sc = solcont_tools()

h       = sc.h
tau5    = sc.tau5
temp    = sc.temp
nel     = sc.nel
nhyd    = sc.nhyd
nprot   = sc.nprot
nHI     = nhyd - nprot

### 2.3
wl     = 500 # nm
wl_A   = wl*10 # Angstrom
tau = pl.zeros(len(tau5)) # init of tau array
ext = sc.exthmin(wl_A, temp, nel)

for i in pl.arange(1, len(tau)):
    tau[i] = tau[i-1] + (0.5*nHI[i]*(ext[i] + ext[i-1]) + nel[i]*sc.sigma_Thomson)*(h[i-1]-h[i])*1e5

pl.figure()
pl.plot(h, tau, label=r'$\tau_{500}$ Calculational')
pl.plot(h, tau5, '--', label=r'$\tau_{500}$ FALC')
pl.xlabel("Height [km]")
pl.ylabel(r"Optical depth $\tau$")
pl.title(r"Optical dept over height for $\lambda$=500 [nm]")
pl.yscale('log')
pl.legend() 
pl.grid('on')
pl.savefig("figs/2opt_depth.png")

### 2.4

Icont   = sc.Icont
wav     = sc.wav
wl      = 0.5

hmean, contfunc, hint, intt = sc.hmean_contfunc_hint_intt(wl)
print 'computed continuum intensity wl = %g : %g erg s^-1 cm^-2 ster^-1 cm^-1' % (wl, intt)

w = pl.where(wav==wl)
print 'observed continuum intensity wl = %g : %g erg s^-1 cm^-2 ster^-1 cm^-1' % (wav[w], Icont[w]*1e10*1e4)

pl.figure()
pl.plot(h, contfunc/max(contfunc), 'g')
pl.axvline(x=hmean)
pl.legend(["Contribution function", "Mean height of formation"])
pl.xlim([-100, 500])
pl.xlabel("Height [km]")
pl.ylabel(r"Contribution function")
pl.title(r"Peak normalized contribution function for $\lambda$=500 [nm]")
pl.grid('on')
pl.savefig("figs/2emer_peaknormcont_vs_h.png")

hmean_list = [hmean]

### more wl's

pl.figure()
cls = ['b', 'g', 'r', 'c']
pl.plot(h, contfunc/max(contfunc), label=r"$\lambda$=%g nm, h$_{mean}$=%g" % (wl*1e3, hmean), color=cls[0]) # keeping the old one
pl.axvline(x=hmean, color=cls[0], linestyle=':', linewidth=4)

wls = [1, 1.6, 5]
for i in pl.arange(len(wls)):
    hmean, contfunc, hint, intt = sc.hmean_contfunc_hint_intt(wls[i])
    pl.plot(h, contfunc/max(contfunc), label=r"$\lambda$=%g nm, h$_{mean}$=%g" % (wls[i]*1e3, hmean), color=cls[i+1])
    pl.axvline(x=hmean, color=cls[i+1], linestyle=':', linewidth=4)
    hmean_list.append(hmean)

pl.legend()
pl.xlim([-100, 500])
pl.ylim([0, 1.05])
pl.xlabel("Height [km]")
pl.ylabel(r"Contribution functions")
pl.title(r"Peak normalized contribution functions for different $\lambda$")
pl.grid('on')
pl.savefig("figs/2emer_peaknormcont_vs_h_wls.png")


### LTE and E-B approx

