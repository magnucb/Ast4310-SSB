from ssbfunctions import solcont_tools
import pylab as pl
import sys
from matplotlib import rc
rc('font',**{'family':'serif'})

sc      = solcont_tools()
Icont   = sc.Icont
wav     = sc.wav

basefalcInt = pl.zeros(len(wav))
for i in pl.arange(len(wav)):
    hmean, contfunc, hint, intt = sc.hmean_contfunc_hint_intt(wav[i])
    basefalcInt[i] = intt/1e14

mus     = [0.1, 0.3, 0.5, 0.7, 0.9, 1.] # pl.linspace(0.1, 1., 6)
calcfalcInt = pl.zeros((len(mus), len(wav)))
pl.figure()
for j in pl.arange(len(mus)):
    for i in pl.arange(len(wav)):
        hmean, contfunc, hint, intt = sc.hmean_contfunc_hint_intt(wav[i], mus[j])
        calcfalcInt[j,i] = intt/1e14
    # pl.plot(wav, calcfalcInt[j]/basefalcInt, label=r"I$_\lambda$ at $\mu$=%g" % mus[j])
    pl.plot(wav, calcfalcInt[j], label=r"I$_\lambda$ at $\mu$=%g" % mus[j])

pl.legend()
pl.grid('on')
pl.title("Intensity variation over limb")
pl.xlabel(r"Wavelength [$\mu$m]")
pl.ylabel(r"I$_{\lambda}$(0,$\mu$) [10$^{14}$ erg / (s cm$^2$ ster cm)]")
# pl.ylabel(r"I$_{\lambda}$(0,$\mu$)/I$_{\lambda}$(0, 1)") # [10$^{14}$ erg / (s cm$^2$ ster cm)]")
pl.savefig("figs/2disk_limbdark.png")
# pl.savefig("figs/2disk_limbdark_ratio.png")

"""
wls = [0.2, 0.5, 1.6, 2.5, 5.0]
pl.figure()
for i in pl.arange(len(wls)):
    hmean, contfunc, hint, intt = sc.hmean_contfunc_hint_intt(wls[i], mus[j])
    calcfalcInt[j,i] = intt/1e14 

pl.legend()
pl.grid('on')
pl.title("Intensity variation over limb")
pl.xlabel(r"Wavelength [$\mu$m]")
pl.ylabel(r"Limb intensity [10$^{14}$ erg / (s cm$^2$ ster cm)]")
pl.savefig("figs/2disk_limbdark.png")
"""