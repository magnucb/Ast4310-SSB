from ssbfunctions import solcont_tools
import pylab as pl
import sys
from matplotlib import rc
rc('font',**{'family':'serif'})

sc = solcont_tools()

nhyd  = sc.nhyd
nprot = sc.nprot
nHI   = nhyd - nprot
h     = sc.h
temp  = sc.temp
eldens= sc.nel
wav   = 0.5 #micrometers
# exthmin: wav in Angstrom, temp in K, eldens in electrons/cm^3
wav_A = wav*1e4 # wavelengths in Angstrom, from micrometers
Hmin_ext = sc.exthmin(wav_A, temp, eldens)
alpha_   = Hmin_ext*nHI
alpha_e  = eldens*sc.sigma_Thomson
alpha_T  = alpha_ + alpha_e

pl.figure()
pl.plot(h, alpha_, label=r'$\alpha$(H$^-$)')
pl.yscale('log')
pl.title(r"$\alpha$(H$^-$) extinctions for $\lambda$=500 over height")
pl.xlabel(r"Height $h$ [km]")
pl.ylabel(r"Extinctions [ cm$^2$/(cm$^3$)]") # 10$^{-24}$
pl.grid('on')
pl.savefig("figs/2cont_ext_over_h.png")

pl.plot(h, alpha_e, label=r'$\alpha$($e$)')
pl.plot(h, alpha_T, label=r'$\alpha$(H$^-$)+$\alpha$($e$)')
pl.legend()
pl.savefig("figs/2cont_ext_over_h_thom.png")