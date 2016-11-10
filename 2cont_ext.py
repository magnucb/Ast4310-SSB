from ssbfunctions import solcont_tools
import pylab as pl
import sys
from matplotlib import rc
rc('font',**{'family':'serif'})

sc = solcont_tools()

temp  = sc.temp[pl.where(sc.h==0)][0]
eldens= sc.nel[pl.where(sc.h==0)][0]
wav   = sc.wav

# exthmin: wav in Angstrom, temp in K, eldens in electrons/cm^3
wav_A = wav*1e4 # wavelengths in Angstrom, from micrometers
Hmin_ext = sc.exthmin(wav_A, temp, eldens)

pl.figure()
pl.plot(wav, Hmin_ext*1e24)
pl.title(r"H$^-$ extinction in FALC at h = 0 km")
pl.xlabel(r"Wavelength $\lambda$ [$\mu$m]")
pl.ylabel(r"H$^-$ extinction [ 10$^{-24}$ cm$^2$/H-atom]")
pl.xlim([0, 2.0])
pl.ylim([0, 2.5])
pl.grid('on')
pl.savefig("figs/2cont_ext_gray.png")


Hmin_ext_mani = sc.exthmin(wav_A, temp, eldens)
pl.figure()
pl.plot(wav, (Hmin_ext*1e24)**-1)
pl.title(r"H$^-$ extinction in FALC at h = 0 km")
pl.xlabel(r"Wavelength $\lambda$ [$\mu$m]")
pl.ylabel(r"(H$^-$ extinction [ 10$^{-24}$ cm$^2$/H-atom])$^{-1}$")
pl.grid('on')
pl.savefig("figs/2cont_ext_Tb.png")