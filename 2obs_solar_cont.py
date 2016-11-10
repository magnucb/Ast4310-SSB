from ssbfunctions import solcont_tools
import pylab as pl
from matplotlib import rc
rc('font',**{'family':'serif'})

sc = solcont_tools()

wav   = sc.wav
Fwav  = sc.Fwav
Fcont = sc.Fcont
Isurf = sc.Isurf
Icont = sc.Icont

print 'max(Ic) =', max(Icont), '*10**10 erg/(s*micrometer*ster*cm**2), at', wav[pl.where(Icont == max(Icont))][0], 'micrometers'

pl.figure()
pl.plot(wav, Fwav, label=r'F$_{\lambda}$')
pl.plot(wav, Fcont, label=r"F$_{\lambda}$'")
pl.plot(wav, Isurf, label=r'I$_{\lambda}$')
pl.plot(wav, Icont, label=r"I$_{\lambda}$'")
pl.legend()
pl.grid('on')
pl.title("Table 5: Fluxes and intensities variations")
pl.ylabel(r"10$^{10}$ erg/(s cm$^2$ $\mu$m ster)")
pl.xlabel(r"Wavelength [$\mu$m]")
pl.savefig("figs/2obs_sol_cont.png")

convfac = sc.wavtohz(wav)

wav_nu   = wav*convfac*1e10
Fwav_nu  = Fwav*convfac*1e10
Fcont_nu = Fcont*convfac*1e10
Isurf_nu = Isurf*convfac*1e10
Icont_nu = Icont*convfac*1e10

pl.figure()
pl.plot(wav, Fwav_nu*1e5, label=r'F$_{\nu}$')
pl.plot(wav, Fcont_nu*1e5, label=r"F$_{\nu}$'")
pl.plot(wav, Isurf_nu*1e5, label=r'I$_{\nu}$')
pl.plot(wav, Icont_nu*1e5, label=r"I$_{\nu}$'")
pl.legend()
pl.grid('on')
pl.title("Table 5: Fluxes and intensities variations")
pl.ylabel(r"10$^{-5}$ erg/(s cm$^2$ Hz ster)")
pl.xlabel(r"Wavelength [$\mu$m]")
pl.savefig("figs/2obs_sol_cont_freq.png")

print 'max(Ic_nu) =', max(Icont_nu), 'erg/(s*Hz*ster*cm**2), at', wav[pl.where(Icont_nu == max(Icont_nu))][0], 'micrometers'