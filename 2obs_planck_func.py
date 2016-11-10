from ssbfunctions import solcont_tools
import pylab as pl
from matplotlib import rc
rc('font',**{'family':'serif'})

sc = solcont_tools()

wav   = sc.wav
Icont = sc.Icont

T = 6450
planking = sc.Planckfunc(T, wav*1e-4)
pl.figure()
pl.plot(wav, Icont, label=r"I$_{\lambda}$'")
pl.plot(wav, planking*1e-14, label=r"B$_{\lambda}$'")
pl.legend()
pl.grid('on')
pl.title("Intensity and Planck function variations at temp. 6450 K")
pl.ylabel(r"10$^{14}$ erg/(s cm$^2$ ster)")
pl.xlabel(r"Wavelength [$\mu$m]")
pl.savefig("figs/2obs_planck.png")