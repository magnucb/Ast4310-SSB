from ssbfunctions import solcont_tools
import astropy.constants as ac
import astropy.units as au
import pylab as pl
from matplotlib import rc
rc('font',**{'family':'serif'})

sc = solcont_tools()

wav   = sc.wav # micrometer, not cm
Icont = sc.Icont
Tb = sc.findTb(wav*1e-4, Icont*1e10)

pl.figure()
pl.plot(wav, Tb, label=r"T$_{b}$'")
pl.grid('on')
pl.title("Brightness temperature over wavelength")
pl.ylabel(r"T$_b$ [K]")
pl.xlabel(r"Wavelength $\lambda$ [$\mu$m]")
pl.savefig("figs/2obs_Tb.png")