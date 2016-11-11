from ssbfunctions import specline, solcont_tools
import pylab as pl
from scipy.signal import argrelextrema
from matplotlib import rc
rc('font',**{'family':'serif'})

sl = specline()


wvn = sl.wvn # cm**-1
wls = 1e8/wvn # cm to Angstrom
wvl_air      = sl.wvl_vac2air(wls)
IntUncor     = sl.IntScaleSunUncorEarth

extremalargs = argrelextrema(IntUncor, pl.less)[0]
extremalargs = extremalargs[IntUncor[extremalargs]<0.2]
minima       = IntUncor[extremalargs]
wavs = wvl_air[extremalargs]


wvl_range = pl.linspace(wavs[0]-2, wavs[0]+2, 1001)
#range around center in angstrom
final_spectra = sl.Na_spectra(wavs[0], wvl_range)


pl.figure()
pl.plot(wvl_air, IntUncor, '--', label=r'Observed intensity')
pl.plot(wvl_range, final_spectra/max(final_spectra), '-', label=r'NaID$_1$ line profile')
pl.legend(loc='lower left')
pl.grid('on')
pl.xlim([wvl_range[0], wvl_range[-1]])
pl.title(r"NaID$_1$ intensity line, computed and data")
pl.xlabel(r"Wavelength [$\AA$]")
pl.ylabel(r"Relative intensity")
# ax = pl.gca()   #This disables the offset
# ax.ticklabel_format(useOffset=False)    #disable offset

pl.savefig("figs/3NaDspectrum.png")
