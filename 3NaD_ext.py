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
wavs = wvl_air[extremalargs]*1e-8    # wavelengths of line in cm

wvl_range = pl.linspace(wavs[0]-2e-8, wavs[0]+2e-8, 1001) #range around center
final_spectra = sl.Na_spectra(wavs[0], wvl_range)


pl.figure()
pl.plot(wvl_air, IntUncor, '--')
pl.plot(wvl_range*1e8, final_spectra/max(final_spectra), '-')
# pl.legend(loc='upper left')
pl.grid('on')
pl.xlim([wvl_range[0], wvl_range[-1]])
pl.title("Solar disk intensity atlas")
pl.xlabel(r"Wavelength [cm]")
pl.ylabel(r"Relative intensity")
# ax = pl.gca()   #This disables the offset
# ax.ticklabel_format(useOffset=False)    #disable offset
pl.savefig("figs/3NaDspectrum.png")

