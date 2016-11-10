from ssbfunctions import falc_tools
import pylab as pl
from matplotlib import rc
rc('font',**{'family':'serif'})

falc = falc_tools()

h = falc.h
nhyd = falc.nhyd
mH = 1.67352e-24 # g
denshyd = nhyd*mH
dens = falc.dens
densHratio = denshyd/dens

pl.figure()
pl.plot(h, densHratio)
pl.ylim([0.71,0.72])
pl.title("Hydrogen density/total density vs. height")
pl.xlabel("Height [km]")
pl.ylabel(r"$\rho_{H}/\rho$")
pl.grid('on')
pl.savefig("figs/ratio_Hmassdens_of_dens_vs_h.png")

nHe = 0.1*nhyd
mHe = 3.97*mH
densHe = nHe*mHe
densHHeratio = (denshyd + densHe)/dens

# pl.xkcd()
pl.figure()
pl.plot(h, densHHeratio)
pl.ylim([0.989,1])
pl.title("H+He density/total density vs. height")
pl.xlabel("Height [km]")
pl.ylabel(r"$(\rho_{H} + \rho_{He})/\rho$")
pl.grid('on')
# pl.savefig("latex/frontpage.png")
pl.savefig("figs/ratio_HHemassdens_of_dens_vs_h.png")
pl.show()