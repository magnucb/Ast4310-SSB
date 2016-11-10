from ssbfunctions import specline, solcont_tools
import pylab as pl
from scipy.signal import argrelextrema
from matplotlib import rc
rc('font',**{'family':'serif'})

sl = specline()

h            = sl.h
temp         = sl.temp
vturb        = sl.vturb
nhyd         = sl.nhyd
nprot        = sl.nprot
nel          = sl.nel
m_H          = sl

ionstage     = 1
level        = 2
b_l          = 1.
b_u          = 1.
A_Na         = 1.8*1e-6 # Sodium abundance
f_lu         = [0.318,0.631]
pgas         = sl.pgasptot*sl.ptot

wvn = sl.wvn # cm**-1
wls = 1e8/wvn # cm to Angstrom
wvl_air      = sl.wvl_vac2air(wls)
IntUncor     = sl.IntScaleSunUncorEarth

extremalargs = argrelextrema(IntUncor, pl.less)[0]
extremalargs = extremalargs[IntUncor[extremalargs]<0.2]
minima       = IntUncor[extremalargs]

wav = wvl_air[extremalargs]*1e-8    # wavelengths of line in cm
wvl = pl.linspace(wav[0]-2e-8, wav[0]+2e-8, 1001) #range around center

# doppler material
doppler = sl.dopplerwidth(wav[0], temp, vturb, sl.m_Na) #values for NaID1

# voigt profile
gamma        = sl.gammavdw_NaD(temp, pgas, level)    #van der waal damping
voigtprofile = pl.zeros((len(h), len(wvl)))

a = pl.zeros((len(h),len(wvl)))
v = pl.zeros((len(h),len(wvl)))
for j in pl.arange(len(h)):
    for i in pl.arange(len(wvl)):
        a[j][i] = (wvl[i]**2)*gamma[j]/(doppler[j]*4*pl.pi*sl.c)
        v[j][i] = (wvl[i]-wav[0])/doppler[j]
        voigtprofile[j][i] = sl.voigt(a[j][i],v[j][i])*doppler[j]*(pl.pi**0.5)


exthminID    = pl.zeros((len(wvl), len(h)))
extcont      = pl.zeros((len(wvl), len(h)))
exttotal     = pl.zeros((len(wvl), len(h)))
ext_elec     = sl.sigma_Thomson*sl.nel

# NaDI extinction
NaID_ext     = pl.zeros((len(wvl), len(h)))
for j in pl.arange(len(wvl)):
    for i in pl.arange(len(h)):
        fac1 = (pl.pi**0.5)*pl.exp(2)*wvl[j]*wvl[i]/ (b_l*sl.m_e*sl.c**2)
        fac2 = sl.sahabolt_Na(temp[i],nel[i],ionstage,level)
        fac3 = nhyd[i]*A_Na*f_lu[0]
        fac4 = (pl.pi**0.5)*voigtprofile[i,j]
        fac5 = 1. - (b_u/b_l)*pl.exp(-sl.hcons*sl.c/(wav[0]*sl.kerg*temp[i]))
        # print type(fac1),type(fac2),type(fac3),type(fac4),type(fac5), i, j
        # print 
        NaID_ext[j][i]  = fac1*fac2*fac3*fac4*fac5
        
        # Continuum extinction
        exthminID[j][i] = sl.exthmin(wvl[j]*1e8, temp[i], nel[i])*(nhyd[i] - nprot[i])
        extcont[j][i]   = exthminID[j][i] + ext_elec[i]
        exttotal[j][i]  = NaID_ext[j][i] + extcont[j][i]


pl.figure()
pl.plot(h,NaID_ext[pl.where(wvl==wav[0])][0],'-',label='NaID')
pl.plot(h,extcont[pl.where(wvl==wav[0])][0],'--',label='Continuum')
pl.xlim([min(h),2000])
pl.title(r'Extinction at line center')
pl.xlabel(r'Height [km]')
pl.ylabel(r'Extinction $\alpha_\lambda^l$ [cm$^{-1}$]')
pl.yscale('log')
pl.legend(loc='best')
pl.ylim([1e-14,1e-2])
# pl.show()

deltawav = (wvl - wav[0])*1e8
NaID_exttrans = pl.transpose(NaID_ext)
pl.figure()
pl.plot(deltawav,NaID_exttrans[pl.where(h==0)][0],'-.',label='NaI h=0')
# pl.plot(deltawav,extcont_swapped[pl.where(h==0)][0],'--',label='Cont h=0')
# pl.plot(deltawav,exttotal_swapped[pl.where(h==0)][0],label=r'Total h=0')
pl.plot(deltawav,NaID_exttrans[pl.where(h==560)][0],'-.',label='NaI h=560')
# pl.plot(deltawav,extcont_swapped[pl.where(h==560)][0],'--',label='Cont h=560')
# pl.plot(deltawav,exttotal_swapped[pl.where(h==560)][0],label=r'Total h=560')
pl.yscale('log')
pl.xlim([min(deltawav),max(deltawav)])
pl.title(r'Extinction')
pl.ylabel(r'Extinction $\alpha_\lambda^l$ [cm$^{-1}$]')
pl.xlabel(r'$\Delta\lambda$[$\AA$]')
pl.legend(loc='best')

### for the final contribution to falc calc data
sc = solcont_tools()
Icont   = sc.Icont
wav     = sc.wav
h       = sc.h
# IntUncor = S_spec

pl.figure()
calcfalcInt = pl.zeros(len(wav))

for i in pl.arange(len(wav)):

    hmean, contfunc, hint, intt = sc.hmean_contfunc_hint_intt(wav[i])
    calcfalcInt[i] = intt/1e14

pl.plot(wav, calcfalcInt, label="FALC calculation")
pl.plot(wav, Icont, label="Observed data")
pl.legend()
pl.grid('on')
pl.title("Observed and computed continuum intensity")
pl.xlabel(r"Wavelength [$\mu$m]")
pl.ylabel(r"Intensity [10$^{14}$ erg / (s cm$^2$ ster cm)]")
pl.savefig("figs/2disk_center_int.png")

tau5        = sl.tau5


tau         = pl.zeros(len(tau5))
calc        = pl.zeros(len(tau5))
integrand   = pl.zeros(len(tau5))
contfunc    = pl.zeros(len(tau5))


hmean, contfunc, hint, intt = hmean_contfunc_hint_intt





























pl.show()