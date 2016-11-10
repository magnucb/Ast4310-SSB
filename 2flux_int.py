from ssbfunctions import solcont_tools
import pylab as pl
import sys
from matplotlib import rc
rc('font',**{'family':'serif'})

xgauss=[-0.7745966692,0.0000000000,0.7745966692]
wgauss=[ 0.5555555555,0.8888888888,0.5555555555]

sc          = solcont_tools()
wav         = sc.wav
tau5        = sc.tau5
nel         = sc.nel
nhyd        = sc.nhyd
nprot       = sc.nprot
h           = sc.h
temp        = sc.temp
Fcont       = sc.Fcont
fluxspec    = pl.zeros(len(wav))
integrand   = pl.zeros(len(tau5))
intmu       = pl.zeros((3, len(wav)))
intt        = 0.


for imu in pl.arange(3):
    mu = 0.5 + xgauss[imu] / 2.
    wg = wgauss[imu] / 2.

    for iw in pl.arange(len(wav)):
        wl = wav[iw]
        ext = pl.zeros(len(tau5))
        tau = pl.zeros(len(tau5))
        intt = 0.
        for i in pl.arange(1, len(tau5)):
            ext[i] = sc.exthmin(wl*1e4, temp[i], nel[i])*(nhyd[i]-nprot[i]) + sc.sigma_Thomson*nel[i]
            tau[i] = tau[i-1] + 0.5*(ext[i] + ext[i-1])*(h[i-1] - h[i])*1e5
            integrand[i] = sc.Planckfunc(temp[i], wl*1e-4)*pl.exp(-tau[i]/mu)
            intt += 0.5*(integrand[i] + integrand[i-1])*(tau[i]-tau[i-1])/mu
        
        intmu[imu,iw] = intt
        fluxspec[iw] = fluxspec[iw] + wg*intmu[imu, iw]*mu
fluxspec *= 2

pl.figure()
pl.plot(wav, fluxspec*1e-14, label='computed from FALC')
pl.plot(wav, Fcont, label='observed (Allen 1978)')
pl.legend(loc='upper right')
pl.title('Observed and computed continuum flux')
pl.xlabel(r'Wavelength [$\mu$m]')
pl.ylabel(r'Astrophysical flux [10$^{14}$ erg / (s cm$^{2}$ ster cm)]')
pl.grid(True)
pl.savefig("figs/2flux.png")