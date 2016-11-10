import pylab as pl
from scipy import special

class falc_tools:    
    def __init__(self):
        self.h, self.tau5, self.colm, self.temp, self.vturb, self.nhyd, self.nprot, self.nel, self.ptot, self.pgasptot, self.dens = pl.loadtxt('data/falc.dat', usecols=(0,1,2,3,4,5,6,7,8,9,10), unpack=True)
        self.keV    = 8.61734E-5 
        self.kerg   = 1.380658e-16 # boltzmann's constant in ergs
        self.erg2eV = 6.242e11
        self.c      = 2.99792e10   # cm /s
        self.cmm    = 2.99792e14   # mum/s
        self.hcons  = 6.62607e-27  # erg s 
        self.sigma_Thomson = 6.648e-25   # Thomson cross section [cm^2]
        self.m_e    = 9.109390E-28  #Electron mass in grams
        self.m_H    = 1.67352E-24   #Hydrogen mass in grams
        self.m_He   = 3.98*self.m_H  #Helium mass in grams
        self.m_p    = 1.67262E-24  #Proton mass in grams
        self.m_u    = 1.660538e-24  #Atomic mass unit in grams

    def nphot_lower(self, T):
        # returns in no.dens/cm^3
        return 20*(T**3)

    def nphot_high(self, T):
        # returns in no.dens/cm^3
        return 20*(T**3)/(2*pl.pi)

class earth_tools(falc_tools):
    def __init__(self):
        falc_tools.__init__(self)
        self.h, self.logP, self.temp, self.logrho, self.logN = pl.loadtxt('data/earth.dat', usecols=(0,1,2,3,4), unpack=True)


class solcont_tools(falc_tools):
    def __init__(self):
        falc_tools.__init__(self)
        self.wav, self.Fwav, self.Fcont, self.Isurf, self.Icont = pl.loadtxt('data/solspect.dat', usecols=(0,1,2,3,4), unpack=True)

    def wavtohz(self, wav):
        # converts from wavelentgh dependency to frequency dependency
        return (wav**2)/self.cmm

    def exthmin(self, wav, temp, eldens):
        # H-minus extinction, from Gray 1992
        # input:
        # wav = wavelength [Angstrom] (float or float array)
        # temp = temperature [K]
        # eldens = electron density [electrons cm-3]
        # output:
        # H-minus bf+ff extinction [cm^2 per neutral hydrogen atom]
        # assuming LTE ionization H/H-min
        # physics constants in cgs (all cm)
        k = self.kerg   # 1.380658e-16  # Boltzmann constant [erg/K]
        h = self.hcons  # 6.626076e-27  # Planck constant [erg s]
        c = self.c      # 2.997929e10   # velocity of light [cm/s]

        # other parameters
        theta   = 5040./temp
        elpress = eldens*k*temp

        # evaluate H-min bound-free per H-min ion = Gray (8.11)
        # his alpha = my sigma in NGSB/AFYC (per particle without stimulated)
        sigmabf = (1.99654 - 1.18267E-5*wav + 2.64243E-6*wav**2 - 4.40524E-10*wav**3 + 3.23992E-14*wav**4 - 1.39568E-18*wav**5 + 2.78701E-23*wav**6)
        sigmabf *= 1e-18 # cm^2 per H-min ion

        if pl.size(wav) > 1:
            sigmabf[pl.where(wav > 16444)] = 0 # H-min ionization limit at lambda=1.6444 micron
        elif (pl.size(wav) == 1):
            if wav > 16444:
                sigmabf = 0

        # convert into bound-free per neutral H atom assuming Saha = Gray p135
        # units: cm2 per neutral H atom in whatever level (whole stage)
        
        graysaha = 4.158E-10*elpress*theta**2.5*10.**(0.754*theta) # Gray (8.12)
        kappabf = sigmabf*graysaha # per neutral H atom
        kappabf = kappabf*(1.-pl.exp(-h*c/(wav*1E-8*k*temp)))

        # check Gray's Saha-Boltzmann with AFYC (edition 1999) p168
        # logratio=-0.1761-pl.log10(elpress)+pl.log10(2.)+2.5*pl.log10(temp)-theta*0.754
        # print 'Hmin/H ratio=',1/(10.**logratio) # OK, same as Gray factor SB
        # evaluate H-min free-free including stimulated emission = Gray p136

        lwav = pl.log10(wav)
        f0 = - 2.2763 - 1.6850*lwav + 0.76661*lwav**2 - 0.0533464*lwav**3
        f1 = 15.2827 - 9.2846*lwav + 1.99381*lwav**2 - 0.142631*lwav**3
        f2 = - 197.789 + 190.266*lwav - 67.9775*lwav**2 + 10.6913*lwav**3 - 0.625151*lwav**4

        ltheta = pl.log10(theta)
        kappaff = 1e-26*elpress*10**(f0 + f1*ltheta + f2*ltheta**2)
        return kappabf + kappaff    

    def Planckfunc(self, temp, wvl):
        var1 = 2.*self.hcons*(self.c**2)/(wvl**5)
        var2 = self.hcons*self.c/(self.kerg*temp*wvl)
        return var1/(pl.exp(var2)-1)

    def findTb(self, wav, Icont):
        # cgs units all around
        var1 = self.hcons*self.c / (wav*self.kerg)
        var2 = 2.*self.hcons*self.c*self.c / (Icont*wav**5) * 1e-4 # !this expletive factor!
        return var1/pl.log(var2 + 1)

    def hmean_contfunc_hint_intt(self, wl, mu=1.):
        # wl in nanometers
        ext     = pl.zeros(len(self.tau5))
        tau     = pl.zeros(len(self.tau5))
        intt    = 0.
        hint    = 0.
        integrand = pl.zeros(len(self.tau5))
        contfunc  = pl.zeros(len(self.tau5))

        for i in pl.arange(1, len(self.tau5)):
            ext[i] = self.exthmin(wl*1e4, self.temp[i], self.nel[i])*(self.nhyd[i] - self.nprot[i]) + self.sigma_Thomson*self.nel[i]
            tau[i] = tau[i-1] + 0.5*(ext[i] + ext[i-1])*(self.h[i-1] - self.h[i])*1e5

            integrand[i] = self.Planckfunc(self.temp[i], wl*1e-4)*pl.exp(-tau[i]/mu)
            intt += 0.5*(integrand[i] + integrand[i-1])*(tau[i] - tau[i-1])/mu
            hint += self.h[i]*0.5*(integrand[i] + integrand[i-1])*(tau[i] - tau[i-1])
            contfunc[i] = integrand[i]*ext[i]

        hmean = hint/intt
        return hmean, contfunc, hint, intt

class specline(solcont_tools):
    def __init__(self):
        solcont_tools.__init__(self)
        self.wvn, self.IntScaleEarth, self.IntScaleSunUncorEarth, self.IntScaleSunCorEarth = pl.loadtxt('data/int_nad.dat', usecols=(0,1,2,3), unpack=True)
        self.chiion_Na = pl.array([5.139, 47.29, 71.64])
        self.m_Na      = 22.99*1.6605*1e-24     # Sodium mass in grams
        self.wvl       = 1e8/self.wvn           # Angstrom

    def wvl_vac2air(self, wvl):
        # wvl in Angstrom
        return 0.99972683*wvl + 0.0107 - 196.25/wvl

    def partfunc_Na(self, temp):
        # partition functions Na
        # input: temp (K)
        # output: float array(3) = partition functions U1,U2,U3
        u = pl.zeros(3)
        # partition function Na I: follow Appendix D of Gray 1992
        # log(U1(T)) = c0 + c1 * log(theta) + c2 * log(theta)^2 +
        # c3 *log(theta)^3 + c4 log(theta)^4
        # with theta=5040./T
        theta = 5040./temp
        # partition function Na I : Appendix D of Gray (1992)

        c0 =  0.30955
        c1 = -0.17778
        c2 =  1.10594
        c3 = -2.42847
        c4 =  1.70721
        logU1 = c0 + c1*pl.log10(theta) + c2*pl.log10(theta)**2 + c3*pl.log10(theta)**3 + c4*pl.log10(theta)**4
        u[0]  = 10**logU1
        # partition function Na II and Na III: approximate by the
        # statistical weights of the ion ground states
        u[1]  = 1 # from Allen 1976
        u[2]  = 6 # from Allen 1976
        return u


    def saha_Na(self, temp, eldens, ionstage):
        kevT   = self.keV * temp
        kergT  = self.kerg * temp

        u = self.partfunc_Na(temp)
        u = pl.append(u, 2)
        sahaconst = (2./eldens)*(2.*pl.pi*self.m_e*kergT/(self.hcons**2))**1.5

        nstage = pl.zeros(4)
        nstage[0] = 1.
        for r in range(len(nstage)-1):
            nstage[r+1] = nstage[r]*sahaconst*(u[r+1]/u[r])*pl.exp(-self.chiion_Na[r]/kevT)

        ntotal = pl.sum(nstage)
        nstagerel = nstage/ntotal
        return nstagerel[ionstage - 1]


    def boltz_Na(self, temp, r, s):
        "Boltzmann distribution n_r,s/N_r"
        En1 = self.hcons*self.c*self.erg2eV / 5895.92e-8 
        En2 = self.hcons*self.c*self.erg2eV / 5889.95e-8 
        u   = self.partfunc_Na(temp)
        chi = [0, En1, En2]
        g   = [2, 2, 4]
        relnrs = g[s]*pl.exp(-(chi[s])/(self.keV*temp)) / u[r]
        return relnrs

    def sahabolt_Na(self, temp, eldens, ionstage, level):
        return self.saha_Na(temp, eldens, ionstage)*self.boltz_Na(temp, ionstage, level)

    def dopplerwidth(self, wav,temp,v_t,m):
        # Takes in central wavelength in cm, temperature in K, v_t in km/s, 
        # and m in grams and returns dopplerwidth in cm.
        return wav/self.c*pl.sqrt(2.*self.kerg*temp/m + v_t*v_t*1e10)

    def voigt(self, a, v):
        z = v + 1j*a
        # print "a=", a," ; v=", v
        V = special.wofz(z).real
        return V

    # voigt_NaD = voigt(a_voigt, v_voigt) / dopplerwidth

    def gammavdw_NaD(self, temp, pgas, s):
        # Van der Waals broadening for Na D1 and Na D2
        # s=2 : Na D1
        # s=3 : Na D2
        # classical Unsold recipe
        rsq_u = self.rsq_NaD(s)
        rsq_l = self.rsq_NaD(1) # lower level D1 and D2 lines are ground state s=1
        loggvdw = 6.33 + 0.4*pl.log10(rsq_u - rsq_l) + pl.log10(pgas) - 0.7*pl.log10(temp)
        return 10**loggvdw

    def rsq_NaD(self, s):
        # compute mean square radius of level s of Na D1 and Na D2 transitions
        # -> needed for van der Waals broadening in SSB
        # s=1 : ground state, angular momentum l=0
        # s=2 : Na D1 upper level l=1
        # s=3 : Na D2 upper level l=1
        h = self.hcons
        c = self.c
        erg2ev = 1/1.60219e-12 # erg to eV conversion
        E_ionization = 5.139   # [eV] ionization energy
        E_n = pl.zeros(3)      # energy level: E_n[0]=0:ground state
        E_n[1] = h*c*erg2ev/5895.94e-8 # Na D1: 2.10285 eV
        E_n[2] = h*c*erg2ev/5889.97e-8 # Na D2: 2.10498 eV
        Z = 1.                 # ionization stage, neutral Na: NaI
        Rydberg = 13.6         # [eV] Rydberg constant
        l = [0., 1., 1.]       # angular quantum number
        nstar_sq = Rydberg*(Z**2)/(E_ionization - E_n[s-1])
        rsq = nstar_sq*(5*nstar_sq + 1 - 3*l[s-1]*(l[s-1] + 1))/(2*Z**2)
        return rsq