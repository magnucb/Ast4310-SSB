import pylab as pl

def parfunc_Na(self, temp):
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

def voigt(self, gamma, x):
    z = x + 1j*gamma
    V = special.wofz(z).real()
    return V

# voigt_NaD = voigt(a_voigt, v_voigt) / dopplerwidth

def gammavdw_NaD(self, temp, pgas, s):
    # Van der Waals broadening for Na D1 and Na D2
    # s=2 : Na D1
    # s=3 : Na D2
    # classical Unsold recipe
    rsq_u = rsq_NaD(s)
    rsq_l = rsq_NaD(1) # lower level D1 and D2 lines are ground state s=1
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