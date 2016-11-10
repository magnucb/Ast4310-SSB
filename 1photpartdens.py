from ssbfunctions import falc_tools
import pylab as pl
from matplotlib import rc
rc('font',**{'family':'serif'})

fl = falc_tools()

h = fl.h
t = fl.temp
nhyd = fl.nhyd

nphot_l = fl.nphot_lower(t[h==min(h)][0])
print "nphot_lower          =", nphot_l, "[1/cm**3]"

nphotnhyd_l = nphot_l/nhyd[h==min(h)][0]
print "nphot_low/nhyd_low   = %g" % nphotnhyd_l


nphot_h = fl.nphot_high(5770)
print "nphot_higher         =", nphot_h, "[1/cm**3]"

nphotnhyd_h = nphot_h/nhyd[h==max(h)][0]
print "nphot_high/nhyd_high = %g" % nphotnhyd_h

