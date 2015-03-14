"""Need critical Rayleigh # (Ra_c) for the case where beta=5e5 but
C=0.

Use Brummell & Hart (GApFD 1993) formula:

Ra_c = k**6/a**2 + beta**2/(4 k**2)

with k**2 = n**2 * pi**2 + a**2

where state variables are A \propto exp(i*a*x) sin(n*pi*y).

We set n = 1.

"""
import pylab as P
import numpy as np

n = 1
a = np.linspace(0,100,1000)
beta = 5e5

k2 = n**2 * np.pi**2 + a**2
Rac = k2**3/a**2 + beta**2/(4*k2)

Rac_min = np.min(Rac)
P.plot(a, Rac, 'k')
P.axhline(Rac_min,color='k',alpha=0.4)
P.text(5,0.31e8,r'$\mathrm{Ra_c} = %5.3e$' % Rac_min,fontsize=18)
P.ylim(0,1e8)
P.xlabel(r"$\alpha$",fontsize=18)
P.ylabel(r"$\mathrm{Ra_c}$",fontsize=18)
P.savefig("../figs/rac_linear_zeroC.png")
