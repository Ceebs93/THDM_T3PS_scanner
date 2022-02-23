import matplotlib.pyplot as plt
import pylab as P
import numpy as np

P.rc('text', usetex=True)
P.rc('text.latex', preamble='\usepackage{amsmath}\usepackage{color}')
font = {'size' : 16}
P.rc('font', **font)
P.rc('grid', linewidth=1,color='#666666')

data = np.loadtxt("HSgga.dat")
dataT = zip(*data)

fig, ax = plt.subplots()

chi2min = min(dataT[3])
minindex = dataT[3].index(chi2min)

# 1/2sigma filled regions:
ax.tricontourf(dataT[2],dataT[1],map(lambda d: d-chi2min, dataT[3]),levels=[0,2.3,5.99],colors=['green','yellow','none'])
# 1/2sigma contours:
ax.tricontour(dataT[2],dataT[1],map(lambda d: d-chi2min, dataT[3]),levels=[2.3,5.99],colors=['k','k'],linewidths=2)
# SM point:
ax.plot(1.,1.,'+',markersize=10,color='k') 
# BF point:
ax.plot(dataT[2][minindex],dataT[1][minindex],'*',markersize=10,color='w') 
# H->gammagamma rate contours:
Hgaga = ax.tricontour(dataT[2],dataT[1],dataT[6],levels=[0.6,0.8,1,1.2,1.4],colors=['b'],linewidths=1)
ax.clabel(Hgaga, fontsize=8, inline=1)

ax.set_xlim([min(dataT[2]),max(dataT[2])])
ax.set_ylim([min(dataT[1]),max(dataT[1])])

plt.grid()

plt.xlabel(r'$\kappa_\gamma$')
plt.ylabel(r'$\kappa_g$')

plt.savefig('HSgga.pdf')

