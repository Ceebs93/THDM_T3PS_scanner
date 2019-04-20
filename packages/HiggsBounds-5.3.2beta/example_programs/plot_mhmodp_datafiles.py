import matplotlib.pyplot as plt
import pylab as P
import numpy as np

P.rc('text', usetex=True)
P.rc('text.latex', preamble='\usepackage{amsmath}\usepackage{color}')
font = {'size' : 16}
P.rc('font', **font)
P.rc('grid', linewidth=1,color='#666666')


data = np.loadtxt("../example_data/mhmodplus/mhmod+_HiggsBounds_results.dat")

dataT = zip(*data)

fig, ax = plt.subplots()

ax.scatter(dataT[9],dataT[10], c=dataT[6],marker='s', alpha=1,edgecolors='none')
ax.tricontour(dataT[9],dataT[10],dataT[7],levels=[1.],colors=['red'],linewidths=2)


ax.set_xlim([min(dataT[9]),max(dataT[9])])
ax.set_ylim([min(dataT[10]),max(dataT[10])])


plt.xlabel(r'$M_A~[\mathrm{GeV}]$')
plt.ylabel(r'$\tan\beta$')

plt.savefig('mhmodplus.pdf')

