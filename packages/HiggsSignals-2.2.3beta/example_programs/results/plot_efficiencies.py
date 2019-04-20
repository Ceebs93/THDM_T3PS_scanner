import matplotlib.pyplot as plt
import pylab as P
import numpy as np

P.rc('text', usetex=True)
P.rc('text.latex', preamble='\usepackage{amsmath}\usepackage{color}')
font = {'size' : 16}
P.rc('font', **font)
P.rc('grid', linewidth=1,color='#666666')

files = ["HS_efficiencies_kv0p9.dat",\
         "HS_efficiencies_kv1p0.dat",\
         "HS_efficiencies_kv1p1.dat"]

ls = [ ['red','solid'], ['green','dashed'],['blue','dotted']] 
labels = ['$\kappa_V =0.9$', '$\kappa_V =1.0$', '$\kappa_V =1.1$']        

fig, ax = plt.subplots()

for i,f in enumerate(files):         
	data = np.loadtxt(f)
	dataT = zip(*data)

	chi2min = min(dataT[2])
	print chi2min
# 	minindex = dataT[5].index(chi2min)

	ax.plot(dataT[1],dataT[2],color=ls[i][0],linestyle=ls[i][1], linewidth=2, label=labels[i])

ax.set_xlim([min(dataT[1]),max(dataT[1])])
ax.set_ylim([40.,100.])

# ax.text(0.02,1.015,r"$\Delta m_H^\text{theo} = "+format(dataT[1][0],'4.2f')+"~\mathrm{GeV}$, $\Lambda = "+format(dataT[7][0],'3.1f')+"$",transform=ax.transAxes)
# ax.text(0.42,1.02,r",transform=ax.transAxes)
plt.grid()

leg = plt.legend(fancybox = True, fontsize=14, loc=3)
leg.get_frame().set_alpha(0.5)

plt.xlabel(r'$\xi_V$')
plt.ylabel(r'$\chi^2$')

plt.savefig('HS_efficiencies.pdf')

