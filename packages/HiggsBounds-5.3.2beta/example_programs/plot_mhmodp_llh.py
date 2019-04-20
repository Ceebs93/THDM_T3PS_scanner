import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pylab as P
import numpy as np

P.rc('text', usetex=True)
P.rc('text.latex', preamble='\usepackage{amsmath}\usepackage{color}')
font = {'size' : 16}
P.rc('font', **font)
P.rc('grid', linewidth=1,color='#666666')

plotsize = (7,4.8)

data = np.loadtxt("mhmodp_HBwithLHClikelihood.dat")
dataT = zip(*data)

# ATLobs = np.loadtxt("auxiliary_data/ATL_observed.dat")
# ATLobsT = zip(*ATLobs)
# 
# ATLexp = np.loadtxt("auxiliary_data/ATL_expected.dat")
# ATLexpT = zip(*ATLexp)

CMSobs = np.loadtxt("auxiliary_data/CMS_17020_observed.dat")
CMSobsT = zip(*CMSobs)

CMSexp = np.loadtxt("auxiliary_data/CMS_17020_expected.dat")
CMSexpT = zip(*CMSexp)



fig, ax = plt.subplots(figsize=plotsize)


#----- ATLAS observed exclusion --> find the minimum of -2lnL -----#
# 
# ATL_chi2min = min(dataT[21])
# ATL_minindex = dataT[21].index(ATL_chi2min) 
# 
# print "Best fit point: ", ATL_chi2min, ATL_minindex, dataT[3][ATL_minindex], dataT[4][ATL_minindex]
# 
#-----  ATLAS expected exclusion --> find the minimum of -2lnL ------#
#
# ATL_chi2min = min(dataT[22])
# ATL_minindex = dataT[22].index(ATL_chi2min) 
# 
# print "Best fit point: ", ATL_chi2min, ATL_minindex, dataT[3][ATL_minindex], dataT[4][ATL_minindex]
#
#----- CMS observed exclusion --> find the minimum of -2lnL -----#

# CMS_chi2min = min(dataT[19])
# CMS_minindex = dataT[19].index(CMS_chi2min) 
# 
# print "Best fit point: ", CMS_chi2min, CMS_minindex, dataT[3][CMS_minindex], dataT[4][CMS_minindex]

CMS_chi2min = min(dataT[17])
CMS_minindex = dataT[17].index(CMS_chi2min) 

print "Best fit point: ", CMS_chi2min, CMS_minindex, dataT[3][CMS_minindex], dataT[4][CMS_minindex]


#-----  CMS expected exclusion --> find the minimum of -2lnL ------#

# CMS_chi2min = min(dataT[20])
# CMS_minindex = dataT[20].index(CMS_chi2min) 
# 
# print "Best fit point: ", CMS_chi2min, CMS_minindex, dataT[3][CMS_minindex], dataT[4][CMS_minindex]#

#-----  Naive Combination --> the sum of all available likelihoods ------#
#
# llh_comb_obs = map(lambda l1,l2,l3 : l1+l2+l3, dataT[17],dataT[19],dataT[21])
# llh_comb_exp = map(lambda l1,l2,l3 : l1+l2+l3, dataT[18],dataT[20],dataT[22])
#
# comb_chi2min = min(llh_comb_obs)
# comb_minindex = llh_comb_obs.index(comb_chi2min)
# print "Best fit point: ", comb_chi2min, comb_minindex, dataT[3][comb_minindex], dataT[4][comb_minindex]
#
#------

#------ plot ATLAS observed exclusion likelihood  -----#
# 
# cplot = ax.scatter(dataT[3],dataT[4], c=map(lambda d: d-ATL_chi2min, dataT[21]),marker='s',vmin=0.,vmax=20., alpha=1,edgecolors='none')
# tcplot = ax.tricontour(dataT[3],dataT[4],map(lambda d: d-ATL_chi2min, dataT[21]),levels=[5.9],colors=['k'],linewidths=2)
# ax.plot(ATLobsT[0],ATLobsT[1], linewidth=2, color='r', linestyle="--",label=r'$95\%~\mathrm{C.L.}~\mathrm{excl.}~(\mathrm{official})$')
# plottitle = r'$h/H/A\to\tau^+\tau^-~\mathrm{observed~exclusion~likelihood~(ATLAS)}$'
# 
#------ plot ATLAS expected exclusion likelihood  -----#

# cplot = ax.scatter(dataT[3],dataT[4], c=map(lambda d: d-ATL_chi2min, dataT[22]),marker='s',vmin=0.,vmax=20., alpha=1,edgecolors='none')
# tcplot = ax.tricontour(dataT[3],dataT[4],map(lambda d: d-ATL_chi2min, dataT[22]),levels=[5.9],colors=['k'],linewidths=2)
# ax.plot(ATLexpT[0],ATLexpT[1], linewidth=2, color='r', linestyle="--", label=r'$95\%~\mathrm{C.L.}~\mathrm{excl.}~(\mathrm{official})$')
# plottitle = r'$h/H/A\to\tau^+\tau^-~\mathrm{expected~exclusion~likelihood~(ATLAS)}$'

#------ plot CMS observed exclusion likelihood  -----#

# cplot = ax.scatter(dataT[3],dataT[4], c=map(lambda d: d-CMS_chi2min, dataT[19]),marker='s',vmin=0.,vmax=20., alpha=1,edgecolors='none')
# tcplot = ax.tricontour(dataT[3],dataT[4],map(lambda d: d-CMS_chi2min, dataT[19]),levels=[5.99],colors=['k'],linewidths=2)
# ax.plot(CMSobsT[0],CMSobsT[1], linewidth=2, color='r', linestyle="--",label=r'$95\%~\mathrm{C.L.}~\mathrm{excl.}~(\mathrm{official})$')
# plottitle = r'$h/H/A\to\tau^+\tau^-~\mathrm{observed~exclusion~likelihood~(CMS)}$'

cplot = ax.scatter(dataT[3],dataT[4], c=map(lambda d: d-CMS_chi2min, dataT[17]),marker='s',vmin=0.,vmax=20., alpha=1,edgecolors='none')
tcplot = ax.tricontour(dataT[3],dataT[4],map(lambda d: d-CMS_chi2min, dataT[17]),levels=[5.99],colors=['k'],linewidths=2)
ax.plot(CMSobsT[0],CMSobsT[1], linewidth=2, color='r', linestyle="--",label=r'$95\%~\mathrm{C.L.}~\mathrm{excl.}~(\mathrm{official})$')
plottitle = r'$h/H/A\to\tau^+\tau^-~\mathrm{observed~exclusion~likelihood~(CMS)}$'


#------ plot CMS expected exclusion likelihood  -----#

# cplot = ax.scatter(dataT[3],dataT[4], c=map(lambda d: d-CMS_chi2min, dataT[20]),marker='s',vmin=0.,vmax=20., alpha=1,edgecolors='none')
# tcplot = ax.tricontour(dataT[3],dataT[4],map(lambda d: d-CMS_chi2min, dataT[20]),levels=[5.9],colors=['k'],linewidths=2)
# ax.plot(CMSexpT[0],CMSexpT[1], linewidth=2, color='r', linestyle="--", label=r'$95\%~\mathrm{C.L.}~\mathrm{excl.}~(\mathrm{official})$')
# plottitle = r'$h/H/A\to\tau^+\tau^-~\mathrm{expected~exclusion~likelihood~(CMS)}$'


#------ plot combined observed exclusion likelihood  -----#

# cplot = ax.scatter(dataT[3],dataT[4], c=map(lambda d: d-comb_chi2min, llh_comb_obs),marker='s',vmin=0.,vmax=20., alpha=1,edgecolors='none')
# tcplot = ax.tricontour(dataT[3],dataT[4], map(lambda d: d-comb_chi2min, llh_comb_obs),levels=[5.9],colors=['k'],linewidths=2)

cbar = plt.colorbar(cplot)
cbar.set_label(r'$-2\,\mathrm{ln}(L)$')

tclabel = r'$95\%~\mathrm{C.L.}~\mathrm{excl.}~(\mathrm{from}~-2\mathrm{ln}(L))$'
tcplot.collections[0].set_label(tclabel)

ax.set_xlim([min(dataT[3]),max(dataT[3])])
ax.set_ylim([min(dataT[4]),max(dataT[4])])

legend = ax.legend(loc='lower right', fontsize= 12)

plt.grid()

plt.xlabel(r'$M_A~[\mathrm{GeV}]$')
plt.ylabel(r'$\tan\beta$')
plt.title(plottitle, fontsize=14)

plt.tight_layout()

# plt.savefig('mhmodp_llh_ATLAS_observed.pdf')
# plt.savefig('mhmodp_llh_ATLAS_expected.pdf')
plt.savefig('mhmodp_llh_CMS_observed.pdf')
# plt.savefig('mhmodp_llh_CMS_expected.pdf')
# plt.savefig('mhmodp_llh_combined_observed.pdf')

