import matplotlib.pyplot as plt
import pylab as P
import numpy as np

P.rc('text', usetex=True)
P.rc('text.latex', preamble='\usepackage{amsmath}\usepackage{color}')
font = {'size' : 14}
P.rc('font', **font)
P.rc('grid', linewidth=1,color='#666666')

data1 = np.loadtxt("2Higgses_pdf1.dat")
data1T = zip(*data1)

data2 = np.loadtxt("2Higgses_pdf2.dat")
data2T = zip(*data2)

data3 = np.loadtxt("2Higgses_pdf3.dat")
data3T = zip(*data3)


plt.close('all')

fig, axes = plt.subplots(nrows = 2, ncols = 3,sharex='col',sharey='row')

# fig.clf()
# for a in axes:
# 	for aa in a:
# 		aa.axis('off')
# plt.tick_params(axis='both', left='off', top='off', right='off', bottom='off', labelleft='off', labeltop='off', labelright='off', labelbottom='off')
# plt.show()


plts = [[0,0,data1T,8,"Box pdf (1)"],\
        [0,1,data2T,8,"Gaussian pdf (2)"],\
        [0,2,data3T,8,"Box+Gaussian pdf (3)"],\
        [1,0,data1T,11,""],\
        [1,1,data2T,11,""],\
        [1,2,data3T,11,""]]

# plts = [[0,0,data1T,7,"Box pdf (1)"],\
#         [0,1,data2T,7,"Gaussian pdf (2)"],\
#         [0,2,data3T,7,"Box+Gaussian pdf (3)"],\
#         [1,0,data1T,11,""],\
#         [1,1,data2T,11,""],\
#         [1,2,data3T,11,""]]


for p in plts:
# 	print p
	
# 	axes[p[0],p[1]].cla()
	
# 	axes[p[0],p[1]].scatter(p[2][0],p[2][1])
	
	chi2min = min(p[2][p[3]])
	minindex = p[2][p[3]].index(chi2min)

	print chi2min,  p[2][0][minindex],p[2][1][minindex]
	
# 	if p[0] == 1:
# 		make these tick labels invisible
# 		axes[p[0],p[1]].cla()
# 		plt.setp(axes[p[0],p[1]].get_xticklabels(), visible=False)
	
	# 1/2sigma filled regions:
	axes[p[0],p[1]].tricontourf(p[2][0],p[2][1],map(lambda d: d-chi2min, p[2][p[3]]),levels=np.arange(0,10,step=0.5),cmap=plt.cm.YlOrRd)#,colors=['green','yellow','none']
	# 1/2sigma contours:
	axes[p[0],p[1]].tricontour(p[2][0],p[2][1],map(lambda d: d-chi2min, p[2][p[3]]),levels=[2.3,5.99],colors=['k','k'],linewidths=2)
	axes[p[0],p[1]].plot(p[2][0][minindex],p[2][1][minindex],'*',markersize=10,color='w') 

	axes[p[0],p[1]].set_xlim([min(p[2][0]),max(p[2][0])])
	axes[p[0],p[1]].set_ylim([min(p[2][1]),max(p[2][1])])

	axes[p[0],p[1]].set_xticks(np.arange(min(p[2][0]),max(p[2][0])+1,step=1.))
	axes[p[0],p[1]].set_yticks(np.arange(min(p[2][1]),max(p[2][1])+1,step=1.))	
# 	axes[p[0],p[1]].relim()
	axes[p[0],p[1]].grid()

 	if p[0] == 1:
 		axes[p[0],p[1]].set_xlabel(r'$M_1~[\mathrm{GeV}]$')
 	if p[1] == 0:
 		axes[p[0],p[1]].set_ylabel(r'$M_2~[\mathrm{GeV}]$')
 	
 	if p[0] == 0:
 		axes[p[0],p[1]].set_title(p[4],fontsize=14)	
# 		make these tick labels invisible
# 		plt.setp(axes[p[0],p[1]].get_xticklabels(), visible=False)
# 	axes[p[0],p[1]].label_outer()
		
plt.subplots_adjust(hspace=0.10, wspace=0.20)		
# 
# fig.tight_layout()



plt.savefig('2Higgses.pdf')

