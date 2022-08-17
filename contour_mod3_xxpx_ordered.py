import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from pylab import *
from matplotlib.cm import ScalarMappable
from matplotlib import rcParams
X00=np.loadtxt("../seqa/dihedral_ree_local_4th_dih.dat")
X0=np.loadtxt("../seqb/dihedral_ree_local_rightP_at_3rdpos.dat")
X1= np.loadtxt("../seqc_modify/dihedral_ree_local.dat")
X2=np.loadtxt("../seqd/dihedral_ree_local_rightP_at_3rdpos.dat")
X3=np.loadtxt("../seqe/dihedral_ree_local_rightP_at_3rdpos.dat")
X4=np.loadtxt("../seqx/dihedral_ree_local_rightP_at_3rdpos.dat")
X5=np.loadtxt("../seqy/dihedral_ree_local_rightP_at_3rdpos.dat")
X6=np.loadtxt("../seqz/dihedral_ree_local_rightP_at_3rdpos.dat")
###########################################################
mappables=[]
mappables1=[]
mappables2=[]
H00,xedges,yedges = np.histogram2d(X00[:,0], X00[:,1], bins=25, range=[[-3.2, 3.2],[0.5, 1.1]])#, density='true')
mappables.append(log(0.00001+H00/len(X00[:,1])))
mappables1.append(log(0.00001+H00/len(X00[:,1])))
H0,xedges,yedges = np.histogram2d(X0[:,0], X0[:,1], bins=25, range=[[-3.2, 3.2],[0.5, 1.1]])#, density='true')
mappables.append(log(0.00001+H0/len(X0[:,1])))
mappables1.append(log(0.00001+H0/len(X0[:,1])))
H1, xedges,yedges = np.histogram2d(X1[:,0], X1[:,1], bins=25, range=[[-3.2, 3.2],[0.5, 1.1]])#, density ='true')
mappables.append(log(0.00001+H1/len(X1[:,1])))
mappables2.append(log(0.00001+H1/len(X1[:,1])))
H2, xedges,yedges = np.histogram2d(X2[:,0], X2[:,1], bins=25, range=[[-3.2, 3.2],[0.5, 1.1]])#, density ='true')
mappables.append(log(0.00001+H2/len(X2[:,1])))
mappables2.append(log(0.00001+H2/len(X2[:,1])))
H3, xedges,yedges = np.histogram2d(X3[:,0], X3[:,1], bins=25, range=[[-3.2, 3.2],[0.5, 1.1]])#, density ='true')
mappables.append(log(0.00001+H3/len(X3[:,1])))
mappables2.append(log(0.00001+H3/len(X3[:,1])))
H4, xedges,yedges = np.histogram2d(X4[:,0], X4[:,1], bins=25, range=[[-3.2, 3.2],[0.5, 1.1]])#, density ='true')
mappables.append(log(0.00001+H4/len(X4[:,1])))
mappables2.append(log(0.00001+H4/len(X4[:,1])))
H5, xedges,yedges = np.histogram2d(X5[:,0], X5[:,1], bins=25, range=[[-3.2, 3.2],[0.5, 1.1]])#, density ='true')
mappables.append(log(0.00001+H5/len(X5[:,1])))
mappables2.append(log(0.00001+H5/len(X5[:,1])))
H6, xedges,yedges = np.histogram2d(X6[:,0], X6[:,1], bins=25, range=[[-3.2, 3.2],[0.5, 1.1]])#, density ='true')
mappables.append(log(0.00001+H6/len(X6[:,1])))
mappables2.append(log(0.00001+H6/len(X6[:,1])))
vmin = (np.min(mappables))
vmax = (np.max(mappables))
levels = 20

level_boundaries = np.linspace(vmin, vmax, levels + 1)
print(vmin, vmax)

fig, axes = plt.subplots(2,3, figsize=(10,5))
rc('axes', linewidth=3)
rcParams['ytick.labelsize'] = 16
for ax,H in zip(axes.ravel(),mappables2):
    #ax.set_xlabel("$r_{X-X-P-X}(nm)$", fontsize=16, fontweight='bold')
    #ax.set_ylabel(r'$\phi_{X-X-P-X}$', fontsize=16, fontweight='bold')

    ax.set_xticks((0.6, 0.8,1.0),fontsize=20)
    ax.set_yticks((-2, 0, 2), fontsize=20)
    ax.set_xticklabels((0.6,0.8,1.0), fontsize=20)
    ax.set_yticklabels((-2, 0, 2), fontsize=20)
    ax.set_box_aspect(1)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2)
    im = ax.contourf(H,vmin = vmin, vmax =vmax, extent=[0.5,1.1,-3.2,3.2], levels=levels, cmap='CMRmap_r')
     
fig.tight_layout()
#cb=fig.colorbar(im, ax=axes.ravel())
#cb.boundaries = np.array(0,2,4,6,8,12,14))
#cb.set_ticks((0,4,8,12))
fig.colorbar(
    ScalarMappable(norm=im.norm, cmap=im.cmap),ax=axes.ravel(), ticks=[-10, -8, -6, -4, -2],
    
    boundaries=level_boundaries,
    values=(level_boundaries[:-1] + level_boundaries[1:]) / 2,
)
plt.savefig("xxpx_six_plots.png",bbox_inches = "tight", dpi=300 )
#plt.show()
plt.close()

#################################################
fig, axes = plt.subplots(1,2, figsize=(10,5))
rc('axes', linewidth=3)
rcParams['ytick.labelsize'] = 16
for ax,H in zip(axes.ravel(),mappables1):
    #ax.set_xlabel("$r_{X-X-P-X}(nm)$", fontsize=16, fontweight='bold')
    #ax.set_ylabel(r'$\phi_{X-X-P-X}$', fontsize=16, fontweight='bold')

    ax.set_xticks((0.6, 0.8,1.0),fontsize=20)
    ax.set_yticks((-2, 0, 2), fontsize=20)
    ax.set_xticklabels((0.6,0.8,1.0), fontsize=25)
    ax.set_yticklabels((-2, 0, 2), fontsize=25)
    ax.set_box_aspect(1)
    im = ax.contourf(H,vmin = vmin, vmax =vmax, extent=[0.5,1.1,-3.2,3.2], levels=levels, cmap='CMRmap_r')

fig.tight_layout()
#cb=fig.colorbar(im, ax=axes.ravel())
#cb.boundaries = np.array(0,2,4,6,8,12,14))
#cb.set_ticks((0,4,8,12))
fig.colorbar(
    ScalarMappable(norm=im.norm, cmap=im.cmap),ax=axes.ravel(), ticks=[-10, -8, -6, -4, -2],

    boundaries=level_boundaries,
    values=(level_boundaries[:-1] + level_boundaries[1:]) / 2,
)
plt.savefig("xxpx_ref_plots.png",bbox_inches = "tight", dpi=300 )
#plt.show()
plt.close()
