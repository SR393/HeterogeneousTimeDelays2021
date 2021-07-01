import numpy as np
import scipy.io as io
import matplotlib as mpl
import matplotlib.pyplot as plt

folder = r'D:\Honours Results\Coherence Exemplars\cut'
locs = io.loadmat(folder+r'\locs.mat')['locs']
numnodes = int(locs.shape[1])
x = locs[0, :]
y = locs[1, :]
z = locs[2, :]

x_range = [np.amin(x), np.amax(x)]
x_mid = 0.5*(x_range[1] + x_range[0])
y_range = [np.amin(y), np.amax(y)]
y_mid = 0.5*(y_range[1] + y_range[0])
z_range = [np.amin(z), np.amax(z)]
z_mid = 0.5*(z_range[1] + z_range[0])

nodes = np.arange(0, numnodes)

under = np.where(z < z_mid)[0]
over = np.delete(nodes, under)

classes = np.array(['Synchronised', 'Coherent', 'Partially Coherent', 'Incoherent'])

fig, axes = plt.subplots(2, 2, figsize = (14, 10), gridspec_kw = {'wspace':0, 'hspace':0.2})

cmap = 'twilight'
norm = mpl.colors.Normalize(vmin=-np.pi,vmax=np.pi)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

for i, ax in enumerate(axes.flat):

    phases = io.loadmat(folder+'\\'+classes[i]+'.mat')['phases']

    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.set_title(classes[i], size = 30)
    frame = ax.scatter(-x, y, c = phases[0, :], cmap = cmap, s = 200)

fig.subplots_adjust(right=0.7)
cbar_ax = fig.add_axes([0.75, 0.15, 0.04, 0.7])
cbar = fig.colorbar(sm, cax = cbar_ax)
cbar.set_label(label='Phase (rad)', size=30)
cbar.ax.tick_params(labelsize=25)
plt.savefig(folder+r'\Poster_Coherence_Class_Figure.png', dpi = 400)


