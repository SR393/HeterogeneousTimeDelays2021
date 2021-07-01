
import numpy as np
import scipy.io as io
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as mcm
from mpl_toolkits.mplot3d import Axes3D

folder = r'D:\Honours Results\Flow_Exemplar_flows\flow_modes_exemplars'
locs = io.loadmat(folder+r'\locs.mat')['locs']
x = locs[0, :]
y = locs[1, :]
z = locs[2, :]

flownames = np.array(['Travelling', 'Source', 'Sink', 'Rotating', 'Diverging', 'Complex'])

fig, axes = plt.subplots(2, 3, figsize = (14, 7), subplot_kw={'projection':'3d'}, gridspec_kw = {'wspace':0, 'hspace':0.1})

cmap = mcm.cividis
norm = mpl.colors.Normalize(vmin=0.01,vmax=0.17)
sm = plt.cm.ScalarMappable(cmap='cividis', norm=norm)
sm.set_array([])

for i, ax in enumerate(axes.flat):

    data = io.loadmat(folder+'\\'+flownames[i]+'.mat')['flow_modes']
    data = data[0, 0]['V'][0, 0]
    x_flows = data['vx']
    x_flows_dominant_mode = x_flows[:, 0]
    y_flows = data['vy']
    y_flows_dominant_mode = y_flows[:, 0]
    z_flows = data['vz']
    z_flows_dominant_mode = z_flows[:, 0]
    flow_mags = np.sqrt(x_flows_dominant_mode**2+y_flows_dominant_mode**2+z_flows_dominant_mode**2)

    q = ax.quiver(x, z, y, x_flows_dominant_mode, z_flows_dominant_mode, y_flows_dominant_mode, length = 10, arrow_length_ratio = 0.5, 
                  linewidths = 0.8, normalize = True)

    q.set_color(cmap(norm(flow_mags)))
    ax.set_title(flownames[i], size = 25, y = 0.85)
    ax.grid(False)
    ax.set_axis_off()
    ax.view_init(elev=0, azim=270)
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(np.min(y), np.max(y))
    ax.set_zlim(np.min(z), np.max(z))

fig.subplots_adjust(right=0.7)
cbar_ax = fig.add_axes([0.75, 0.15, 0.04, 0.7])
cbar = fig.colorbar(sm, cax = cbar_ax)
cbar.set_label(label='Flow Speed (m/s)', size=25, labelpad=10)
cbar.ax.tick_params(labelsize=20)
plt.savefig(folder+r'\Poster_Flows_Figure.png', dpi = 400)


