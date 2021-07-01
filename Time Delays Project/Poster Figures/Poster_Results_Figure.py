from matplotlib.colors import ListedColormap
import numpy as np
import scipy.io as io
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

folder = r'D:\Honours Results\Flow_Exemplar_flows\CoupSpeed_flows'
locs = io.loadmat(folder+r'\locs.mat')['locs']
x = locs[0, :]
y = locs[1, :]
z = locs[2, :]

couplings = np.array([8, 7, 6, 5])
speeds = np.array([40, 60, 80, 100, 120, 140])

flownames = []

for C in couplings:
    for sp in speeds:

        flownames.append(r'Sp'+str(sp)+r'Coup0'+str(C)+r'_flows')

flownames = np.array(flownames)

fig, axes = plt.subplots(4, 6, figsize = (14, 8), subplot_kw = {'projection':'3d'}, gridspec_kw = {'wspace':0, 'hspace':0})
plt.gcf().subplots_adjust(bottom = 0.17, top = 0.93)

classes = np.array([1, 1, 1, 2, 2, 2, 
                    0, 1, 2, 2, 2, 2, 
                    0, 0, 1, 2, 2, 1,
                    0, 0, 0, 0, 0, 0])

cmap = ListedColormap(["red", "purple", "blue"])
norm = mpl.colors.Normalize(vmin=0,vmax=2)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
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

    q = ax.quiver(x, z, y, x_flows_dominant_mode, z_flows_dominant_mode, y_flows_dominant_mode, length = 10, arrow_length_ratio = 0.5, 
                  linewidths = 0.3, normalize = True)

    class_colour = classes[i]*np.ones(512)
    q.set_color(cmap(norm(class_colour)))
    ax.grid(False)
    ax.set_axis_off()
    ax.view_init(elev=0, azim=270)
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(np.min(y), np.max(y))
    ax.set_zlim(np.min(z), np.max(z))

ax = fig.add_axes([0.07, 0.15, 0.63, 0.000001])
ax.set_xlabel('Conduction Speed (m/s)', size = 30, labelpad = 15)
ax.axes.get_yaxis().set_visible(False)
ax.set_xticks(np.linspace(1/6, 1, 6))
ax.set_xticklabels(['40', '60', '80', '100', '120', '140'], size = 25)

ax = fig.add_axes([0.1, 0.05, 0.000001, 0.78])
ax.set_ylabel('Connection Strength (a.u.)', size = 30, labelpad = 15, y = 0.6)
ax.axes.get_xaxis().set_visible(False)
ax.set_yticks(np.linspace(0.25, 1, 4))
ax.set_yticklabels(['0.5', '0.6', '0.7', '0.8'], size = 25)

fig.subplots_adjust(right=0.75)
cbar_ax = fig.add_axes([0.775, 0.35, 0.04, 0.3])
cbar = fig.colorbar(sm, cax = cbar_ax)
cbar.set_ticks([1/3, 1, 5/3])
cbar.set_ticklabels(['Complex', 'Diverging', 'Rotating'])
cbar.ax.tick_params(labelsize=25)
#plt.show()
plt.savefig(folder + r'\Poster_Results_Figure.png', dpi = 300)


