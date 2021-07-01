import numpy as np
import scipy.io as io
import matplotlib.pyplot as plt
import matplotlib.patches as patches

folder = r'D:\Honours Results\LarterBreakspear\LastSims'

K = 0.8
speed = 100

folder1 = folder+r'\Coupling = '+str(K)
matlab_files = folder1+r'\Matlab Files'

datafile = io.loadmat(matlab_files+r'\Sp'+str(speed)+r'Coup'+str(K)+r'.mat')
data = datafile['data']
positions = datafile['locs']
positions = positions.transpose()
numnodes = int(positions.shape[1])
x = positions[0]
y = positions[1]
z = positions[2]

x_range = [np.amin(x), np.amax(x)]
x_mid = 0.5*(x_range[1] + x_range[0])
y_range = [np.amin(y), np.amax(y)]
y_mid = 0.5*(y_range[1] + y_range[0])
z_range = [np.amin(z), np.amax(z)]
z_mid = 0.5*(z_range[1] + z_range[0])

nodes = np.arange(0, numnodes)

under = np.where(z < z_mid)[0]
over = np.delete(nodes, under)

plotted_nodes = nodes

#Set up figure axes for 'brain net' 
region_size = 150
no_frames = 4

fig = plt.figure(figsize = (16, 5))
plt.suptitle('"Rotating" Pattern of Neural Activity', size = 35)
grid = plt.GridSpec(1, no_frames, wspace = 0.04, hspace = 0)
im = plt.imshow(data, cmap = 'cividis')

for i in range(no_frames):

    ax = plt.subplot(grid[0, i])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.set_xlabel('t = '+str(i*2)+' ms', size = 30, labelpad = 20)
    ax.scatter(-x[plotted_nodes], y[plotted_nodes], c = data[24032+32*i, plotted_nodes], cmap = 'cividis', s = region_size)
    ax.set_aspect('equal')
    arc = patches.Arc([0, -15], 90, 90, angle = 0, theta1 = 15, theta2 = 345, capstyle = 'round', linestyle = '-', lw = 2, color = 'black')
    ax.add_patch(arc)

    arrowhead_position = (45*np.cos(15*np.pi/180), -15 + 45*np.sin(15*np.pi/180))

    ax.add_patch(patches.RegularPolygon(arrowhead_position, 3, 10, 205*np.pi/180, color = 'black'))

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.01, 0.7])
cbar = fig.colorbar(im, cax = cbar_ax)
cbar.set_label(label = 'Activity Amp. (mV)', size = 30, labelpad = 20)
cbar.ax.tick_params(labelsize=25)
plt.savefig(folder+r'\Rotating_Pattern_Frames.png', dpi = 300)
