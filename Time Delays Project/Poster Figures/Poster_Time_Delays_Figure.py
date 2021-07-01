import numpy as np
import scipy.io as io
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.path as pltPath

image_1_dir = r'D:\Images\White_matter_tracts.jpg'
image_1 = plt.imread(image_1_dir)
image_1 = image_1[:, :360]

#image_2_dir = r''
#image_2 = plt.imread(image_2_dir)

fig = plt.figure(figsize = (10, 10))

ax = fig.add_subplot(111)
ax.set_axis_off()
ax.imshow(image_1)
ax.set_aspect('equal')
pos1 = ax.get_position()
pos2 = [pos1.x0 - 0.05, pos1.y0,  pos1.width, pos1.height] 
ax.set_position(pos2)

#i
region_i = patches.Circle((110, 164), 15, ec = 'black', fill = False, lw = 2.5)
ax.add_patch(region_i)
ax.text(106.5, 169.2, 'i', size = 30, family = 'serif')

#j
region_j = patches.Circle((261, 115), 15, ec = 'black', fill = False, lw = 2.5)
ax.add_patch(region_j)
ax.text(259, 119, 'j', size = 30, family = 'serif')

#tract path 1
vertices_1 = np.array([(101.5, 151.4), (92, 139), (94, 133), (94, 126), (97, 120), (100, 115), (121.5, 103.5), (142, 95.5),  
                       (168.5, 87.5), (185.5, 84.5), (207.5, 90), (230, 101), (245, 106.5), (241.3, 105.3)])

codes_1 = np.array([pltPath.Path.MOVETO, pltPath.Path.CURVE3, pltPath.Path.CURVE3, pltPath.Path.CURVE3, pltPath.Path.CURVE3, 
                    pltPath.Path.CURVE3, pltPath.Path.CURVE3, pltPath.Path.CURVE3, pltPath.Path.CURVE3, pltPath.Path.CURVE3, 
                    pltPath.Path.CURVE3, pltPath.Path.CURVE3, pltPath.Path.CURVE3, pltPath.Path.LINETO])

path_1 = patches.PathPatch(path = pltPath.Path(vertices_1, codes = codes_1), fill = False, ls = '--', lw = 2)
ax.add_patch(path_1)

arrowhead_1_position = (241.3, 105.3)
ax.add_patch(patches.RegularPolygon(arrowhead_1_position, 3, 5, 175*np.pi/180, color = 'black'))

ax.text(177.1, 75, r'$L_{ij}$', size = 35)

#time delay 1
#arrow_1 = patches.FancyArrow(117.9, 160.2, 250.9 - 117.9, 119.4 - 160.2, color = 'black', length_includes_head = True, 
                             #head_width = 4, head_length = 4, lw = 0.8)
#ax1.add_patch(arrow_1)

#k
region_k = patches.Circle((127, 290.5), 15, ec = 'black', fill = False, lw = 2.5)
ax.add_patch(region_k)
ax.text(121, 296.5, 'k', size = 30, family = 'serif')

#l
region_l = patches.Circle((129.7, 17), 15, ec = 'black', fill = False, lw = 2.5)
ax.add_patch(region_l)
ax.text(127.3, 25, 'l', size = 30, family = 'serif')

#tract path 2
vertices_2 = np.array([(129.6, 275.5), (130.1, 276.1), (134.6, 264.3), (137.5, 251.2), (144.8, 240.2), (148.5, 231.6), (150.9, 212.8), (153.8, 196.5),  
                       (159.5, 181.4), (163.2, 166.3), (159.5, 146.3), (155, 123.1), (148.5, 95.7), (143.2, 75.3), (134.6, 43.9), (133.4, 36.8)])

codes_2 = np.array([pltPath.Path.MOVETO, pltPath.Path.CURVE3, pltPath.Path.CURVE3, pltPath.Path.CURVE3, pltPath.Path.CURVE3, 
                    pltPath.Path.CURVE3, pltPath.Path.CURVE3, pltPath.Path.CURVE3, pltPath.Path.CURVE3, pltPath.Path.CURVE3, 
                    pltPath.Path.CURVE3, pltPath.Path.CURVE3, pltPath.Path.CURVE3, pltPath.Path.CURVE3, pltPath.Path.CURVE3, 
                    pltPath.Path.LINETO])

path_2 = patches.PathPatch(path = pltPath.Path(vertices_2, codes = codes_2), fill = False, ls = '--', lw = 2)
ax.add_patch(path_2)

arrowhead_2_position = (133.4, 36.8)
ax.add_patch(patches.RegularPolygon(arrowhead_2_position, 3, 5, 51*np.pi/180, color = 'black'))

ax.text(167, 169, r'$L_{kl}$', size = 35)

#time delay 2
#arrow_2 = patches.FancyArrow(123.6, 290, 126.9 - 123.6, 21.4 - 290, color = 'black', length_includes_head = True, 
                             #head_width = 4, head_length = 4, lw = 0.8)
#ax1.add_patch(arrow_2)

caption_ax = fig.add_axes([0.1, 0.05, 0.6, 0.1])
caption_ax.set_axis_off()
caption_ax.text(0.02, 0.2, r'$L_{ij} < L_{kl}$'+r' '+r'&'+r' '+r'$s_{xy} = s$'+' '+r'$\forall x, y \Longrightarrow \tau_{ij} < \tau_{kl}$', size = 30)

#plt.show()
plt.savefig(r'D:\Honours Results\Poster_Figure\Poster_Time_Delays_Figure.png', dpi = 300)