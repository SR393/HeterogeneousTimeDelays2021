# Collection of functions to make simulation scripts shorter and more 
# easily changeable

from tvb.simulator.lab import *
import pycohcorrgram as pyccg
import os as os
import numpy as np
import scipy.signal as scpsig
import scipy.stats as spst
import scipy.io as io
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.cm as mcm
from matplotlib import colors
from tqdm import tqdm


def trunc_n(x, n):
    return np.trunc(x*10**n)/(10**n)

def set_osc_pars(strength, QV_Max = 1.0, VT = 0.0, delta_V = 0.65):

# Changes default oscillator parameters to ones used in Roberts et al., 2019, allows
# coupling, max V firing rate, V 'threshold' and regional threshold variance to be set

    coupled_C = np.array([strength])
    QV_Max = np.array([QV_Max])
    VT = np.array([VT])
    delta_V = np.array([delta_V])

    oscillator = models.LarterBreakspear(gCa = np.array([1.0]), gNa = np.array([6.7]),
                                        gK= np.array([2.0]), gL = np.array([0.5]),
                                        VCa = np.array([1.0]), VNa = np.array([0.53]),
                                        VK = np.array([-0.7]), VL = np.array([-0.5]),
                                        TCa = np.array([-0.01]), TNa = np.array([0.3]),
                                        TK = np.array([0.0]), d_Ca = np.array([0.15]),
                                        d_Na = np.array([0.15]), d_K = np.array([0.3]),
                                        VT = VT, ZT = np.array([0.0]),
                                        d_V = delta_V, d_Z = np.array([0.65]),
                                        QV_max = QV_Max, QZ_max = np.array([1.0]),
                                        b = np.array([0.1]), phi = np.array([0.7]),
                                        rNMDA = np.array([0.25]), aee = np.array([0.36]),
                                        aie = np.array([2.0]), ane = np.array([1.0]),
                                        aei = np.array([2.0]), ani = np.array([0.4]),
                                        Iext = np.array([0.3]), tau_K = np.array([1.0]),
                                        C = coupled_C, t_scale = np.array([1.0]))

    return oscillator


def set_init_coup(weights, a = 0.5, b = 1.0, midpoint = 0.0, sigma = 0.65):

# Sets initial coupling parameters and connection weights

    a = np.array([a])
    b = np.array([b])
    midpoint = np.array([midpoint])
    sigma = np.array([sigma])

    global coupling

    coupling = coupling.HyperbolicTangent(a = a, midpoint = midpoint,
                                          b = b, sigma = sigma)

    in_strengths = weights.sum(axis=1) 
    weights /= in_strengths

    return coupling, weights


def structural_sequence(ls, Is, tracts):

#sets up fibre tract matrices according to structural sequence defined by Weibull dist.
    
    #Set up for fibre length sequences
    numnodes = tracts.shape[0]
    whole_tracts = tracts.flatten()
    whole_tracts = whole_tracts[whole_tracts != 0]
    f_sequence = np.argsort(np.argsort(whole_tracts))
    x = np.arange(0, len(whole_tracts))/len(whole_tracts)

    fits = np.zeros((len(ls), len(Is), numnodes, numnodes))

    #Create fibre length matrices
    for i in range(len(ls)):
        maxlen = ls[i]*6 #max length when max(p_n) = inf
        for k in range(len(Is)):
            f_i = ls[i]*np.abs(np.log(1-x))**Is[k]
            f_i[f_i > maxlen] = np.full((len(f_i[f_i > maxlen])), maxlen)
            f_i = f_i[f_sequence]
            f_i = np.reshape(f_i, (numnodes, (numnodes-1)))
            for j in range(numnodes):
                fits[i, k, j, 0:j] = f_i[j, 0:j]
                fits[i, k, j, j+1:] = f_i[j, j:]

    return fits

def surrogate_mats_Weibull(tracts, x):

#Creates surrogate fibre length matrices by randomising x% of connections while preserving
#their statistical properties, assumed to follow a Weibull distribution
    if x != 0:

        #1. Pick random set of fibre tracts to be replaced, size = x%*len(flattened tracts matrix)
        numnodes = tracts.shape[0]
        flat_tracts = tracts.flatten()
        num_rand_lengths = int(np.rint((x/100)*len(flat_tracts)))
        rng = np.random.default_rng()
        zeroes = (numnodes+1)*np.arange(0, numnodes)
        inds = rng.choice(np.delete(np.arange(0, len(flat_tracts)), zeroes), num_rand_lengths, replace = False)
        inds = inds.astype(np.int32)
        inds = np.sort(inds)
        
        #2. Determine Weibull dist. parameters of fibre tract sample
        fibre_lengths_sample = flat_tracts[inds]
        shape, loc, scale = spst.weibull_min.fit(fibre_lengths_sample, floc = 0)
        weibull_sample = spst.weibull_min.rvs(shape, loc = 0, scale = scale, size = num_rand_lengths)

        #3. Randomly reassign tract lengths of sample tracts according to Weibull dist. with same parameters
        flat_tracts[inds] = weibull_sample
        new_tracts = np.reshape(flat_tracts, (int(np.sqrt(len(flat_tracts))), int(np.sqrt(len(flat_tracts))))) #assumes square matrix

    else:

        new_tracts = tracts

    return new_tracts

def get_euc_dists(positions):

    numnodes = int(positions.shape[0])
    Euclidean_dist = np.empty((numnodes, numnodes))

    for i in range(numnodes):
        node_i = positions[i]
        for j in range(numnodes):
            if i != j:
                node_j = positions[j]
                diff = node_i - node_j
                dist = np.abs(np.sqrt(np.sum(diff**2)))
                Euclidean_dist[i, j] = dist
            else:
                Euclidean_dist[i, j] = 0
    
    return Euclidean_dist

def coarse_grain(time, data, avg_segment, cut_length):

#Sliding window time average of input data with no overlap. Cuts transient of length cut_length, then divides
#data into segments of size avg_segment to obtain coarse-grained dataset. Returns coarse grained data and time
#array. avg_segment and cut_length are input in time units.

    numnodes = data.shape[1]
    dt = time[1]-time[0]
    
    avg_segment = int(avg_segment/dt)
    cut_length = int(cut_length/dt)
    time = time[cut_length:]
    data = data[cut_length:, :]
    time_indices = avg_segment*np.arange(0, len(time)/avg_segment)
    time_indices = time_indices.astype(np.int32)
    time = time[time_indices]
    TAVG = np.mean(data.reshape(int(len(time)), avg_segment, numnodes), axis = 1)

    return TAVG, time

def make_net_animations(time, data, positions, fps, save_folder):

# Makes (cubic) 'net' animations of input data given region positions to view activity on all sides. Inputs are:
# time: array of timepoints at which data is evaluated
# data: input data
# positions: region locations
# ffmpeg_path: path to ffmpeg directory; must have ffmpeg installed for function to work 
# save_folder: directory to save animation 

    mpl.rcParams['animation.ffmpeg_path'] = r'C:\Users\Sebastian\miniconda3\envs\TVBEnv\ffmpeg-20200807-fab00b0-win64-static\bin\ffmpeg.exe'
    writer = animation.FFMpegWriter(fps=fps)

    positions = positions.transpose()
    numnodes = int(positions.shape[1])
    x = positions[0]
    y = positions[1]
    z = positions[2]

    #determine which 'side' each region is part of by splitting brain into octants
    x_range = [np.amin(x), np.amax(x)]
    x_mid = 0.5*(x_range[1] + x_range[0])
    y_range = [np.amin(y), np.amax(y)]
    y_mid = 0.5*(y_range[1] + y_range[0])
    z_range = [np.amin(z), np.amax(z)]
    z_mid = 0.5*(z_range[1] + z_range[0])

    nodes = np.arange(0, numnodes)

    back = np.where(y < y_mid)[0]
    front = np.delete(nodes, back)
    under = np.where(z < z_mid)[0]
    over = np.delete(nodes, under)
    left = np.where(x < x_mid)[0]
    right = np.delete(nodes, left)

    print('Making Animation')

    #Set up figure axes for 'brain net' 
    fig = plt.figure(figsize = (12, 8))
    grid = plt.GridSpec(3, 4, wspace = 0, hspace = 0)
    ax2_frontview = plt.subplot(grid[0, 1])
    ax2_frontview.set_xticks([])
    ax2_frontview.set_yticks([])
    ax2_frontview.set_title('Front')
    ax2_leftview = plt.subplot(grid[1, 0])
    ax2_leftview.set_xticks([])
    ax2_leftview.set_yticks([])
    ax2_leftview.set_title('Left')
    ax2_overview = plt.subplot(grid[1, 1])
    ax2_overview.set_xticks([])
    ax2_overview.set_yticks([])
    ax2_rightview = plt.subplot(grid[1, 2])
    ax2_rightview.set_xticks([])
    ax2_rightview.set_yticks([])
    ax2_rightview.set_xlabel('Right')
    ax2_backview = plt.subplot(grid[2, 1])
    ax2_backview.set_xticks([])
    ax2_backview.set_yticks([])
    ax2_backview.set_xlabel('Back')
    ax2_underview = plt.subplot(grid[1, 3])
    ax2_underview.set_xticks([])
    ax2_underview.set_yticks([])
    ax2_underview.set_title('Bottom')
    timestamp = plt.subplot(grid[0, 2])
    timest = timestamp.text(0.5, 0.5, '0 ms', transform = timestamp.transAxes, fontsize=15)
    timestamp.set_xticks([])
    timestamp.set_yticks([])

    cmap = mcm.cividis
    norm = colors.Normalize(vmin = np.min(data), vmax = np.max(data))

    #Initial frames
    region_size = 150 #size of circles representing brain region
    frontview = ax2_frontview.scatter(-x[front], z[front], c = data[0, front], cmap = 'cividis', s = region_size)
    backview = ax2_backview.scatter(x[back], z[back], c = data[0, back], cmap = 'cividis', s = region_size)
    leftview = ax2_leftview.scatter(-y[left], z[left], c = data[0, left], cmap = 'cividis', s = region_size)
    rightview = ax2_rightview.scatter(y[right], z[right], c = data[0, right], cmap = 'cividis', s = region_size)
    overview = ax2_overview.scatter(x[over], y[over], c = data[0, over], cmap = 'cividis', s = region_size)
    underview = ax2_underview.scatter(x[under], y[under], c = data[0, under], cmap = 'cividis', s = region_size)

    #make and save animations frame-by-frame
    with writer.saving(fig, save_folder, 100):
        for n in tqdm(range(int(len(time)))):
            data_n = data[n, :]
            frontview.set_color(cmap(norm(data_n[front])))
            backview.set_color(cmap(norm(data_n[back])))
            overview.set_color(cmap(norm(data_n[over])))
            underview.set_color(cmap(norm(data_n[under])))
            leftview.set_color(cmap(norm(data_n[left])))
            rightview.set_color(cmap(norm(data_n[right])))
            timest.set_text(str(time[n])+' ms')
            writer.grab_frame()

    plt.close('all')

def threshold_weights(weights, x):

#Removes all elements of matrix 'weights' smaller than the top x%

    weights_flat = weights.flatten()
    resort_indices = np.argsort(np.argsort(weights_flat))
    weights_flat = np.sort(weights_flat)
    no_weights_below_thresh = int(np.rint(((100-x)/100)*len(weights_flat)))
    weights_flat[0:no_weights_below_thresh] = np.zeros(no_weights_below_thresh)
    weights_flat = weights_flat[resort_indices]
    numnodes = weights.shape[0]
    weights = np.reshape(weights_flat, (numnodes, numnodes))

    return weights

def substitute_artificial_hist(length, oscillator, wm, integrator, what_to_watch):

#Substitutes constant history of length 'length' for each variable at each node,
#values chosen randomly from intervals specific to each variable

    sim = simulator.Simulator(model = oscillator, connectivity = wm,
                            coupling = coupling, integrator = integrator, 
                            monitors = what_to_watch, simulation_length = length,
                            conduction_speed = 100.00).configure()
    
    numnodes = wm.weights.shape[0]

    v_fixed = np.random.default_rng().uniform(-0.6, 0.2, numnodes)

    for i in range(numnodes):

        sim.history.buffer[:, 0, i, 0] = v_fixed[i]*np.ones(sim.history.buffer.shape[0])

    sim.current_state[0, :] = v_fixed.reshape(numnodes, 1)
    sim.current_state[1, :] = np.random.default_rng().uniform(0.0, 0.6, numnodes).reshape(numnodes, 1)
    sim.current_state[2, :] = np.random.default_rng().uniform(-0.03, 0.13, numnodes).reshape(numnodes, 1)

    return sim

def run_history_sim(length, folder, oscillator, wm, integrator, what_to_watch):

#Runs simulation of length 'length' to determine the history of each variable at each node, decoupling
#all nodes so that the history for each is identical to all others.

    sim = simulator.Simulator(model = oscillator, connectivity = wm,
                            coupling = coupling, integrator = integrator, 
                            monitors = what_to_watch, simulation_length = length).configure()

    sim.model.C = np.array([0.0])

    print("Creating Simulation History")

    #Set initial conditions
    #Replace history 
    v_fixed = 0.0
    sim.history.buffer = v_fixed * np.ones(sim.history.buffer.shape)
    w_fixed = 0.0
    z_fixed = 0.0

    sim.current_state[0, ...] = v_fixed
    sim.current_state[1, ...] = w_fixed
    sim.current_state[2, ...] = z_fixed
    (time_trans, data_trans), = sim.run()
    data_trans = data_trans[:, 0, :, 0]

    #plot history to check
    plt.plot(time_trans, data_trans, 'k')
    plt.title('History of All Oscillators, Variable V')
    plt.xlabel('Time (ms)')
    plt.savefig(folder+r'\History.png')
    plt.close('all')

    return sim

def disconnect_lr_hemispheres(weights):

#Sets all interhemispheric weights to 0

    numnodes = weights.shape[0]
    weights[:int(numnodes/2), int(numnodes/2):numnodes] = np.zeros((int(numnodes/2), int(numnodes/2)))
    weights[int(numnodes/2):numnodes, :int(numnodes/2)] = np.zeros((int(numnodes/2), int(numnodes/2)))

    return weights

def disconnect_ap_hemispheres(weights, positions):

#Sets all weights between anterior/posterior hemispheres to 0

    numnodes = weights.shape[0]
    positions = positions.transpose()
    y = positions[1]

    y_range = [np.amin(y), np.amax(y)]
    y_mid = 0.5*(y_range[1] + y_range[0])

    nodes = np.arange(0, numnodes)

    back = np.where(y < y_mid)[0]
    front = np.delete(nodes, back)

    numnodes = weights.shape[0]
    weights[np.ix_(front, back)] = np.zeros((int(len(front)), int(len(back))))
    weights[np.ix_(back, front)] = np.zeros((int(len(back)), int(len(front))))

    return weights

def disconnect_sections(weights, positions, x_pos = np.array([None]), y_pos = np.array([None]), z_pos = np.array([None])):

# Sets all weights between sections defined by x_pos, y_pos and z_pos to 0. Inputs are weights matrix, node
# positions, and arrays with fraction of distance from midpoint of the brain - eg, x = np.array([1/3, -1/3])
# splits network into even thirds along the x direction, x = np.array([0]) splits network in half, etc.

    numnodes = weights.shape[0]
    nodes = np.arange(0, numnodes)
    positions = positions.transpose()

    plane_list = [x_pos, y_pos, z_pos]

    #loop through coordinates
    for i in range(3):

        co = plane_list[i]

        if not all(co == None):

            #calculate position in network coordinates, shifted so midpoint = origin
            coordinate = positions[i]
            print(co)
            range_ = [np.amin(coordinate), np.amax(coordinate)]
            midpoint = 0.5*(range_[0] + range_[1])
            pos = co*(range_[0] - range_[1])/2 + midpoint

            for j in range(len(pos)):

                one_side = np.where(positions[i] < co[j])[0]
                other_side = np.delete(nodes, one_side)

                numnodes = weights.shape[0]
                weights[np.ix_(one_side, other_side)] = np.zeros((int(len(one_side)), int(len(other_side))))
                weights[np.ix_(other_side, one_side)] = np.zeros((int(len(other_side)), int(len(one_side))))

    return weights

def subnetwork_average(data, subset):

#Take nodewwise average activity of subset of nodes

    data = data[:, subset]
    avg_data = np.mean(data, axis = 1)

    return avg_data

def plot(x_axis, data, save_loc, title = None, x_label = None, y_label = None, show = False):

#Bunches up a few plt functions to save space at the cost of utility

    plt.plot(x_axis, data)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)

    if show:
        plt.show()
    
    else:
        plt.savefig(save_loc)
    
    plt.close('all')

def surrogate_mats_shuffle(tracts, positions, no_bins = 100):

#Creates surrogate fibre length matrices which preserve statistical properties and background trends
#by detrending, shuffling, then adding background trends back. 'Trends' are general relationships
#between mean and euclidean distance and standard deviation and euclidean distance

    #Bins defined by intervals of euclidean distances
    Euc_dists = get_euc_dists(positions)
    max_euc_dist = np.max(Euc_dists)

    euc_dists_flat = Euc_dists.flatten()
    euc_dists_flat = euc_dists_flat[euc_dists_flat != 0]
    tracts_flat = tracts.flatten()
    tracts_flat = tracts_flat[tracts_flat != 0]

    numnodes = tracts.shape[0]

    means = np.empty(no_bins)
    stds = np.empty(no_bins)
    subtracted_tracts = np.empty(numnodes*(numnodes-1))

    for i in range(no_bins):

        ni1 = np.nonzero(euc_dists_flat <= (max_euc_dist/no_bins)*(i+1))
        ni2 = np.nonzero(euc_dists_flat > (max_euc_dist/no_bins)*i)
        ni = np.intersect1d(ni1, ni2)
        means[i] = np.mean(tracts_flat[ni])
        stds[i] = np.std(tracts_flat[ni])

        #standardise sample from each bin individually
        subtracted_tracts[ni] = tracts_flat[ni] - means[i]
        subtracted_tracts[ni] = subtracted_tracts[ni]/stds[i]

    #create indices to shuffle tract lengths
    rng = np.random.default_rng()
    inds = rng.choice(numnodes*(numnodes-1), numnodes*(numnodes-1), replace = False)
    inds = inds.astype(np.int32)
    subtracted_tracts = subtracted_tracts[inds]

    #Retrend
    for i in range(no_bins):

        ni1 = np.nonzero(euc_dists_flat <= (max_euc_dist/no_bins)*(i+1))
        ni2 = np.nonzero(euc_dists_flat > (max_euc_dist/no_bins)*i)
        ni = np.intersect1d(ni1, ni2)

        subtracted_tracts[ni] = subtracted_tracts[ni]*stds[i]
        subtracted_tracts[ni] = subtracted_tracts[ni] + means[i]

    #new tracts array with 0 'diagonal' (after reshaping)
    new_tracts = np.zeros(numnodes**2)

    for i in range(numnodes):

        new_tracts[numnodes*i+(i+1):numnodes*(i+1)+(i+1)] = subtracted_tracts[i*numnodes:(i+1)*numnodes]
    
    new_tracts = np.reshape(new_tracts, (numnodes, numnodes))

    return new_tracts

def create_IHCC(data, numnodes, opts, saveloc):

    left = np.arange(0, int(numnodes/2))
    right = np.arange(int(numnodes/2), numnodes)
    [IHCC, ord1, ord2, ct, cl] = pyccg.pycohcorrgram(data, left, right, opts, False) #obtain IHCC
    io.savemat(saveloc+'.mat', {'IHCC':IHCC})

def IHCC_plot(data_loc, numnodes, opts, simlength, dt, saveloc):

    #Inputs to pycohcorrgram
    maxlag = opts['maxlag']
    window_steps = opts["winlen"]/dt
    window_gap = opts["winlen"]*opts["overlapfrac"]/dt
    simsteps = simlength/dt
    n_windows = int(np.ceil((simsteps - window_steps)/(window_steps - window_gap))) + 1

    #For IHCC tick markers
    xticklocs = np.floor(np.linspace(0, n_windows, 6))
    for i in range(6):
        xticklocs[i] = int(xticklocs[i])

    yticklocs = np.floor(np.linspace(0, 2*opts["maxlag"]/dt + 1, 5))
    for i in range(5):
        yticklocs[i] = int(yticklocs[i])
    
    IHCC = io.loadmat(data_loc+'.mat')['IHCC']

    plt.imshow(IHCC, aspect = 'auto', cmap = 'coolwarm')
    plt.xticks(xticklocs, trunc_n(np.linspace(0, simlength, len(xticklocs)), 2))
    plt.yticks(yticklocs, trunc_n(np.linspace(-maxlag, maxlag, len(yticklocs)), 2))
    plt.xlabel('Time (ms)')
    plt.ylabel('Lag (ms)')
    col = plt.colorbar()
    col.set_label('Cross Corr.', rotation = 270)
    plt.savefig(saveloc+'.png')
    plt.close('all')

def identify_transitions_IHCC(data_loc, filename_base, thrsh, dt):

# Identifies state transitions using IHCC as points at which IHCC drops below thrsh

    IHCC = io.loadmat(data_loc)['IHCC']
    transition_times = []

    for i in range(IHCC.shape[1]):

        timepoint = IHCC[:, i]

        if all(np.abs(timepoint) < thrsh):

            transition_times.append(i*dt)

    return np.array(transition_times)

def pyphase_nodes(y):

    phases = np.zeros(y.shape)

    for i in range(y.shape[1]):

        cent = y[:, i] - np.mean(y[:, i])
        phases[:, i] = np.unwrap(np.angle(scpsig.hilbert(cent)))

    return phases

def order_parameters(data, thrshs, positions, saveloc):

    #define euclidean distances between nodes
    #to determine which fall inside radius of inclusion 
    Euc_dist = get_euc_dists(positions)

    simsteps, numnodes = data.shape[0], data.shape[1]

    #store data to be saved
    Rlocals = np.empty(len(thrshs))
    cnt = 0
    for r in tqdm(thrshs):

        groups = []

        #group defined by nodes within r of node i
        for i in range(numnodes):
            groups.append(np.nonzero(Euc_dist[i] < r)[0])
        
        orderparameters = np.empty((simsteps, numnodes))

        #calculate R for every group
        for k in range(len(groups)):

            p = np.exp(1j*pyphase_nodes(data[:, groups[k]]))
            coh = np.absolute(np.mean(p, axis = 1))
            orderparameters[:, k] = coh
        
        #calc. R_local
        R_local = np.mean(np.mean(orderparameters, axis = 1))

        Rlocals[cnt] = R_local

        cnt += 1

    #save to plot later
    io.savemat(saveloc, {'rlocals':Rlocals, 'radii':thrshs})

def rlocal_plots(data, thrshs, saveloc):

    data = data[0]
    data[0] = 1 #First point otherwise NAN
    #line = np.linspace(ordpardata[0], ordpardata[-1], 13) #Line connecting curve ends

    #make plot for each R_local(r) curve
    plt.plot(thrshs, data)
    #plt.plot(thrshs, line, 'k--')
    plt.ylabel(r'$R_{local}$', size = 15)
    plt.xlabel('Radius of Inclusion (mm)', size = 15)
    plt.ylim(0.1, 1.1)
    plt.title('Local Order Parameter vs Radius of Local Group', size = 17)
    plt.savefig(saveloc)

def plot_in_lengths(tracts, positions, saveloc, show = False):

    in_lengths = np.mean(tracts, axis = 0)

    fig1 = plt.figure(figsize = (22, 22))
    length_ax = fig1.add_subplot(111, projection = '3d')
    length_ax.grid(False)
    length_ax.set_xticks([])
    length_ax.set_yticks([])
    length_ax.set_zticks([])
    scat = length_ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2], c = in_lengths, cmap = 'coolwarm', s = 1000)
    scat.set_array(in_lengths)
    col = plt.colorbar(scat)
    col.set_label('Mean In-Length (mm)', size = 20, rotation = 270, labelpad = 30)
    if show:
        plt.show()
    plt.savefig(saveloc)

def inst_freq(data, time):

    phases = pyphase_nodes(data)
    phases_shifted = np.roll(phases, -1, axis = 0)
    dphase = phases_shifted - phases
    dphase = dphase[:len(phases)-1]
    dt = np.roll(time, -1) - time
    dt = dt[:len(time)-1]
    dt = dt.reshape((dt.shape[0], 1))
    inst_freqs = dphase/dt

    return inst_freqs, time[:len(time)-1]