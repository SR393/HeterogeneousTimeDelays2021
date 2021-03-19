from tvb.simulator.lab import *
import function_collection as fc
import os as os
import numpy as np
import scipy.io as io

folder = r'D:\Results\Paper\TimeDelays\Nodes_512\WeightThresholded'

if not os.path.exists(folder):
    os.makedirs(folder)

network_513 = 'connectivity_513.zip'
network_512 = 'geombinsurr_partial_000_Gs.zip'

#set up connectome
wm = connectivity.Connectivity.from_file(network_512)
positions = wm.centres
tracts = wm.tract_lengths
numnodes = tracts.shape[0]

oscillator = fc.set_osc_pars(strength = 0)

#simulation parameters
dt = 0.5
integrator = integrators.Dopri5(dt = dt)
simlength = 40000.00
simsteps = simlength/dt
mon_raw = monitors.Raw()
what_to_watch = (mon_raw, )

#global dynamics parameters
coupling, wm.weights = fc.set_init_coup(wm.weights)
wm.weights = fc.threshold_weights(wm.weights, 20)
np.fill_diagonal(wm.weights, 0)
wm.speed = np.array([5])

coups = np.array([0.6, 0.5])

for c in coups:

    savestring = r'Coupling='+str(c)

    folder1 = folder+'\\'+savestring
    if not os.path.exists(folder1):
        os.makedirs(folder1)

    #store data here
    matlab_files = folder+r'\Matlab Files'
    if not os.path.exists(matlab_files):
            os.makedirs(matlab_files)

    sim = fc.substitute_artificial_hist(10, oscillator, wm, integrator, what_to_watch)

    print('Started Simulation')
    sim.simulation_length = simlength
    sim.model.C = np.array([c])
    (time, data), = sim.run()
    data = data[:, 0, :, 0]

    io.savemat(matlab_files+r'\\'+savestring+'.mat', {'data':data, 'locs':positions, 'time_vec':time})

    cg_time_step = 2
    transient = 5000
    data, time = fc.coarse_grain(time, data, dt, cg_time_step, transient, numnodes)

    io.savemat(matlab_files+r'\\'+savestring+'TAVG.mat', {'data':data, 'locs':positions, 'time_vec':time})

    fc.make_net_animations(time, data, positions, 100, r'C:\Users\Sebastian\miniconda3\envs\TVBEnv\ffmpeg-20200807-fab00b0-win64-static\bin\ffmpeg.exe', 
                         folder1 + r'\Anim_tavg'+savestring+'NET.mp4')




