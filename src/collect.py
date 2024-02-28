from matplotlib import pyplot as plt
import numpy as np
import cloudpickle
import argparse
import os

plt.rcParams.update({'font.size': 15})

# read in arguments and important variables
parser = argparse.ArgumentParser()
parser.add_argument("-mp", "--modelpath", help="specify the model path", type=str)
parser.add_argument("-i", "--imgpath", help="specify the img path to save plots", type=str)
args = parser.parse_args()

MODEL_PATH = args.modelpath
IMG_PATH = args.imgpath

PI = 3.141592654
MU, R = float(os.environ['MU']), float(os.environ['R'])

LOAD = float(os.environ['LOAD'])
LOAD_TIME = float(os.environ['LOAD_TIME'])
P = float(os.environ['PRESSURE'])

KS = float(os.environ['KS_START'])


print('Loading Functions.')
with open(MODEL_PATH+'UTotal.pkl', mode='rb') as file:
    UTotal_lambda = cloudpickle.load(file)

with open(MODEL_PATH+'Us.pkl', mode='rb') as file:
    Us_lambda = cloudpickle.load(file)

with open(MODEL_PATH+'Ut.pkl', mode='rb') as file:
    Ut_lambda = cloudpickle.load(file)

with open(MODEL_PATH+'cf1.pkl', mode='rb') as file:
    cf1_lambda = cloudpickle.load(file)

with open(MODEL_PATH+'cf2.pkl', mode='rb') as file:
    cf2_lambda = cloudpickle.load(file)

with open(MODEL_PATH+'cfs.pkl', mode='rb') as file:
    cfs_lambda = cloudpickle.load(file)




POST_PATH_1000 = '/code/post/maintain/4560_490_6000_1000/'  
POST_PATH_N1000 = '/code/post/maintain/4560_490_6000_-1000/'  

# load data
with open(POST_PATH_1000+'solMAINTAIN.pkl', mode='rb') as file:
    soln1000 = cloudpickle.load(file)

with open(POST_PATH_N1000+'solMAINTAIN.pkl', mode='rb') as file:
    solnN1000 = cloudpickle.load(file)


with open(POST_PATH_1000+'init_dof.pkl', 'rb') as ff:
    init_dof_1000 = cloudpickle.load(ff)

with open(POST_PATH_N1000+'init_dof.pkl', 'rb') as ff:
    init_dof_N1000 = cloudpickle.load(ff)

with open(POST_PATH_1000+'init_offset.pkl', 'rb') as ff:
    init_offset_1000 = cloudpickle.load(ff)

with open(POST_PATH_N1000+'init_offset.pkl', 'rb') as ff:
    init_offset_N1000 = cloudpickle.load(ff)

with open(POST_PATH_1000+'NODDYH.pkl', 'rb') as ff:
    H_1000 = cloudpickle.load(ff)

with open(POST_PATH_N1000+'NODDYH.pkl', 'rb') as ff:
    H_N1000 = cloudpickle.load(ff)



init_dofs = [init_dof_1000, init_dof_N1000]
init_offsets = [init_offset_1000, init_offset_N1000]
h_arrays = [H_1000, H_N1000]


# make hydrostatic plot
plt.figure()
labelList = ['1000', '-1000']
labelListR = [r'$\dot{k}^s$ = 1000', r'$\dot{k}^s$ = -1000']
for ss, soln in enumerate([soln1000, solnN1000]):
    times, soln_x1, soln_y1, soln_y2, soln_vx1, soln_vy1, soln_vy2  = soln

    init_dof = init_dofs[ss] # somehow find a way to update these
    init_offset = init_offsets[ss]
    h_array = h_arrays[ss]

    if float(labelList[ss]) > 0:
        THETA_LIMIT = 0.
    else:
        THETA_LIMIT = 60.

    plt.plot(h_array[0], h_array[1]/h_array[1][0], label=labelListR[ss])
    plt.xlabel('Time (s)')
    plt.ylabel(r'Hydrostatic Stress Ratio, $\frac{\sigma_h}{\sigma^0_h}$')

plt.legend()
plt.tight_layout()
plt.savefig('/code/hydrostatic.png')

# 
