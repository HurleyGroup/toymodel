from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
import numpy as np
import cloudpickle
import argparse
import os, sys

# read in arguments and important variables
parser = argparse.ArgumentParser()
parser.add_argument("-p", "--postpath", help="specify the post path (output of stress rise)", type=str)
parser.add_argument("-a", "--ca", help="specify the configuration angle", type=float)
parser.add_argument("-i", "--imgpath", help="specify the img path to save plots", type=str)
args = parser.parse_args()

POST_PATH = args.postpath
IMG_PATH = args.imgpath
THETA = args.ca

PI = 3.141592654

MU, R = float(os.environ['MU']), float(os.environ['R'])

print('Loading solution.')
with open(POST_PATH+'sol.pkl', mode='rb') as file:
    soln = cloudpickle.load(file)

theta = THETA*PI/180.

# obtain particle positions and velocities
times = soln.t # times are same for all

soln_x1, soln_y1, soln_y2  = soln.y[0], soln.y[1], soln.y[2]
soln_vx1, soln_vy1, soln_vy2  = soln.y[3], soln.y[4], soln.y[5]

# we will only look at what is going on when velocity is <= 0 for y velocities and >= 0 for x velocities
mask_x1 = soln_vx1>=0
mask_y1 = soln_vy1<=0
mask_y2 = soln_vy2<=0

# now apply the mask accordingly
times_x1 = times[mask_x1]
times_y1 = times[mask_y1]
times_y2 = times[mask_y2]

soln_x1, soln_y1, soln_y2 = soln_x1[mask_x1], soln_y1[mask_y1], soln_y2[mask_y2]
soln_vx1, soln_vy1, soln_vy2 = soln_vx1[mask_x1], soln_vy1[mask_y1], soln_vy2[mask_y2]

print('x1 last time = ', times_x1[-1])
print('y1 last time = ', times_y1[-1])
print('y2 last time = ', times_y2[-1])

earliest_time = min(times_x1[-1], times_y1[-1], times_y2[-1])
global_times = times[times<=earliest_time]
NUM_STEPS = len(global_times)

# potential energy over time
print('Loading Functions.')
with open(POST_PATH+'UTotal.pkl', mode='rb') as file:
    UTotal_lambda = cloudpickle.load(file)

with open(POST_PATH+'Us.pkl', mode='rb') as file:
    Us_lambda = cloudpickle.load(file)

with open(POST_PATH+'Ut.pkl', mode='rb') as file:
    Ut_lambda = cloudpickle.load(file)

with open(POST_PATH+'cf1.pkl', mode='rb') as file:
    cf1_lambda = cloudpickle.load(file)

with open(POST_PATH+'cf2.pkl', mode='rb') as file:
    cf2_lambda = cloudpickle.load(file)

with open(POST_PATH+'cfs.pkl', mode='rb') as file:
    cfs_lambda = cloudpickle.load(file)


def compute_pe(funcT, funcS, funcTan, x1, y1, y2):
    Ut = np.zeros(NUM_STEPS)
    Us = np.zeros(NUM_STEPS)
    Utan = np.zeros(NUM_STEPS)
    for t in range(NUM_STEPS):
        Ut[t] = funcT(x1[t], y1[t], y2[t])
        Us[t] = funcS(x1[t], y1[t], y2[t])
        Utan[t] = funcTan(x1[t], y1[t], y2[t])
    return Ut, Us, Utan

def compute_stress(fc1, fc2, fc3, x1, y1, y2):
    c1 = [-R*np.sin(theta), -R*np.cos(theta)]
    c2 = [R, 0]
    c3 = [-R*np.sin(theta), R*np.cos(theta)]

    return (np.outer(fc1, c1) + np.outer(fc2, c2) + np.outer(fc3, c3))/(PI*(R**2))


def compute_force(funcCF1, funcCF2, funcCS, x1, y1, y2):
    cf1 = np.zeros([NUM_STEPS, 2])
    cf2 = np.zeros([NUM_STEPS, 2])
    cfs = np.zeros([NUM_STEPS, 2])
    stress = np.zeros([NUM_STEPS, 2, 2])
    for t in range(NUM_STEPS):
        cf1[t] = np.asarray(funcCF1(x1[t], y1[t], y2[t])).flatten()
        cf2[t] = np.asarray(funcCF2(x1[t], y1[t], y2[t])).flatten()
        cfs[t] = np.asarray(funcCS(x1[t], y1[t], y2[t])).flatten()

        # contact force 1 multiplied by negative 1 since compression
        stress[t] = compute_stress(-1*cf1[t], cfs[t], cf2[t], x1[t], y1[t], y2[t])

    return stress


UTotal, Us, Ut = compute_pe(UTotal_lambda, Us_lambda, Ut_lambda,
                            soln_x1, soln_y1, soln_y2)

stresses = compute_force(cf1_lambda, cf2_lambda, cfs_lambda,
                         soln_x1, soln_y1, soln_y2)

final_state = np.asarray([soln_x1[:NUM_STEPS], soln_y1[:NUM_STEPS], soln_y2[:NUM_STEPS],
                          soln_vx1[:NUM_STEPS], soln_vy1[:NUM_STEPS], soln_vy2[:NUM_STEPS]])



# find overlap, if any
def find_overlap(Ulat, Ustrong, first=1):
    if np.all(Ustrong[first:]>Ulat[first:]):
        return 'strong'
    elif np.all(Ulat[first:]>Ustrong[first:]):
        return 'confused'
    else:
        # work backwards to find intersection
        
        # check after rates change
        mask = np.diff(Ulat)>np.diff(Ustrong)
        start_arr = int(sum(~mask)/2) # be conservative
        Ulat_search = Ulat[start_arr:]
        Ustr_search = Ustrong[start_arr:]

        plt.figure()
        plt.plot(Ulat)
        plt.plot(Ustrong)
        plt.savefig('/code/test.png')

        # find where it changes


    print(Ulat.shape)
    print(Ustrong.shape)
    print(global_times.shape)



find_overlap(Us, (UTotal-Us)/2)
sys.exit()




# Save final state
np.save(POST_PATH+'final_state.npy', final_state, allow_pickle=True)






## Plotting
fig, (ax1, ax2) = plt.subplots(1,2, figsize=(10,5), constrained_layout=True)
ax1.plot(global_times, (UTotal-Us)/2, label='strong')
ax1.plot(global_times, Us, label='confused')
ax1.set_xlabel('Time')
ax1.set_ylabel('Potential Energy (per contact)')
ax1.legend()

ax3 = ax2.twinx()
ax2.plot(global_times, -1*stresses[:,1,1]*2*R, label=r"$\sigma_{yy}$", color='green')
ax3.plot(global_times, np.diff(-1*stresses[:,1,1]*2*R), label=r"$\dot{\sigma_{yy}}$", color='red')
ax2.set_xlabel('Time')
ax2.set_ylabel('Stress Units')
ax3.set_ylabel('Stress Change with Time')
lines, labels = ax2.get_legend_handles_labels()
lines2, labels2 = ax3.get_legend_handles_labels()
ax2.legend(lines + lines2, labels + labels2, loc=4)

plt.savefig(IMG_PATH+'PEnStress.png')







#
