from matplotlib.animation import FuncAnimation, FFMpegWriter
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
from matplotlib import pyplot as plt
import numpy as np
import cloudpickle
import argparse
import os, sys

# read in arguments and important variables
parser = argparse.ArgumentParser()
parser.add_argument("-mp", "--modelpath", help="specify the model path", type=str)
parser.add_argument("-p", "--postpath", help="specify the post path (output of stress rise)", type=str)
parser.add_argument("-i", "--imgpath", help="specify the img path to save plots", type=str)
parser.add_argument('--maintain', action=argparse.BooleanOptionalAction, help="whether to constrain hydrostatic pressure as constant or not")
args = parser.parse_args()

MODEL_PATH = args.modelpath
POST_PATH = args.postpath
IMG_PATH = args.imgpath
MAINTAIN = args.maintain # If true, then we apply constraint

PI = 3.141592654

MU, R = float(os.environ['MU']), float(os.environ['R'])


print('Loading solution.')
if MAINTAIN:
    with open(POST_PATH+'solMAINTAIN.pkl', mode='rb') as file:
        soln = cloudpickle.load(file)
else:
    with open(POST_PATH+'solNOMAINTAIN.pkl', mode='rb') as file:
        soln = cloudpickle.load(file)

print('Loading Constants')
LOAD = float(os.environ['LOAD'])
LOAD_TIME = float(os.environ['LOAD_TIME'])
P = float(os.environ['PRESSURE'])

KS = float(os.environ['KS_START'])
KS_RATE_RISE = float(os.environ['KS_RATE_RISE'])
KS_TIME_RISE = float(os.environ['KS_TIME_RISE'])
KS_RATE_DROP = float(os.environ['KS_RATE_DROP'])

with open(POST_PATH+'init_dof.pkl', 'rb') as ff:
    init_dof = cloudpickle.load(ff)

with open(POST_PATH+'init_offset.pkl', 'rb') as ff:
    init_offset = cloudpickle.load(ff)

x10, y10, y20 = init_dof # somehow find a way to update these


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

def save_animation(positions):
    # Define the number of frames and time interval between frames
    num_frames = len(positions[0])  # Assuming all arrays have the same length
    # fps = 1/np.mean(np.diff(times))  # 1 / 60 (60 fps)

    # Create a figure and axis for the animation
    fig, ax = plt.subplots()

    # Create empty circles for each object
    circles = [plt.Circle((0, 0), R, fc='r'),
               plt.Circle((positions[0][0], positions[1][0]), R, fc='r'),
               plt.Circle((0, positions[2][0]), R, fc='r')]

    # Add circles to the axis
    for circle in circles:
        ax.add_patch(circle)

    # Set axis limits
    ax.set_xlim(-0.05, 0.15)
    ax.set_ylim(-0.05, 0.2)
    ax.set_aspect('equal')

    # Function to update the positions of the circles in each frame
    def update(frame):
        circles[1].center = [positions[0][frame], positions[1][frame]]
        circles[2].center = [0.0, positions[2][frame]]
        return circles

    # Create the animation
    ani = FuncAnimation(fig, update, frames=num_frames, blit=True)

    # Specify the output format and create the writer (60 fps)
    writer = FFMpegWriter(fps=30, metadata=dict(artist='Adyota Gupta'), bitrate=1800)

    # Save the animation as an MP4 file
    if MAINTAIN:
        ani.save('/code/anim_MAINTAIN.mp4', writer=writer)
    else:
        ani.save('/code/anim_NOMAINTAIN.mp4', writer=writer)

    # Close the figure to release resources
    plt.close(fig)

# changes in PE
def compute_pe(funcT, funcS, funcTan, times, x1, y1, y2):
    Ut = np.zeros(NUM_STEPS)
    Us = np.zeros(NUM_STEPS)
    Utan = np.zeros(NUM_STEPS)
    for t in range(NUM_STEPS):
        Ut[t] = funcT(times[t], x1[t], y1[t], y2[t], LOAD, LOAD_TIME, KS, KS_RATE_RISE, KS_TIME_RISE, KS_RATE_DROP, P, x10, y10, y20, init_offset)
        Us[t] = funcS(times[t], x1[t], y1[t], y2[t], LOAD, LOAD_TIME, KS, KS_RATE_RISE, KS_TIME_RISE, KS_RATE_DROP, P, x10, y10, y20, init_offset)
        Utan[t] = funcTan(times[t], x1[t], y1[t], y2[t], LOAD, LOAD_TIME, KS, KS_RATE_RISE, KS_TIME_RISE, KS_RATE_DROP, P, x10, y10, y20, init_offset)
    return Ut-Ut[0], Us-Us[0], Utan-Utan[0]

# compute forces
def compute_force(funcCF1, funcCF2, funcCS, times, x1, y1, y2):
    cf1 = np.zeros([NUM_STEPS, 2])
    cf2 = np.zeros([NUM_STEPS, 2])
    cfs = np.zeros([NUM_STEPS, 2])
    stress = np.zeros([NUM_STEPS, 2, 2])
    for t in range(NUM_STEPS):
        cf1[t] = np.asarray(funcCF1(times[t], x1[t], y1[t], y2[t], LOAD, LOAD_TIME, KS, KS_RATE_RISE, KS_TIME_RISE, KS_RATE_DROP, P, x10, y10, y20, init_offset)).flatten()
        cf2[t] = np.asarray(funcCF2(times[t], x1[t], y1[t], y2[t], LOAD, LOAD_TIME, KS, KS_RATE_RISE, KS_TIME_RISE, KS_RATE_DROP, P, x10, y10, y20, init_offset)).flatten()
        cfs[t] = np.asarray(funcCS(times[t], x1[t], y1[t], y2[t], LOAD, LOAD_TIME, KS, KS_RATE_RISE, KS_TIME_RISE, KS_RATE_DROP, P, x10, y10, y20, init_offset)).flatten()

    return cf1, cf2, cfs

# compute triplet configuration angle
def getConfigAngle(x1, y1, y2):
    returned = np.zeros(NUM_STEPS)

    for t in range(NUM_STEPS):
        if np.sqrt(x1[t]**2 + y1[t]**2) > 2*R:
            continue
        if np.sqrt(x1[t]**2 + (y2[t]-y1[t])**2) > 2*R:
            continue

        returned[t] =  (PI/2 - np.arctan2(y1[t], x1[t]))*180./PI

    return returned

# compute stresses [triplet frame]
def compute_stress(fc1, fc2, fc3, x1, y1, y2):
    returned = np.zeros([NUM_STEPS, 2, 2])

    for t in range(NUM_STEPS):
        theta = (PI/2 - np.arctan2(y1[t], x1[t]))*180./PI
        c1 = [-R*np.sin(theta), -R*np.cos(theta)]
        c2 = [R, 0]
        c3 = [-R*np.sin(theta), R*np.cos(theta)]

        returned[t] = (np.outer(fc1[t], c1) + np.outer(fc2[t], c2) + np.outer(fc3[t], c3))/(PI*(R**2))

    return returned

def rotate_stress(stresses, ALPHA):
    ORIENT= ALPHA*PI/180.
    Q = np.array([[np.cos(ORIENT), np.sin(ORIENT)],[-np.sin(ORIENT), np.cos(ORIENT)]])
    oriented_stresses = np.zeros(stresses.shape)

    for t in range(len(stresses)):
        oriented_stresses[t] = np.einsum('ij,im,jn', stresses[t], Q, Q)

    return oriented_stresses

# obtain particle positions and velocities
times, soln_x1, soln_y1, soln_y2, soln_vx1, soln_vy1, soln_vy2  = soln
NUM_STEPS = len(times)

# create animation of particles buckling/fluttering
save_animation([soln_x1, soln_y1, soln_y2])

UTotal, Us, Ut = compute_pe(UTotal_lambda, Us_lambda, Ut_lambda,
                            times, soln_x1, soln_y1, soln_y2)

cf1, cf2, cfs =  compute_force(cf1_lambda, cf2_lambda, cfs_lambda,
                            times, soln_x1, soln_y1, soln_y2)

stresses = compute_stress(cf1, cfs, cf2, soln_x1, soln_y1, soln_y2)

oriented_stresses = rotate_stress(stresses, 9)


# plot stresses over time
plt.figure()
plt.plot(times, 0.5*(oriented_stresses[:,0,0]+oriented_stresses[:,1,1]))
if MAINTAIN:
    plt.savefig(IMG_PATH+'testH_MAINTAIN.png')
else:
    plt.savefig(IMG_PATH+'testH_NOMAINTAIN.png')

plt.figure()
plt.plot(times, oriented_stresses[:,0,0])
plt.plot(times, oriented_stresses[:,1,1])
plt.plot(times, oriented_stresses[:,0,1]+oriented_stresses[:,1,0])
if MAINTAIN:
    plt.savefig(IMG_PATH+'testStress_MAINTAIN.png')
else:
    plt.savefig(IMG_PATH+'testStress_NOMAINTAIN.png')


# plot configuration angle over time
theta_t = getConfigAngle(soln_x1, soln_y1, soln_y2)

plt.figure()
plt.plot(times, theta_t)
if MAINTAIN:
    plt.savefig(IMG_PATH+'testTheta_MAINTAIN.png')
else:
    plt.savefig(IMG_PATH+'testTheta_NOMAINTAIN.png')

# plot potential energies over time
Ustr = (UTotal-Us)/2
fig, (ax1, ax2) = plt.subplots(1,2, figsize=(10,5), constrained_layout=True)
ax1.plot(times, Ustr, label='strong', color='blue')
ax1.plot(times, Us, label='confused', color='orange')
ax1.set_xlabel('Time')
ax1.set_ylabel('Potential Energy Change (per contact)')
ax1.legend()

ax2.plot(times, -1*oriented_stresses[:,1,1]+oriented_stresses[0,1,1], label=r"$\sigma_{yy}$", color='green')
ax2.plot(times, -1*oriented_stresses[:,0,0]+oriented_stresses[0,0,0], label=r"$\sigma_{xx}$", color='red')
ax2.plot(times, oriented_stresses[:,0,1]-oriented_stresses[0,0,1], label=r"$\sigma_{xy}$", color='blue')
ax2.legend()

ax2.set_xlabel('Time')
ax2.set_ylabel('Stress Units')

lines, labels = ax2.get_legend_handles_labels()

if MAINTAIN:
    plt.savefig(IMG_PATH+'PEnStress_MAINTAIN.png')
else:
    plt.savefig(IMG_PATH+'PEnStress_NOMAINTAIN.png')





#
