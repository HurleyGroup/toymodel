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

with open(POST_PATH+'hD.pkl', mode='rb') as file:
   hD_lambda = cloudpickle.load(file)

def hD_lambdaN(a,b,c):
    g = hD_lambda(a,b,c)
    g1, g2 = g[0], g[1]
    g3, g4 = g[2], g[3]
    g5, g6 = g[4], g[5]
    return np.array([g1, g2, g3, g4, g5, g6]).reshape((3,2)).astype(np.float64)


theta = THETA*PI/180.

# obtain particle positions and velocities
times = soln.t # times are same for all

soln_x1, soln_y1, soln_y2  = soln.y[0], soln.y[1], soln.y[2]
soln_u1, soln_u2  = soln.y[3], soln.y[4]


soln_vx1, soln_vy1, soln_vy2  = np.zeros(len(times)),  np.zeros(len(times)),  np.zeros(len(times))

for t in range(len(times)):
    soln_vx1[t], soln_vy1[t], soln_vy2[t]  = hD_lambdaN(soln_x1[t], soln_y1[t], soln_y2[t]).dot(np.array([soln_u1[t], soln_u2[t]]).reshape(2,1))


global_times = times
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


# create animation of particles buckling/fluttering

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
    ani.save('/code/anim.mp4', writer=writer)

    # Close the figure to release resources
    plt.close(fig)

save_animation([soln_x1, soln_y1, soln_y2])


def compute_pe(funcT, funcS, funcTan, x1, y1, y2):
    Ut = np.zeros(NUM_STEPS)
    Us = np.zeros(NUM_STEPS)
    Utan = np.zeros(NUM_STEPS)
    for t in range(NUM_STEPS):
        Ut[t] = funcT(x1[t], y1[t], y2[t])
        Us[t] = funcS(x1[t], y1[t], y2[t])
        Utan[t] = funcTan(x1[t], y1[t], y2[t])
    return Ut-Ut[0], Us-Us[0], Utan-Utan[0]

def compute_stress(fc1, fc2, fc3, x1, y1, y2):
    c1 = [-R*np.sin(theta), -R*np.cos(theta)]
    c2 = [R, 0]
    c3 = [-R*np.sin(theta), R*np.cos(theta)]

    return (np.outer(fc1, c1) + np.outer(fc2, c2) + np.outer(fc3, c3))/(PI*(R**2))


def compute_force(funcCF1, funcCF2, funcCS, x1, y1, y2, vx1, vy1, vy2):
    cf1 = np.zeros([NUM_STEPS, 2])
    cf2 = np.zeros([NUM_STEPS, 2])
    cfs = np.zeros([NUM_STEPS, 2])
    stress = np.zeros([NUM_STEPS, 2, 2])
    for t in range(NUM_STEPS):
        cf1[t] = np.asarray(funcCF1(x1[t], y1[t], y2[t], vx1[t], vy1[t], vy2[t])).flatten()
        cf2[t] = np.asarray(funcCF2(x1[t], y1[t], y2[t], vx1[t], vy1[t], vy2[t])).flatten()
        cfs[t] = np.asarray(funcCS(x1[t], y1[t], y2[t], vx1[t], vy1[t], vy2[t])).flatten()


        # contact force 1 multiplied by negative 1 since compression
        stress[t] = compute_stress(cf1[t], cfs[t], cf2[t], x1[t], y1[t], y2[t])

    plt.figure()
    plt.plot(cfs[:])
    plt.savefig('/code/testCFS.png')

    plt.figure()
    plt.plot(cf2[:])
    plt.savefig('/code/testCF2.png')

    plt.figure()
    plt.plot(cf1[:])
    plt.savefig('/code/testCF1.png')

    return stress


UTotal, Us, Ut = compute_pe(UTotal_lambda, Us_lambda, Ut_lambda,
                            soln_x1, soln_y1, soln_y2)

stresses = compute_force(cf1_lambda, cf2_lambda, cfs_lambda,
                         soln_x1, soln_y1, soln_y2,
                         soln_vx1, soln_vy1, soln_vy2 )

final_state = np.asarray([soln_x1[NUM_STEPS-1], soln_y1[NUM_STEPS-1], soln_y2[NUM_STEPS-1],
                          soln_vx1[NUM_STEPS-1], soln_vy1[NUM_STEPS-1], soln_vy2[NUM_STEPS-1]])

plt.figure()
plt.plot(0.5*(stresses[:,0,0]+stresses[:,1,1]))
plt.savefig('/code/testH.png')


# find overlap, if any
def intersect(line1, line2, init_guess):
    intersection_func = lambda x: line1(x) - line2(x)
    intersection_point = fsolve(intersection_func, init_guess)

    return intersection_point



def find_overlap(Ulat, Ustr):
    if Ulat[-1] >= Ustr[-1]:
        sign = -1 # lateral larger than strong (at the end!)
    else:
        sign = 1 # strong larger than lateral (at the end!)

    if np.all(Ustr>Ulat) or np.all(Ulat>=Ustr): # no crossing
        return [global_times[-1], sign*np.absolute(stresses[-1,1,1]), False] # last stress (confused larger, negative; strong larger, positive
    else: # there's a crossing
        # work backwards to find intersection
        index = NUM_STEPS-1

        if sign == -1:
            while Ulat[index] >= Ustr[index]:
                index -= 1
        else:
            while Ustr[index] >= Ulat[index]:
                index -= 1

        # with crossing, interpolate to find stress!
        Ustr_interp_func = interp1d(global_times[index:(index+2)],Ustr[index:(index+2)])
        Ulat_interp_func = interp1d(global_times[index:(index+2)],Ulat[index:(index+2)])
        sigmayy_interp_func = interp1d(global_times[index:(index+2)],stresses[index:(index+2),1,1])

        itime = intersect(Ulat_interp_func, Ustr_interp_func, 0.5*np.diff(global_times[index:(index+2)])+global_times[index])

        return [itime, sign*np.absolute(sigmayy_interp_func(itime)), True]


def getConfigAngle(x1, y1, y2):
    PI = 3.141592654
    if np.sqrt(x1**2 + y1**2) > 2*R:
        return 0
    if np.sqrt(x1**2 + (y2-y1)**2) > 2*R:
        return 0


    return  (PI/2 - np.arctan2(y1, x1))*180./PI


theta_t = np.zeros(len(global_times))
for i in range(len(theta_t)):
    theta_t[i] = getConfigAngle(soln_x1[i], soln_y1[i], soln_y2[i])


plt.figure()
plt.plot(theta_t)
plt.savefig(IMG_PATH+'theta.png')


# find where it changes
# ipoint = find_overlap(Us, (UTotal-Us)/2)


# Save final state and intersection point
np.save(POST_PATH+'final_state.npy', final_state, allow_pickle=True)

# with open(POST_PATH+'ipoint.pkl', mode='wb') as file:
#     cloudpickle.dump(ipoint, file)


## Plotting
Ustr = (UTotal-Us)/2
fig, (ax1, ax2) = plt.subplots(1,2, figsize=(10,5), constrained_layout=True)
ax1.plot(global_times, Ustr, label='strong', color='blue')
ax1.plot(global_times, Us, label='confused', color='orange')
# if ipoint[-1]:
#     index = int(sum(ipoint[0] > global_times))-1
#     Ustr_interp_func = interp1d(global_times[index:(index+2)],Ustr[index:(index+2)])
#     ax1.plot(ipoint[0], Ustr_interp_func(ipoint[0]), marker="o", color='black')

ax1.set_xlabel('Time')
ax1.set_ylabel('Potential Energy (per contact)')
ax1.legend()

#ax3 = ax2.twinx()
ax2.plot(global_times, -1*stresses[:,1,1], label=r"$\sigma_{yy}$", color='green')
ax2.plot(global_times, -1*stresses[:,0,0], label=r"$\sigma_{xx}$", color='red')
ax2.plot(global_times, stresses[:,0,1], label=r"$\sigma_{xy}$", color='blue')
# if ipoint[-1]:
#     ax2.plot(ipoint[0], np.absolute(ipoint[1]), marker="o", color='black')

#ax3.plot(global_times[:-1], np.diff(-1*stresses[:,1,1]), label=r"$\dot{\sigma}_{yy}$", color='red')
ax2.set_xlabel('Time')
ax2.set_ylabel('Stress Units')
#ax3.set_ylabel('Stress Change with Time')
lines, labels = ax2.get_legend_handles_labels()
#lines2, labels2 = ax3.get_legend_handles_labels()
#ax2.legend(lines + lines2, labels + labels2, loc=4)

plt.savefig(IMG_PATH+'PEnStress.png')





#
