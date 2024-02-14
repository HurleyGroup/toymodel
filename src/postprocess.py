from matplotlib.animation import FuncAnimation, FFMpegWriter
from matplotlib import pyplot as plt
import mpl_axes_aligner
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
    writer = FFMpegWriter(fps=120, metadata=dict(artist='Adyota Gupta'), bitrate=1800)

    # Save the animation as an MP4 file
    if MAINTAIN:
        ani.save(IMG_PATH+'/anim_MAINTAIN.mp4', writer=writer)
    else:
        ani.save(IMG_PATH+'/anim_NOMAINTAIN.mp4', writer=writer)

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
    return Ut, Us, Utan

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
        theta = (PI/2 - np.arctan2(y1[t], x1[t]))
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

def herons(s1, s2, s3):
    semi = (s1 + s2 + s3)/2.
    return np.sqrt(semi*(semi-s1)*(semi-s2)*(semi-s3))

def distance(point1, point2):
    # Calculate the Euclidean distance between two points
    return np.sqrt((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2)

def get_normals(x1, x2, x3, ccw):
    # x1, x2, x3 = points of the triangle
    # ccw can be either 1 or -1 to ensure the right normal is obtained (outwards pointing) [ccw means we HAVE to make it ccw]
    # returned is just the returned value (3x2 array -> 3 'a' vectors)
    returned = np.zeros((3,2))

    returned[0,:] = (1-2*ccw)*(x3-x2)
    returned[1,:] = (1-2*ccw)*(x1-x3)
    returned[2,:] = (1-2*ccw)*(x2-x1)

    # pointing outwards, normal (opposite reciprocal)
    returned[0,0], returned[0,1] = returned[0,1], -returned[0,0]
    returned[1,0], returned[1,1] = returned[1,1], -returned[1,0]
    returned[2,0], returned[2,1] = returned[2,1], -returned[2,0]

    # normalize vector in the correct direction
    returned[0,:] /= np.linalg.norm(returned[0,:])
    returned[1,:] /= np.linalg.norm(returned[1,:])
    returned[2,:] /= np.linalg.norm(returned[2,:])

    # compute b, by making its magnitue the length of 'face'
    returned[0,:] *= np.linalg.norm(x3-x2)
    returned[1,:] *= np.linalg.norm(x1-x3)
    returned[2,:] *= np.linalg.norm(x2-x1)

    # compute a, since a = -b/D, where D = 2 for 2D case
    returned[0,:] *= -0.5
    returned[1,:] *= -0.5
    returned[2,:] *= -0.5

    return returned # returns the three a vectors


# obtain particle positions and velocities
times, soln_x1, soln_y1, soln_y2, soln_vx1, soln_vy1, soln_vy2  = soln

def get_mask_before_b(arr, b):
    # Find the index where the array first crosses the threshold value b
    crossing_index = np.where(np.diff(np.sign(arr - b)) != 0)[0]

    if len(crossing_index) > 0:
        # If crossing index is found, return the portion of the array up to the crossing point
        mask = np.arange(len(arr)) <= crossing_index[0]
        return mask
    else:
        # If crossing index is not found, return the entire array
        return np.ones_like(arr, dtype=bool)


if KS_RATE_RISE > 0:
    THETA_LIMIT = 0.
else:
    THETA_LIMIT = 60.

NUM_STEPS = len(times)
angleMask = get_mask_before_b(getConfigAngle(soln_x1, soln_y1, soln_y2), THETA_LIMIT)
times = times[angleMask]
soln_x1 = soln_x1[angleMask]
soln_y1 = soln_y1[angleMask]
soln_y2 = soln_y2[angleMask]
soln_vx1 = soln_vx1[angleMask]
soln_vy1 = soln_vy1[angleMask]
soln_vy2 = soln_vy2[angleMask]
NUM_STEPS = len(times)


# create animation of particles buckling/fluttering
print('Creating animation.')
save_animation([soln_x1, soln_y1, soln_y2])


print('Computing Contact Forces.')
cf1, cf2, cfs =  compute_force(cf1_lambda, cf2_lambda, cfs_lambda,
                            times, soln_x1, soln_y1, soln_y2)


print('Computing Strains.')
# bottom particle, middle particle, top particle
xs = np.vstack([np.zeros(len(soln_x1)), np.zeros(len(soln_y1)), soln_x1, soln_y1, np.zeros(len(soln_y2)), soln_y2]).T

# need to use triangles
cp1 = np.zeros((len(xs), 2))
cp2 = np.zeros((len(xs), 2))
cps = np.zeros((len(xs), 2))
areas = np.zeros(len(xs))

for p in range(len(cp1)):
    dc1 = xs[p,0:2]-xs[p,2:4]
    dc2 = xs[p,4:6]-xs[p,2:4]

    dist1 = np.linalg.norm(xs[p,2:4]-xs[p,0:2])
    dist2 = np.linalg.norm(xs[p,2:4]-xs[p,4:6])

    dc1 /= dist1
    dc2 /= dist2

    cp1[p] = 0.5*dist1*dc1+xs[p,2:4]
    cp2[p] = 0.5*dist2*dc2+xs[p,2:4]
    cps[p] = np.array([R, 0.])+xs[p,2:4]

    areas[p] = herons(distance(cp1[p], cp2[p]), distance(cp1[p], cps[p]), distance(cp2[p], cps[p]))


xs = np.hstack([cp1, cps, cp2])
us = np.diff(xs, axis=0)

d1s = xs[:,2:4]-xs[:,0:2]  # compute distances between point 2 and 1
d2s = xs[:,4:6]-xs[:,0:2]  # """ between point 3 and 1
d3s = xs[:,2:4]-xs[:,4:6]  # """ between point 2 and 3

ccw_flags = np.cross(d1s, d2s)<0 # used to get the correct normal when computing strain using complementary area vector

# compute area of triangle
u_avgs = (us[:,0:2]+us[:,2:4]+us[:,4:6])/3. # average displacement (u_0)

# compute a vectors
avecs = np.zeros((len(xs),3,2))
for t in range(len(avecs)):
    avecs[t] = get_normals(xs[t, 0:2], xs[t, 2:4], xs[t, 4:6], ccw_flags[t])

total_strain = np.zeros((len(us),2,2))
for t in range(len(total_strain)):
    for q in range(3):
        qqa, qqb = int(q*2), int((q+1)*2)
        total_strain[t] += np.outer(us[t, qqa:qqb]-u_avgs[t], avecs[t, q])

    total_strain[t] = total_strain[t]/areas[t]
    total_strain[t] = (total_strain[t]+total_strain[t].T)/2. # symmetric part of deformation gradient is strain

total_strain = np.cumsum(total_strain, axis=0)
oriented_strains = rotate_stress(total_strain, 9) #15)




print('Computing U.')
UTotal, Us, Ut = compute_pe(UTotal_lambda, Us_lambda, Ut_lambda,
                            times, soln_x1, soln_y1, soln_y2)



print('Computing Stresses.')
stresses = compute_stress(cf1, cfs, cf2, soln_x1, soln_y1, soln_y2)
oriented_stresses = rotate_stress(stresses, 9) #)


print('Generating Plots.')
# plot areas over time
plt.figure()
plt.plot(times, areas)
if MAINTAIN:
    plt.savefig(IMG_PATH+'A_MAINTAIN.png')
else:
    plt.savefig(IMG_PATH+'A_NOMAINTAIN.png')

# plot hydrostatic stress over time
plt.figure()
plt.plot(times, 0.5*(oriented_stresses[:,0,0]+oriented_stresses[:,1,1]))
#plt.ylim(-5000,0)
if MAINTAIN:
    plt.savefig(IMG_PATH+'H_MAINTAIN.png')
else:
    plt.savefig(IMG_PATH+'H_NOMAINTAIN.png')

# plot changes in stress energy components over time
plt.figure()
plt.plot(times[:-1], areas[:-1]*oriented_stresses[:-1,0,0]*oriented_strains[:,0,0] - areas[0]*oriented_stresses[0,0,0]*oriented_strains[0,0,0]  )
plt.plot(times[:-1], areas[:-1]*oriented_stresses[:-1,1,1]*oriented_strains[:,1,1] - areas[0]*oriented_stresses[0,1,1]*oriented_strains[0,1,1])
plt.plot(times[:-1], areas[:-1]*oriented_stresses[:-1,0,1]*oriented_strains[:,0,1]+areas[:-1]*oriented_stresses[:-1,1,0]*oriented_strains[:,1,0] - areas[0]*( oriented_stresses[0,0,1]*oriented_strains[0,0,1]+oriented_stresses[0,1,0]*oriented_strains[0,1,0] ) )
if MAINTAIN:
    plt.savefig(IMG_PATH+'StressEnergy_MAINTAIN.png')
else:
    plt.savefig(IMG_PATH+'StressEnergy_NOMAINTAIN.png')

# plot stress components over time
plt.figure()
plt.plot(times, oriented_stresses[:,0,0] - oriented_stresses[0,0,0], label=r"$\Delta\sigma_{xx}$")
plt.plot(times, oriented_stresses[:,1,1] - oriented_stresses[0,1,1], label=r"$\Delta\sigma_{yy}$")
plt.plot(times, oriented_stresses[:,0,1]+oriented_stresses[:,1,0] - (oriented_stresses[0,0,1]+oriented_stresses[0,1,0]), label=r"$2\Delta\sigma_{xy}$")
plt.legend()
if MAINTAIN:
    plt.savefig(IMG_PATH+'Stress_MAINTAIN.png')
else:
    plt.savefig(IMG_PATH+'Stress_NOMAINTAIN.png')


# plot strain components over time
plt.figure()
plt.plot(times[:-1], oriented_strains[:,0,0], label=r"$\epsilon_{xx}$")
plt.plot(times[:-1], oriented_strains[:,1,1], label=r"$\epsilon_{yy}$")
plt.plot(times[:-1], oriented_strains[:,0,1]+oriented_strains[:,1,0], label=r"$2\epsilon_{xy}$")
plt.legend()
if MAINTAIN:
    plt.savefig(IMG_PATH+'Strain_MAINTAIN.png')
else:
    plt.savefig(IMG_PATH+'Strain_NOMAINTAIN.png')


# plot configuration angle over time
theta_t = getConfigAngle(soln_x1, soln_y1, soln_y2)

plt.figure()
plt.plot(times, theta_t)
if MAINTAIN:
    plt.savefig(IMG_PATH+'Theta_MAINTAIN.png')
else:
    plt.savefig(IMG_PATH+'Theta_NOMAINTAIN.png')

# plot potential energies over time
fig, ax1  = plt.subplots(1, figsize=(7,7), constrained_layout=True)
ax2 = ax1.twinx()
Ustr = (UTotal-Us)/2.0
ax2.plot(times, UTotal-UTotal[0], label=r'$\Delta$PE', color='black')
ax2.set_ylabel(r'Change in PE')

ax1.plot(times, oriented_stresses[:,1,1]-oriented_stresses[0,1,1], label=r"$\Delta\sigma_{yy}$", color='green')
ax1.plot(times, oriented_stresses[:,0,0]-oriented_stresses[0,0,0], label=r"$\Delta\sigma_{xx}$", color='red')
ax1.plot(times, oriented_stresses[:,0,1]-oriented_stresses[0,0,1], label=r"$\Delta\sigma_{xy}$", color='blue')
ax1.set_xlabel('Time')
ax1.set_ylabel('Change in Stress Components')

handles1, labels1 = ax1.get_legend_handles_labels()
handles2, labels2 = ax2.get_legend_handles_labels()

# Combine handles and labels
handles = handles1 + handles2
labels = labels1 + labels2

ax1.legend(handles, labels, loc=2)

mpl_axes_aligner.align.yaxes(ax1, 0, ax2, 0, 0.5)

ax1.set_box_aspect(1)
ax2.set_box_aspect(1)

if MAINTAIN:
    plt.savefig(IMG_PATH+'PEnStress_MAINTAIN.png')
else:
    plt.savefig(IMG_PATH+'PEnStress_NOMAINTAIN.png')





#
