from matplotlib.colors import Normalize, LinearSegmentedColormap
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from sympy import *
import numpy as np
import cloudpickle
import argparse
import os, sys

PLOT = False

# read in arguments and important variables
parser = argparse.ArgumentParser()
parser.add_argument("-p", "--postpath", help="specify the post path (output of stress rise)", type=str)
parser.add_argument("-a", "--ca", help="specify the configuration angle", type=float)
parser.add_argument("-d", "--delta", help="specify the initial overlap between particles", type=float)
args = parser.parse_args()

POST_PATH = args.postpath
THETA = args.ca
DELTAN0 = args.delta

PI = 3.141592654

MU, R = float(os.environ['MU']), float(os.environ['R'])
INT_TIME = float(os.environ['INT_TIME'])

print('Loading Functions.')
with open(POST_PATH+'f1.pkl', mode='rb') as file:
   f1_lambda = cloudpickle.load(file)

with open(POST_PATH+'f2.pkl', mode='rb') as file:
   f2_lambda = cloudpickle.load(file)

with open(POST_PATH+'f3.pkl', mode='rb') as file:
   f3_lambda = cloudpickle.load(file)

theta = Float(THETA*PI/180.)

x10 = 2*(R - DELTAN0)*sin(theta)
y10 = 2*(R - DELTAN0)*cos(theta)
y20 = 4*(R - DELTAN0)*cos(theta)
x1d0 = Float(0.0)
y1d0 = Float(0.0)
y2d0 = Float(0.0)

# Now perform solving of system
def helper(t, vars):
    return [vars[3], # velocities
            vars[4],
            vars[5],
            -f1_lambda(vars[0],vars[1],vars[2]), # accelerations
            -f2_lambda(vars[0],vars[1],vars[2]),
            -f3_lambda(vars[0],vars[1],vars[2])]



print('Integrating.')
start = 0
end = INT_TIME
points = 100
tspan = np.linspace(start, end, points)

sol = solve_ivp(helper, [tspan[0], tspan[-1]],
                        [x10, y10, y20, x1d0, y1d0, y2d0],
                        t_eval=tspan)

print('Saving Solution.')
with open(POST_PATH+'sol.pkl', mode='wb') as file:
   cloudpickle.dump(sol, file)

if PLOT:
    print('Plotting.')
    cmap = plt.get_cmap('inferno')
    norm = Normalize(vmin=tspan.min(), vmax=tspan.max())

    fig, axes = plt.subplots(3, 2, figsize=(10, 20), constrained_layout=True)
    ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = axes

    # phase plots
    for i in range(len(sol.t) - 1):
        ax1.plot(sol.t[i:i+2], sol.y[0][i:i+2], '.-', color=cmap(norm(sol.t[i])))
        ax2.plot(sol.y[0][i:i+2], sol.y[3][i:i+2], '.-', color=cmap(norm(sol.t[i])))
        ax3.plot(sol.t[i:i+2], sol.y[1][i:i+2], '.-', color=cmap(norm(sol.t[i])))
        ax4.plot(sol.y[1][i:i+2], sol.y[4][i:i+2], '.-', color=cmap(norm(sol.t[i])))
        ax5.plot(sol.t[i:i+2], sol.y[2][i:i+2], '.-', color=cmap(norm(sol.t[i])))
        ax6.plot(sol.y[2][i:i+2], sol.y[5][i:i+2], '.-', color=cmap(norm(sol.t[i])))


    cbar = plt.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), ax=axes.ravel().tolist())

    plt.savefig('/code/test.png')

    # middle particle trajectory
    fig = plt.figure()
    for i in range(len(sol.t) - 1):
        plt.plot(sol.y[0][i:i+2], sol.y[1][i:i+2], '.-', color=cmap(norm(sol.t[i])))

    cbar = plt.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), ax=fig.get_axes())

    plt.savefig('/code/test1.png')

    # middle particle 3D trajectory phase plot
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    for i in range(len(sol.t) - 1):
        ax.plot(sol.y[0][i:i+2], sol.y[1][i:i+2], np.sqrt(sol.y[3][i:i+2]**2 + sol.y[4][i:i+2]**2), '.-', color=cmap(norm(sol.t[i])))

    plt.savefig('/code/test2.png')
else:
    print('Skipping Plot.')



#
