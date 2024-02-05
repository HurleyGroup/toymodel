from scipy.integrate._ivp.base import OdeSolver
from tqdm import tqdm

# save the old methods - we still need them
old_init = OdeSolver.__init__
old_step = OdeSolver.step

# define our own methods
def new_init(self, fun, t0, y0, t_bound, vectorized, support_complex=False):

    # define the progress bar
    self.pbar = tqdm(total=t_bound - t0, unit='ut', initial=t0, ascii=True, desc='IVP')
    self.last_t = t0

    # call the old method - we still want to do the old things too!
    old_init(self, fun, t0, y0, t_bound, vectorized, support_complex)


def new_step(self):
    # call the old method
    old_step(self)

    # update the bar
    tst = self.t - self.last_t
    self.pbar.update(tst)
    self.last_t = self.t

    # close the bar if the end is reached
    if self.t >= self.t_bound:
        self.pbar.close()


# overwrite the old methods with our customized ones
OdeSolver.__init__ = new_init
OdeSolver.step = new_step

#################################################################################################

from matplotlib.colors import Normalize, LinearSegmentedColormap
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from sympy import *
import numpy as np
import cloudpickle
import argparse
import os, sys

PLOT = True

# read in arguments and important variables
parser = argparse.ArgumentParser()
parser.add_argument("-mp", "--modelpath", help="specify the model path", type=str)
parser.add_argument("-p", "--postpath", help="specify the post path", type=str)
parser.add_argument("-i", "--imgpath", help="specify the image path", type=str)
parser.add_argument('--maintain', action=argparse.BooleanOptionalAction, help="whether to constrain hydrostatic pressure as constant or not")
args = parser.parse_args()

MODEL_PATH = args.modelpath
POST_PATH = args.postpath
IMG_PATH = args.imgpath
MAINTAIN = args.maintain # If true, then we apply constraint

PI = 3.141592654
R = float(os.environ['R'])

print('Loading Initial Conditions.')
with open(POST_PATH+'init_dof.pkl', 'rb') as ff:
    init_dof = cloudpickle.load(ff)

with open(POST_PATH+'init_offset.pkl', 'rb') as ff:
    init_offset = cloudpickle.load(ff)

x10, y10, y20 = init_dof # somehow find a way to update these
xd0 = np.zeros((3,1))

# Load relevant variables
THETA = float(os.environ['CA'])
LOAD = float(os.environ['LOAD'])
LOAD_TIME = float(os.environ['LOAD_TIME'])
P = float(os.environ['PRESSURE'])

KS = float(os.environ['KS_START'])
KS_RATE_RISE = float(os.environ['KS_RATE_RISE'])
KS_TIME_RISE = float(os.environ['KS_TIME_RISE'])
KS_RATE_DROP = float(os.environ['KS_RATE_DROP'])
KS_TIME_DROP = float(os.environ['KS_TIME_DROP'])

INT_TIME = LOAD_TIME + KS_TIME_RISE + KS_TIME_DROP


print('Loading Functions.')
with open(MODEL_PATH+'f1.pkl', mode='rb') as file:
   f1_lambda = cloudpickle.load(file)

with open(MODEL_PATH+'f2.pkl', mode='rb') as file:
   f2_lambda = cloudpickle.load(file)

with open(MODEL_PATH+'f3.pkl', mode='rb') as file:
   f3_lambda = cloudpickle.load(file)


print('Loading holonomic constraint.')
with open(MODEL_PATH+'hC.pkl', mode='rb') as file:
   hC_lambda = cloudpickle.load(file)

with open(MODEL_PATH+'hCdot.pkl', mode='rb') as file:
   hCdot_lambda = cloudpickle.load(file)

with open(MODEL_PATH+'hD.pkl', mode='rb') as file:
   hD_lambda = cloudpickle.load(file)

with open(MODEL_PATH+'hDdot.pkl', mode='rb') as file:
   hDdot_lambda = cloudpickle.load(file)

with open(MODEL_PATH+'hv.pkl', mode='rb') as file:
   hv_lambda = cloudpickle.load(file)


def hC_lambdaN(t,a,b,c, LOAD, LOAD_TIME, KS, KS_RATE_RISE, KS_TIME_RISE, KS_RATE_DROP, P, x10, y10, y20, init_offset):
    g = hC_lambda(t,a,b,c, LOAD, LOAD_TIME, KS, KS_RATE_RISE, KS_TIME_RISE, KS_RATE_DROP, P, x10, y10, y20, init_offset)
    g1, g2, g3 = g[0], g[1], g[2]
    return np.array([g1, g2, g3]).reshape((1,3)).astype(np.float64)

def hCdot_lambdaN(t,a,b,c,d,e,f, LOAD, LOAD_TIME, KS, KS_RATE_RISE, KS_TIME_RISE, KS_RATE_DROP, P, x10, y10, y20, init_offset):
    g = hCdot_lambda(t,a,b,c,d,e,f, LOAD, LOAD_TIME, KS, KS_RATE_RISE, KS_TIME_RISE, KS_RATE_DROP, P, x10, y10, y20, init_offset)
    g1, g2, g3 = g[0], g[1], g[2]
    return np.array([g1, g2, g3]).reshape((1,3)).astype(np.float64)


def hD_lambdaN(t,a,b,c, LOAD, LOAD_TIME, KS, KS_RATE_RISE, KS_TIME_RISE, KS_RATE_DROP, P, x10, y10, y20, init_offset):
    g = hD_lambda(t,a,b,c, LOAD, LOAD_TIME, KS, KS_RATE_RISE, KS_TIME_RISE, KS_RATE_DROP, P, x10, y10, y20, init_offset)
    g1, g2 = g[0], g[1]
    g3, g4 = g[2], g[3]
    g5, g6 = g[4], g[5]
    return np.array([g1, g2, g3, g4, g5, g6]).reshape((3,2)).astype(np.float64)


def hDdot_lambdaN(t,a,b,c,d,e,f, LOAD, LOAD_TIME, KS, KS_RATE_RISE, KS_TIME_RISE, KS_RATE_DROP, P, x10, y10, y20, init_offset):
    g = hDdot_lambda(t,a,b,c,d,e,f, LOAD, LOAD_TIME, KS, KS_RATE_RISE, KS_TIME_RISE, KS_RATE_DROP, P, x10, y10, y20, init_offset)
    g1, g2 = g[0], g[1]
    g3, g4 = g[2], g[3]
    g5, g6 = g[4], g[5]
    return np.array([g1, g2, g3, g4, g5, g6]).reshape((3,2)).astype(np.float64)

def hv_lambdaN(t,a,b,c, LOAD, LOAD_TIME, KS, KS_RATE_RISE, KS_TIME_RISE, KS_RATE_DROP, P, x10, y10, y20, init_offset):
    g = hv_lambda(t,a,b,c, LOAD, LOAD_TIME, KS, KS_RATE_RISE, KS_TIME_RISE, KS_RATE_DROP, P, x10, y10, y20, init_offset)
    return float(g)


if MAINTAIN:
    D0 = hD_lambdaN(0, x10, y10, y20, LOAD, LOAD_TIME, KS, KS_RATE_RISE, KS_TIME_RISE, KS_RATE_DROP, P, x10, y10, y20, init_offset)
    C0 = hC_lambdaN(0, x10, y10, y20, LOAD, LOAD_TIME, KS, KS_RATE_RISE, KS_TIME_RISE, KS_RATE_DROP, P, x10, y10, y20, init_offset)
    v0 = hC_lambdaN(0, x10, y10, y20, LOAD, LOAD_TIME, KS, KS_RATE_RISE, KS_TIME_RISE, KS_RATE_DROP, P, x10, y10, y20, init_offset)

    u0 = np.linalg.inv(D0.T.dot(D0)).dot(D0.T).dot(xd0 - C0.T*v0)

    state_vector_0 = [x10, y10, y20, u0[0,0], u0[1,0]]
else:
    state_vector_0 = [x10, y10, y20, xd0[0,0], xd0[1,0], xd0[2,0]]


# Now perform solving of system (v is 0)
def helperMAINTAIN(t, vars, LOAD, LOAD_TIME, KS, KS_RATE_RISE, KS_TIME_RISE, KS_RATE_DROP, P, x10, y10, y20, init_offset): # x1, y1, y2, u1, u2
    ## actual stuff
    C = hC_lambdaN(t, vars[0], vars[1], vars[2], LOAD, LOAD_TIME, KS, KS_RATE_RISE, KS_TIME_RISE, KS_RATE_DROP, P, x10, y10, y20, init_offset)
    D = hD_lambdaN(t, vars[0], vars[1], vars[2], LOAD, LOAD_TIME, KS, KS_RATE_RISE, KS_TIME_RISE, KS_RATE_DROP, P, x10, y10, y20, init_offset)
    v = hv_lambdaN(t, vars[0], vars[1], vars[2], LOAD, LOAD_TIME, KS, KS_RATE_RISE, KS_TIME_RISE, KS_RATE_DROP, P, x10, y10, y20, init_offset)
    u = np.array([vars[3], vars[4]]).reshape((2,1)).astype(np.float64)

    xdot = D.dot(u).flatten() + v*C.flatten()

    f1 = f1_lambda(t, vars[0], vars[1], vars[2], xdot[0], xdot[1], xdot[2], LOAD, LOAD_TIME, KS, KS_RATE_RISE, KS_TIME_RISE, KS_RATE_DROP, P, x10, y10, y20, init_offset)
    f2 = f2_lambda(t, vars[0], vars[1], vars[2], xdot[0], xdot[1], xdot[2], LOAD, LOAD_TIME, KS, KS_RATE_RISE, KS_TIME_RISE, KS_RATE_DROP, P, x10, y10, y20, init_offset)
    f3 = f3_lambda(t, vars[0], vars[1], vars[2], xdot[0], xdot[1], xdot[2], LOAD, LOAD_TIME, KS, KS_RATE_RISE, KS_TIME_RISE, KS_RATE_DROP, P, x10, y10, y20, init_offset)
    f = np.array([f1, f2, f3]).reshape((3,1)).astype(np.float64)

    Ddot = hDdot_lambdaN(t, vars[0], vars[1], vars[2], xdot[0], xdot[1], xdot[2], LOAD, LOAD_TIME, KS, KS_RATE_RISE, KS_TIME_RISE, KS_RATE_DROP, P, x10, y10, y20, init_offset)
    Cdot = hCdot_lambdaN(t, vars[0], vars[1], vars[2], xdot[0], xdot[1], xdot[2], LOAD, LOAD_TIME, KS, KS_RATE_RISE, KS_TIME_RISE, KS_RATE_DROP, P, x10, y10, y20, init_offset)

    udot = np.linalg.inv(D.T.dot(D)).dot(D.T).dot( f -  Ddot.dot(u) - Cdot.T*v).flatten()

    return [xdot[0], xdot[1], xdot[2], udot[0], udot[1]]

def helperNOMAINTAIN(t, vars, LOAD, LOAD_TIME, KS, KS_RATE_RISE, KS_TIME_RISE, KS_RATE_DROP, P, x10, y10, y20, init_offset): # x1, y1, y2, x1d, y1d, y2d
    x = [vars[0], vars[1], vars[2]]
    xdot = [vars[3], vars[4], vars[5]]

    f1 = f1_lambda(t, x[0], x[1], x[2], xdot[0], xdot[1], xdot[2], LOAD, LOAD_TIME, KS, KS_RATE_RISE, KS_TIME_RISE, KS_RATE_DROP, P, x10, y10, y20, init_offset)
    f2 = f2_lambda(t, x[0], x[1], x[2], xdot[0], xdot[1], xdot[2], LOAD, LOAD_TIME, KS, KS_RATE_RISE, KS_TIME_RISE, KS_RATE_DROP, P, x10, y10, y20, init_offset)
    f3 = f3_lambda(t, x[0], x[1], x[2], xdot[0], xdot[1], xdot[2], LOAD, LOAD_TIME, KS, KS_RATE_RISE, KS_TIME_RISE, KS_RATE_DROP, P, x10, y10, y20, init_offset)

    return [xdot[0], xdot[1], xdot[2], f1, f2, f3]


print('Integrating.')
start = 0
end = INT_TIME
points = 500
tspan = np.linspace(start, end, points)

if MAINTAIN:
    sol = solve_ivp(helperMAINTAIN, [tspan[0], tspan[-1]],
                        state_vector_0,
                        t_eval=tspan,
                        method='LSODA',
                        jac=None,
                        args=[LOAD, LOAD_TIME, KS, KS_RATE_RISE, KS_TIME_RISE, KS_RATE_DROP, P, x10, y10, y20, init_offset])
else:
    sol = solve_ivp(helperNOMAINTAIN, [tspan[0], tspan[-1]],
                        state_vector_0,
                        t_eval=tspan,
                        args=[LOAD, LOAD_TIME, KS, KS_RATE_RISE, KS_TIME_RISE, KS_RATE_DROP, P, x10, y10, y20, init_offset])



print('Saving Solution.')
soln_x1, soln_y1, soln_y2  = np.zeros(len(sol.t)),  np.zeros(len(sol.t)),  np.zeros(len(sol.t))
soln_vx1, soln_vy1, soln_vy2  = np.zeros(len(sol.t)),  np.zeros(len(sol.t)),  np.zeros(len(sol.t))

print('Positions:')
for t in tqdm(range(len(sol.t))):
    soln_x1[t], soln_y1[t], soln_y2[t]  = sol.y[0][t], sol.y[1][t], sol.y[2][t]

print('Velocities:')
if MAINTAIN:
    for t in tqdm(range(len(sol.t))):
        D = hD_lambdaN(sol.t[t], soln_x1[t], soln_y1[t], soln_y2[t], LOAD, LOAD_TIME, KS, KS_RATE_RISE, KS_TIME_RISE, KS_RATE_DROP, P, x10, y10, y20, init_offset)
        C = hC_lambdaN(sol.t[t], soln_x1[t], soln_y1[t], soln_y2[t], LOAD, LOAD_TIME, KS, KS_RATE_RISE, KS_TIME_RISE, KS_RATE_DROP, P, x10, y10, y20, init_offset)
        v = hv_lambdaN(sol.t[t], soln_x1[t], soln_y1[t], soln_y2[t], LOAD, LOAD_TIME, KS, KS_RATE_RISE, KS_TIME_RISE, KS_RATE_DROP, P, x10, y10, y20, init_offset)
        u = np.array([sol.y[3][t], sol.y[4][t]]).reshape(2,1)

        soln_vx1[t], soln_vy1[t], soln_vy2[t]  = D.dot(u) + v*C.T
else:
    for t in tqdm(range(len(sol.t))):
        soln_vx1[t], soln_vy1[t], soln_vy2[t]  = sol.y[3][t], sol.y[4][t], sol.y[5][t]

soln = [sol.t, soln_x1, soln_y1, soln_y2, soln_vx1, soln_vy1, soln_vy2]


if MAINTAIN:
    with open(POST_PATH+'solMAINTAIN.pkl', mode='wb') as file:
       cloudpickle.dump(soln, file)
else:
    with open(POST_PATH+'solNOMAINTAIN.pkl', mode='wb') as file:
       cloudpickle.dump(soln, file)

if PLOT:
    print('Plotting.')
    cmap = plt.get_cmap('inferno')
    norm = Normalize(vmin=tspan.min(), vmax=tspan.max())

    fig, axes = plt.subplots(3, 2, figsize=(10, 20), constrained_layout=True)
    ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = axes

    # phase plots
    for i in range(len(sol.t) - 1):
        ax1.plot(sol.t[i:i+2], soln_x1[i:i+2], '.-', color=cmap(norm(sol.t[i])))
        ax2.plot(soln_x1[i:i+2], soln_vx1[i:i+2], '.-', color=cmap(norm(sol.t[i])))
        ax3.plot(sol.t[i:i+2], soln_y1[i:i+2], '.-', color=cmap(norm(sol.t[i])))
        ax4.plot(soln_y1[i:i+2], soln_vy1[i:i+2], '.-', color=cmap(norm(sol.t[i])))
        ax5.plot(sol.t[i:i+2], soln_y2[i:i+2], '.-', color=cmap(norm(sol.t[i])))
        ax6.plot(soln_y2[i:i+2], soln_vy2[i:i+2], '.-', color=cmap(norm(sol.t[i])))


    cbar = plt.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), ax=axes.ravel().tolist())

    if MAINTAIN:
        plt.savefig(IMG_PATH+'phasePlotsMAINTAIN.png')
    else:
        plt.savefig(IMG_PATH+'phasePlotsNOMAINTAIN.png')

    # middle particle trajectory
    fig = plt.figure()
    for i in range(len(sol.t) - 1):
        plt.plot(soln_x1[i:i+2], soln_y1[i:i+2], '.-', color=cmap(norm(sol.t[i])))

    cbar = plt.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), ax=fig.get_axes())

    if MAINTAIN:
        plt.savefig(IMG_PATH+'middleTrajectoryMAINTAIN.png')
    else:
        plt.savefig(IMG_PATH+'middleTrajectoryNOMAINTAIN.png')

else:
    print('Skipping Plot.')



#
