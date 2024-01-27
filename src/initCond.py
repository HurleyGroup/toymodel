from scipy import optimize
from sympy import *
import numpy as np
import cloudpickle
import argparse
import os


# read in arguments and important variables
parser = argparse.ArgumentParser()
parser.add_argument("-p", "--postpath", help="specify the post path (output of stress rise)", type=str)
parser.add_argument("-a", "--ca", help="specify the configuration angle", type=float)
parser.add_argument("-hp", "--pressure", help="specify the hydrostatic pressure", type=float)
args = parser.parse_args()

# given these configurations, we need to get particle positions we have right pressure.
KN = float(os.environ['KN'])
KS = float(os.environ['KS_RISE'])
R = float(os.environ['R'])
POST_PATH = args.postpath
THETA = args.ca # in degrees!
P = args.pressure
PI = 3.141592654

R = Float(R)
THETA = Float(THETA*PI/180.)

## First establish coordinate system
e1 = Matrix([1, 0])
e2 = Matrix([0, 1])

## create variables of interest
# one variable really changes since THETA kept constant. Call this dn
dn = symbols('dn', real=True)

# now write DOF as a function of dn [we assume symmetry of triplet during initial condition!]
x1 = (2*R-dn)*sin(THETA)
y1 = (2*R-dn)*cos(THETA)
y2 = 2*y1

# now get vectors of interest
d1 = x1*e1 + y1*e2
d2 = (-x1)*e1 + (y2-y1)*e2
d1hat = d1/sqrt(d1.dot(d1))
d2hat = d2/sqrt(d2.dot(d2))


## Now compute force as a function of DOF
# First find displacements
deltan1 = d1 - 2*R*d1hat
deltan2 = d2 - 2*R*d2hat

# now get forces
fn1 = -KN*deltan1
fn2 = KN*deltan2
fns = -(fn1.dot(e1) + fn2.dot(e1))*e1 # for balance


# now get locations for forces (relative to middle particle!)
xf1 = -R*sin(THETA)*e1 - R*cos(THETA)*e2
xf2 = -R*sin(THETA)*e1 + R*cos(THETA)*e2
xfs = R*e1


# assemble constraint function
f = lambdify(dn, P + (xf1.dot(fn1) + xf2.dot(fn2) + xfs.dot(fns))/(2*PI*R*R))

## Now find dn for which the hydrostatic constraint will be met
init_dn = optimize.root(f, 1E-4).x[0]


# Finally, find initial conditions!
x10, y10, y20 = x1.subs(dn, init_dn), y1.subs(dn, init_dn), y2.subs(dn, init_dn)
init_dof = [x10, y10, y20]
init_offset = (-fns[0]/KS).subs(dn, init_dn)

# fn1 = fn1.subs(dn, init_dn)
# fn2 = fn2.subs(dn, init_dn)
# fns = fns.subs(dn, init_dn)
# xf1 = xf1.subs(dn, init_dn)
# xf2 = xf2.subs(dn, init_dn)
# xfs = xfs.subs(dn, init_dn)
#
# stress = (np.outer(xf1,fn1) + np.outer(xf2,fn2) + np.outer(xfs,fns))/(2*PI*R*R)
#
# print(np.dot(stress,[0,1])*2*R)
# print(init_offset)
#
# sys.exit()

with open(POST_PATH+'init_dof.pkl', 'wb') as ff:
    cloudpickle.dump(init_dof, ff)

with open(POST_PATH+'init_offset.pkl', 'wb') as ff:
    cloudpickle.dump(init_offset, ff)





#
