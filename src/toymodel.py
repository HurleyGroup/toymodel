from sympy.physics.vector import dynamicsymbols
from sympy import *
import cloudpickle
import argparse
import os, sys

# read in arguments and important variables
parser = argparse.ArgumentParser()
parser.add_argument("-p", "--postpath", help="specify the post path (output of stress rise)", type=str)
parser.add_argument("-l", "--load", help="specify the load imposed on triplet", type=float)
parser.add_argument("-a", "--ca", help="specify the configuration angle", type=float)
parser.add_argument("-hp", "--pressure", help="specify the hydrostatic pressure", type=float)
args = parser.parse_args()

PI = 3.141592654

KN, KT, KS = float(os.environ['KN']), float(os.environ['KT']), float(os.environ['KS_RISE'])
MU, R, M = float(os.environ['MU']), float(os.environ['R']), float(os.environ['M'])

POST_PATH = args.postpath
LOAD = args.load
THETA = args.ca # in degrees!
P = args.pressure


# load initial conditions [needed for ]
with open(POST_PATH+'init_dof.pkl', 'rb') as ff:
    init_dof = cloudpickle.load(ff)

with open(POST_PATH+'init_offset.pkl', 'rb') as ff:
    init_offset = cloudpickle.load(ff)

x10, y10, y20 = init_dof # somehow find a way to update these


## First establish coordinate system
print('Defining coordinate system.')
e1 = Matrix([1, 0])
e2 = Matrix([0, 1])

## Get initial branch vectors
print('Getting initial conditions.')
d10 = x10*e1 + y10*e2
d20 = -x10*e1 + (y20-y10)*e2
d10hat = d10/sqrt(d10.dot(d10))
d20hat = d20/sqrt(d20.dot(d20))

## create variables of interest
print('Defining variables of interest.')
# degrees of freedom
x1 = dynamicsymbols('x1', real=True)
y1, y2 = dynamicsymbols('y1', real=True), dynamicsymbols('y2', real=True)
x1d = dynamicsymbols('x1', 1, real=True)
y1d, y2d = dynamicsymbols('y1', 1, real=True), dynamicsymbols('y2', 1, real=True)


# material and geometric properties
kn, kt, ks = Float(KN), Float(KT), Float(KS)
R = Float(R)
mu = Float(MU)
m = Float(M)
THETA = Float(THETA*PI/180.)

t = symbols('t')


# loads imposed
F1x, F1y = Float(0.0), Float(0.0)
F2x, F2y = Float(0.0), Float(-1*LOAD)
F1 = Matrix([F1x, F1y])
F2 = Matrix([F2x, F2y])


## Create vectors of interest
print('Defining vectors of interest.')
# normal contact vectors
d1 = x1*e1 + y1*e2
d2 = (-x1)*e1 + (y2-y1)*e2
d3 = y2*e2
d1hat = d1/sqrt(d1.dot(d1))
d2hat = d2/sqrt(d2.dot(d2))
d3hat = d3/sqrt(d3.dot(d3))

# tangential contact vectors
t1 = y1*e1 - x1*e2
t2 = (y2-y1)*e1 + x1*e2
t1hat = t1/sqrt(t1.dot(t1))
t2hat = t2/sqrt(t2.dot(t2))

# define theta
theta = acos(d3hat.dot(d1hat))

## define displacements
print('Finding displacement expressions.')
# normal displacements
deltan1 = d1 - 2*R*d1hat
deltan2 = d2 - 2*R*d2hat

# tangential displacements
deltat1 = (d1-d10) - ((d1-d10).dot(d10hat)*d10hat)
deltat2 = (d2-d20) - ((d2-d20).dot(d20hat)*d20hat)


# lateral support displacements
deltas = (x1 - x10 + init_offset)*e1

## Force expressions
print('Finding monogenic expressions.')
# return spring stiffness accordingly
def getK(k, cond):
    return Piecewise((k, cond >= 0),
                     (0.00001, True), evaluate=False)

kn1 = getK(kn, 2*R - d1.norm())
kn2 = getK(kn, 2*R - d2.norm())
kt1 = getK(kt, 2*R - d1.norm())
kt2 = getK(kt, 2*R - d2.norm())
ks1 = getK(ks, x1-x10 + init_offset)

# normal forces
fn1 = kn1*deltan1
fn2 = kn2*deltan2

# tangential forces
slidecheck1 = kt1*sqrt(deltat1.dot(deltat1)) - mu*sqrt(fn1.dot(fn1))
slidecheck2 = kt2*sqrt(deltat2.dot(deltat2)) - mu*sqrt(fn2.dot(fn2))

sign1 = Piecewise((1, deltat1.dot(t1hat) >= 0),
                  (-1, True))

sign2 = Piecewise((1, deltat2.dot(t2hat) >= 0),
                  (-1, True))

ft1 = Heaviside(-slidecheck1)*kt*deltat1 + sign1*Heaviside(slidecheck1)*mu*sqrt(fn1.dot(fn1))*t1hat
ft2 = Heaviside(-slidecheck2)*kt*deltat2 + sign2*Heaviside(slidecheck2)*mu*sqrt(fn2.dot(fn2))*t2hat


# lateral forces
fs = ks1*deltas # RIGID WALL


# nonconservative forces on top two particles
def Qext(dof):
    return F1.dot(d1.diff(dof)) + F2.dot((d1+d2).diff(dof))


## Kinetic energies
print('Finding kinetic energies.')
T = 0.5*m*(x1d)**2 + 0.5*m*(y1d)**2 + 0.5*m*(y2d)**2


## Potential energies
print('Finding potential energies.')
def f2U(k, f):
    return f.dot(f)/(2.*k)


Un1 = f2U(kn1, fn1)
Un2 = f2U(kn2, fn2)

Ut1 = f2U(kt1, ft1)
Ut2 = f2U(kt2, ft2)

Un = Un1 + Un2
Ut = Ut1 + Ut2
Us = fs.dot(fs)/(2*ks)

U = Un + Ut + Us


## Lagrangian
print('Finding system of equations.')
L = T - U


def lagrangian(dofs):
    assert len(dofs)==2, 'Order must be q and qdot'
    q, qdot = dofs[0], dofs[1]

    return L.diff(qdot).diff(t) - L.diff(q) - Qext(q)

## Obtain system of 4 nonlinear, second-order ODEs
# Actual equations that equal 0
firsteq = lagrangian([x1, x1d])
secondeq = lagrangian([y1, y1d])
thirdeq = lagrangian([y2, y2d])


print('Finding holonomic constraint.')
P2x1 = -R*sin(theta)*e1 - R*cos(theta)*e2
P2x2 = -R*sin(theta)*e1 + R*cos(theta)*e2
P2xs = R*e1

# corrected orientation for middle particle
cf1 = -(fn1+ft1)
cf2 = fn2+ft2
cfs = -fs

g = (Float(P) + (P2x1.dot(cf1) + P2x2.dot(cf2) + P2xs.dot(cfs))/(2*PI*R*R))/m

# C
C1 = g.diff(x1)
C2 = g.diff(y1)
C3 = g.diff(y2)

noddyC = Matrix([C1, C2, C3])
noddyv = -1*(noddyC.dot(noddyC)**-1)*g.diff(t)
noddyv = noddyv.subs(x1d, 0).subs(y1d, 0).subs(y2d, 0) # this ensure that it's a partial derivative wrt time

# we need C dot to compute D dot.
noddyCdot = Matrix([C1.diff(t), C2.diff(t), C3.diff(t)])

# we will find D and D dot as well.
D1 = Matrix([-C2, C1, Float(0.0) ])
D2 = noddyC.cross(D1)
noddyD = Matrix([D1[0], D2[0], D1[1], D2[1], D1[2], D2[2]])
noddyDdot = Matrix([ noddyD[0].diff(t), noddyD[1].diff(t), noddyD[2].diff(t), noddyD[3].diff(t), noddyD[4].diff(t), noddyD[5].diff(t) ])


# function, f, for ODE solvers; (acceleration terms removed!)
print('Finding fns.')
f1 = (firsteq - collect(firsteq, x1d.diff(t),evaluate=False, exact=False)[x1d.diff(t)]*x1d.diff(t))/m
f2 = (secondeq - collect(secondeq, y1d.diff(t),evaluate=False, exact=False)[y1d.diff(t)]*y1d.diff(t))/m
f3 = (thirdeq - collect(thirdeq, y2d.diff(t),evaluate=False, exact=False)[y2d.diff(t)]*y2d.diff(t))/m


# hacky way of fixing sympy's execution fails
def fixMe(f):
    returned = f.subs(Derivative(x1, x1), 1.0)
    returned = returned.subs(Derivative(y1, y1), 1.0)
    returned = returned.subs(Derivative(y2, y2), 1.0)

    return returned


f1 = fixMe(f1)
f2 = fixMe(f2)
f3 = fixMe(f3)
U = fixMe(U)
Us = fixMe(Us)
Ut = fixMe(Ut)
noddyC = fixMe(noddyC)
noddyCdot = fixMe(noddyCdot)
noddyD = fixMe(noddyD)
noddyDdot = fixMe(noddyDdot)
noddyv = fixMe(noddyv)



## lambdify f accordingly
print('Lambdify f1.')
f1_lambda = lambdify((x1, y1, y2, x1d, y1d, y2d), f1, modules=['sympy'])
with open(POST_PATH+'f1.pkl', mode='wb') as file:
   cloudpickle.dump(f1_lambda, file)

print('Lambdify f2.')
f2_lambda = lambdify((x1, y1, y2, x1d, y1d, y2d), f2, modules=['sympy'])
with open(POST_PATH+'f2.pkl', mode='wb') as file:
   cloudpickle.dump(f2_lambda, file)

print('Lambdify f3.')
f3_lambda = lambdify((x1, y1, y2, x1d, y1d, y2d), f3, modules=['sympy'])
with open(POST_PATH+'f3.pkl', mode='wb') as file:
   cloudpickle.dump(f3_lambda, file)

## lambdify PE accordingly
print('Lambdify U.')
U_lambda = lambdify((x1, y1, y2), U, modules=['sympy'])
with open(POST_PATH+'UTotal.pkl', mode='wb') as file:
   cloudpickle.dump(U_lambda, file)


print('Lambdify Us.')
Us_lambda = lambdify((x1, y1, y2), Us, modules=['sympy'])
with open(POST_PATH+'Us.pkl', mode='wb') as file:
   cloudpickle.dump(Us_lambda, file)


print('Lambdify UTan.')
Ut_lambda = lambdify((x1, y1, y2), Ut, modules=['sympy'])
with open(POST_PATH+'Ut.pkl', mode='wb') as file:
   cloudpickle.dump(Ut_lambda, file)

## lambdify contact forces on second particle
print('Lambdify force at contact 1.')
cf1_lambda = lambdify((x1, y1, y2), cf1, modules=['sympy'])
with open(POST_PATH+'cf1.pkl', mode='wb') as file:
   cloudpickle.dump(cf1_lambda, file)

print('Lambdify force at contact 2.')
cf2_lambda = lambdify((x1, y1, y2), cf2, modules=['sympy'])
with open(POST_PATH+'cf2.pkl', mode='wb') as file:
   cloudpickle.dump(cf2_lambda, file)

print('Lambdify lateral force.')
cfs_lambda = lambdify((x1, y1, y2), cfs, modules=['sympy'])
with open(POST_PATH+'cfs.pkl', mode='wb') as file:
   cloudpickle.dump(cfs_lambda, file)

## lambdify constraint stuff
print('Lambdify holonomic C.')
noddyC_lambda = lambdify((x1, y1, y2), noddyC, modules=['sympy'])
with open(POST_PATH+'hC.pkl', mode='wb') as file:
   cloudpickle.dump(noddyC_lambda, file)

print('Lambdify holonomic C dot.')
noddyCdot_lambda = lambdify((x1, y1, y2, x1d, y1d, y2d), noddyCdot, modules=['sympy'])
with open(POST_PATH+'hCdot.pkl', mode='wb') as file:
   cloudpickle.dump(noddyCdot_lambda, file)

print('Lambdify holonomic D.')
noddyD_lambda = lambdify((x1, y1, y2), noddyD, modules=['sympy'])
with open(POST_PATH+'hD.pkl', mode='wb') as file:
   cloudpickle.dump(noddyD_lambda, file)

print('Lambdify holonomic D dot.')
noddyDdot_lambda = lambdify((x1, y1, y2, x1d, y1d, y2d), noddyDdot, modules=['sympy'])
with open(POST_PATH+'hDdot.pkl', mode='wb') as file:
   cloudpickle.dump(noddyDdot_lambda, file)

print('Lambdify holonomic v.')
noddyv_lambda = lambdify((x1, y1, y2), noddyv, modules=['sympy'])
with open(POST_PATH+'hv.pkl', mode='wb') as file:
   cloudpickle.dump(noddyv_lambda, file)




#
