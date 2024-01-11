from sympy.physics.vector import dynamicsymbols
from sympy.utilities.codegen import codegen
from sympy import *
import cloudpickle
import argparse
import os, sys

# read in arguments and important variables
parser = argparse.ArgumentParser()
parser.add_argument("-p", "--postpath", help="specify the post path (output of stress rise)", type=str)
parser.add_argument("-l", "--load", help="specify the load imposed on triplet", type=float)
parser.add_argument("-a", "--ca", help="specify the configuration angle", type=float)
args = parser.parse_args()

PI = 3.141592654

KN, KT, KS = float(os.environ['KN']), float(os.environ['KT']), float(os.environ['KS'])
MU, R, M = float(os.environ['MU']), float(os.environ['R']), float(os.environ['M'])

POST_PATH = args.postpath
LOAD = args.load
THETA = args.ca # in degrees!


## First establish coordinate system
print('Defining coordinate system.')
e1 = Matrix([1, 0])
e2 = Matrix([0, 1])


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
theta = Float(THETA*PI/180.)

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
d1hat = d1/sqrt(d1.dot(d1))
d2hat = d2/sqrt(d2.dot(d2))

# tangential contact vectors
t1 = y1*e1 - x1*e2
t2 = (y2-y1)*e1 + x1*e2
t1hat = t1/sqrt(t1.dot(t1))
t2hat = t2/sqrt(t2.dot(t2))

## define displacements
print('Finding displacement expressions.')
# normal displacements
deltan1 = d1 - 2*R*d1hat
deltan2 = d2 - 2*R*d2hat

# tangential displacements
deltat1 = sqrt((x1-2*R*sin(theta))**2 + (y1-2*R*cos(theta))**2) * t1hat
deltat2 = sqrt((-2*R*(sin(theta)-sin(theta)))**2 + (y2-2*R*(cos(theta)+cos(theta)))**2) * t2hat


# lateral support displacements
deltas = (x1 - 2*R*sin(theta))*e1

## Force expressions
print('Finding monogenic expressions.')
# normal forces
contactcheck1 = Piecewise( (1, 2*R - sqrt(d1.dot(d1)) > 0),
                           (0, True), evaluate=False)

contactcheck2 = Piecewise( (1, 2*R - sqrt(d2.dot(d2)) > 0),
                            (0, True), evaluate=False)


fn1 = contactcheck1*kn*deltan1
fn2 = contactcheck2*kn*deltan2

# tangential forces
slidecheck1 = sqrt(deltat1.dot(deltat1)) - mu*sqrt(fn1.dot(fn1))/kt
slidecheck2 = sqrt(deltat2.dot(deltat2)) - mu*sqrt(fn2.dot(fn2))/kt

ft1 = contactcheck1*(Heaviside(-slidecheck1)*kt*deltat1 + Heaviside(slidecheck1)*mu*sqrt(fn1.dot(fn1))*t1hat)
ft2 = contactcheck2*(Heaviside(-slidecheck2)*kt*deltat2 + Heaviside(slidecheck2)*mu*sqrt(fn2.dot(fn2))*t2hat)


# lateral forces
fs = ks*deltas # RIGID WALL


# vertical friction
print('Finding polygenic expressions.')
D = mu*sqrt(fs.dot(fs))*y1d



# nonconservative forces on top two particles
def Qext(dof):
    return F1.dot(d1.diff(dof)) + F2.dot((d1+d2).diff(dof))


## Kinetic energies
print('Finding kinetic energies.')
T = 0.5*m*(x1d)**2 + 0.5*m*(y1d)**2 + 0.5*m*(y2d)**2


## Potential energies
print('Finding potential energies.')
Un = (1/(2*kn)) * (fn1.dot(fn1) + fn2.dot(fn2))
Ut = (1/(2*kt)) * (ft1.dot(ft1) + ft2.dot(ft2))
Us = (1/(2*ks)) * (fs.dot(fs))
U = Un + Ut + Us



## Lagrangian
print('Finding system of equations.')
L = T - U


def lagrangian(dofs):
    assert len(dofs)==2, 'Order must be q and qdot'
    q, qdot = dofs[0], dofs[1]

    return L.diff(qdot).diff(t) + D.diff(qdot) - L.diff(q) - Qext(q)

## Obtain system of 4 nonlinear, second-order ODEs
# Actual equations that equal 0
firsteq = lagrangian([x1, x1d])
secondeq = lagrangian([y1, y1d])
thirdeq = lagrangian([y2, y2d])


# function, f, for ODE solvers; (acceleration terms removed!)
print('Finding fns.')
f1 = (firsteq - collect(firsteq, x1d.diff(t),evaluate=False, exact=False)[x1d.diff(t)]*x1d.diff(t))/m
f2 = (secondeq - collect(secondeq, y1d.diff(t),evaluate=False, exact=False)[y1d.diff(t)]*y1d.diff(t))/m
f3 = (thirdeq - collect(thirdeq, y2d.diff(t),evaluate=False, exact=False)[y2d.diff(t)]*y2d.diff(t))/m


# hacky way of fixing sympy's execution fails
f1 = f1.subs(Derivative(x1, x1), 1.0)
f2 = f2.subs(Derivative(x1, x1), 1.0)
f3 = f3.subs(Derivative(x1, x1), 1.0)

f1 = f1.subs(Derivative(y1, y1), 1.0)
f2 = f2.subs(Derivative(y1, y1), 1.0)
f3 = f3.subs(Derivative(y1, y1), 1.0)

f1 = f1.subs(Derivative(y2, y2), 1.0)
f2 = f2.subs(Derivative(y2, y2), 1.0)
f3 = f3.subs(Derivative(y2, y2), 1.0)

## lambdify f accordingly
print('Lambdify f1.')
f1_lambda = lambdify((x1, y1, y2), f1, modules=['sympy'])
with open(POST_PATH+'f1.pkl', mode='wb') as file:
   cloudpickle.dump(f1_lambda, file)

print('Lambdify f2.')
f2_lambda = lambdify((x1, y1, y2), f2, modules=['sympy'])
with open(POST_PATH+'f2.pkl', mode='wb') as file:
   cloudpickle.dump(f2_lambda, file)

print('Lambdify f3.')
f3_lambda = lambdify((x1, y1, y2), f3, modules=['sympy'])
with open(POST_PATH+'f3.pkl', mode='wb') as file:
   cloudpickle.dump(f3_lambda, file)

print('Lambdify U.')
U = U.subs(Derivative(x1, x1), 1.0)
U = U.subs(Derivative(y1, y1), 1.0)
U = U.subs(Derivative(y2, y2), 1.0)

U_lambda = lambdify((x1, y1, y2), U, modules=['sympy'])
with open(POST_PATH+'UTotal.pkl', mode='wb') as file:
   cloudpickle.dump(U_lambda, file)


print('Lambdify Us.')
Us = Us.subs(Derivative(x1, x1), 1.0)
Us = Us.subs(Derivative(y1, y1), 1.0)
Us = Us.subs(Derivative(y2, y2), 1.0)

Us_lambda = lambdify((x1, y1, y2), Us, modules=['sympy'])
with open(POST_PATH+'Us.pkl', mode='wb') as file:
   cloudpickle.dump(Us_lambda, file)

print('Lambdify Ut.')
Ut = Ut.subs(Derivative(x1, x1), 1.0)
Ut = Ut.subs(Derivative(y1, y1), 1.0)
Ut = Ut.subs(Derivative(y2, y2), 1.0)

print('Lambdify UTan.')
Ut_lambda = lambdify((x1, y1, y2), Ut, modules=['sympy'])
with open(POST_PATH+'Ut.pkl', mode='wb') as file:
   cloudpickle.dump(Ut_lambda, file)

print('Lambdify force at contact 1.')
cf1 = fn1+ft1
cf1_lambda = lambdify((x1, y1, y2), cf1, modules=['sympy'])
with open(POST_PATH+'cf1.pkl', mode='wb') as file:
   cloudpickle.dump(cf1_lambda, file)


print('Lambdify force at contact 2.')
cf2 = fn2+ft2
cf2_lambda = lambdify((x1, y1, y2), cf2, modules=['sympy'])
with open(POST_PATH+'cf2.pkl', mode='wb') as file:
   cloudpickle.dump(cf2_lambda, file)

print('Lambdify lateral force.')
cfs =  mu*sqrt(fs.dot(fs))*e2 + fs
cfs_lambda = lambdify((x1, y1, y2), cfs, modules=['sympy'])
with open(POST_PATH+'cfs.pkl', mode='wb') as file:
   cloudpickle.dump(cfs_lambda, file)



#
