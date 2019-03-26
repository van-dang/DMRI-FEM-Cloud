# This demo solves the Bloch-Torrey equation applied to computational diffusion MRI using 
# the finite element method coupled with the theta-method for the spatial discretization.

# Requirements: FEniCS 1.6.0, 2016.x.x, 2017.1.0, 2018.1.0

# The scope of usage: 
# (1) one domain, (2) pure homogeneous Neumann, (3) Allow surface diffusion 

# Execute: python Theta_solver_BT_one_comp.py -N 100 -d 20000 -D 50000 -b 100 -m mesh.xml -v 1 0 0

# Copyright (C) 2017 Van-Dang Nguyen

# This file is part of DOLFIN.

# DOLFIN is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# DOLFIN is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.

# First added:  2017-10-10
# Last changed: 2019-03-16

import math
from dolfin import *
from mshr import *

import logging; logging.getLogger('UFL').setLevel(logging.ERROR)

import numpy as np
import numpy.linalg as la

import getopt, sys

sys.path.append("./")

# set_log_active(False)

dolfin_version = dolfin.__version__
print('dolfin version:', dolfin_version)

## default parameters
mfile = 'mesh.xml';
g0, g1, g2 = 0, 1, 0;
kcoeff = 2.4e-3;
porder = 1;
Nsteps = 100;
bvalue = 1000;
delta = 40000;
Delta = 40000;
k = 100;
gnorm = sqrt(bvalue)/sqrt(delta*delta*(Delta-delta/3.0));
g_ratio = 2.675e8;
qvalue = gnorm/g_ratio*1e12;
nskip = 5;
is_save = 0;
is_dt = 0;
is_input_b = 0;
is_input_q = 0;
## end default parameters

try:
## Input parameters from command lines
    for i in range(0, len(sys.argv)):
        arg = sys.argv[i];
        if arg=='-N':
            Nsteps = int(sys.argv[i+1]);
            print('Nsteps:', Nsteps)
        if arg=='-b':  
            is_input_b = 1;  
            bvalue = float(sys.argv[i+1]);
            print('bvalue:', bvalue)
        if arg=='-q':  
            is_input_q = 1;  
            qvalue = float(sys.argv[i+1]);
            print('qvalue:', bvalue)
        if arg=='-m':
            mfile = sys.argv[i+1];
            print('mesh:',mfile)
        if arg=='-D':
            Delta = float(sys.argv[i+1]);
            print('Delta:', Delta)
        if arg=='-d':
            delta = float(sys.argv[i+1]);
            print('delta:', delta)
        if arg=='-K':
            kcoeff = float(sys.argv[i+1]);
            print('diffusion coefficient:', kcoeff)
        if arg=='-k':
            is_dt = True;
            k = float(sys.argv[i+1]);
            print('time step size:', k)
        if arg=='-v':
            g0 = float(sys.argv[i+1]);
            g1 = float(sys.argv[i+2]);
            g2 = float(sys.argv[i+3]);
            print('(g0, g1, g2):',g0, g1, g2)
except:
    print('Nothing!')

if (Delta-delta/3.0<=0):
    print('Check time sequence!');
    exit(0);



if (is_input_b):
   gnorm = sqrt(bvalue)/sqrt(delta*delta*(Delta-delta/3.0));
   qvalue = gnorm/g_ratio*1e12;
elif (is_input_q):
   gnorm = qvalue*g_ratio*1e-12;
   bvalue = gnorm*gnorm*delta*delta*(Delta-delta/3.0);

mesh = Mesh(mfile);
gdim = mesh.geometry().dim()

t, T = 0, Delta+delta;

if (is_dt==False):
    k = T/Nsteps;

## FUNCTION SPACES
if dolfin_version=='1.6.0':
    V = FunctionSpace(mesh , "CG", porder); # order 1, 2 components
    W = MixedFunctionSpace([V, V])
else:
    # For FEniCS 2016, 2017
    Ve = FiniteElement("CG", mesh.ufl_cell(), porder)
    V = FunctionSpace(mesh,Ve);
    TH = MixedElement([Ve,Ve])
    W = FunctionSpace(mesh, TH)

v = TestFunction(W)
v1r, v1i = v[0], v[1]

u = TrialFunction(W);
u1r, u1i = u[0], u[1]

# Initial conditions
if (gdim==2):
  Dirac_Delta = Expression("x[0]*x[0]+x[1]*x[1]<eps",eps=1e6, domain=mesh, degree=2);
if (gdim==3):
  Dirac_Delta = Expression("x[0]*x[0]+x[1]*x[1]+x[2]*x[2]<eps",eps=1e6, domain=mesh, degree=2);
Dirac_Delta = interpolate(Dirac_Delta, V);

u_0 = Function(W);
assign(u_0.sub(0), Dirac_Delta)
u1r_0, u1i_0 = split(u_0)

if (gdim==2):
    GX=Expression("x[0]*g0+x[1]*g1", g0=g0, g1=g1, domain=mesh, degree=3);
if (gdim==3):
    GX=Expression("x[0]*g0+x[1]*g1+x[2]*g2", g0=g0, g1=g1, g2=g2, domain=mesh, degree=3);

def FT(t, delta, Delta):
    ft1 = 1.0*(t>=0 and t<delta) 
    ft2 = -1.0*(t>=Delta and t<=Delta+delta);
    return ft1 + ft2;  

def iFT(t, delta, Delta): # integrate ft
    ft1 = t*(t>=0 and t<delta) 
    ft2 = delta*(t>=delta and t<Delta) 
    ft3 = (delta - t + Delta)*(t>=Delta and t<=Delta+delta) 
    return ft1 + ft2 + ft3;  
stepcounter = 0;

## Theta method
def FuncF(ft, gnorm, GX, ur, ui, vr, vi, kcoeff):
    Fr = ft*gnorm*GX*ui*vr - kcoeff*inner(grad(ur), grad(vr))
    Fi = - ft*gnorm*GX*ur*vi - kcoeff*inner(grad(ui), grad(vi))
    return Fr + Fi

def ThetaMethod_L(ft, gnorm, GX, u1r_0, u1i_0, u1r, u1i, v1r, v1i,k, kcoeff, theta):
    L1 = (u1r_0/k*v1r +u1i_0/k*v1i+(1-theta)*FuncF(ft, gnorm, GX, u1r_0, u1i_0, v1r, v1i, kcoeff))*dx
    return L1

def ThetaMethod_a(ft, gnorm, GX, u1r, u1i, v1r, v1i,k, kcoeff, theta):
    a1 = (u1r/k*v1r   + u1i/k*v1i  -theta*FuncF(ft, gnorm, GX, u1r  , u1i  , v1r, v1i, kcoeff))*dx
    return a1

def NoTimeMatrices(u1r, u1i, v1r, v1i, kcoeff, GX, theta):
    m1 = (u1r*v1r   + u1i*v1i)*dx
    M = assemble(m1);
    j1 = -GX*(u1i*v1r   - u1r*v1i)*dx
    J = assemble(j1);    
    s1 = kcoeff*( inner(grad(u1r), grad(v1r)) + inner(grad(u1i), grad(v1i)) )*dx
    S = assemble(s1)
    M.ident_zeros();    
    return M, J, S

def ThetaMethod_A(ft, gnorm, theta, k, M, J, S):
    return 1./k*M + ft*gnorm*theta*J + theta*S

theta = 0.5;

hmin = mesh.hmin();

print('mesh.hmin: ', hmin,'mesh.hmax: ',mesh.hmax());

print("\ndiffusion coefficient: %e, Nsteps: %d, dt: %f, delta: %f, Delta: %f, gnorm: %f\n"%(kcoeff, Nsteps, k, delta, Delta, gnorm));

print("Gradient direction: %f %f %f\n"%(g0, g1, g2));
  
M, J, S = NoTimeMatrices(u1r, u1i, v1r, v1i, kcoeff, GX, theta);

stepcounter = 0;

solver = KrylovSolver('gmres', 'ilu')
solver.parameters["absolute_tolerance"] = 1e-6
solver.parameters["relative_tolerance"] = 1e-4
solver.parameters["maximum_iterations"] = 1000

while t < T + k: # Time-stepping loop
    if stepcounter % nskip == 0:
        print('t: %f '%t, 'T: %f'%T, 'dt: %f'%k,'gnorm: %f'%gnorm,'Completed %.2f%%'%(float(t)/float(T+k)*100.0));
    ft = FT(t, delta, Delta);
    ift = iFT(t, delta, Delta);
    L = ThetaMethod_L(ft, gnorm, GX,  u1r_0, u1i_0, u1r, u1i, v1r, v1i, k, kcoeff, theta);
    A = ThetaMethod_A(ft, gnorm, theta, k, M, J, S);
    b = assemble(L);
    u = Function(W)
    solver.solve(A,u.vector(),b);
    u1r_0, u1i_0 = split(u)
    t += k;
    stepcounter += 1;

signal = assemble(u1r_0*dx)/assemble(Dirac_Delta*dx);
print('b:',bvalue, 'Signal: %.3e'%signal,', dt:',k,', hmin:',hmin)
print("b: %f, gnorm: %f, q: %f, gdir: (%f, %f, %f), s: %f\n"%(bvalue, gnorm, qvalue, g0, g1, g2, signal));

