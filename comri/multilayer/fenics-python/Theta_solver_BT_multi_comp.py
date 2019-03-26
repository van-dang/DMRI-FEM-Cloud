# This demo solves the Bloch-Torrey equation applied to diffusion MRI using 
# the Nitsche finite element method coupled with the theta-method 
# for the space discretization.

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
# Last changed: 2017-10-13

import sys
import math
from dolfin import *
from mshr import *

import logging; logging.getLogger('UFL').setLevel(logging.ERROR)

import numpy as np
import numpy.linalg as la

# sys.path.append("../")

from parameters_multi_comp import *

set_log_active(False)

dolfin_version = dolfin.dolfin_version()
print 'dolfin version:', dolfin_version

#########################################################
################ Functions
def ThetaMethod_L(ft, gnorm, GX, u0r, u0i, v0r, v0i, u1r, u1i, v1r, v1i, u0r_0, u0i_0, u1r_0, u1i_0,k, kappa, K, theta, phase):
    L0 = (u0r_0/k*v0r + u0i_0/k*v0i+(1-theta)*FuncF(ft, gnorm, GX, u0r_0, u0i_0, v0r, v0i, K))*(1-phase)*dx
    L1 = (u1r_0/k*v1r +u1i_0/k*v1i+(1-theta)*FuncF(ft, gnorm, GX, u1r_0, u1i_0, v1r, v1i, K))*phase*dx
    L_bc = avg((1-theta)*icondition(kappa, u0r_0, u1r_0, v0r, v1r, u0i_0, u1i_0, v0i, v1i))*abs(jump(phase))*dS;
    return L0+L1-L_bc


def ThetaMethod_a(ft, gnorm, GX, u0r, u0i, v0r, v0i, u1r, u1i, v1r, v1i, u0r_0, u0i_0, u1r_0, u1i_0,k,kappa, K, theta, phase):
    a0 = (u0r/k*v0r   + u0i/k*v0i  -theta*FuncF(ft, gnorm, GX, u0r  , u0i  , v0r, v0i, K))*(1-phase)*dx
    a1 = (u1r/k*v1r   + u1i/k*v1i  -theta*FuncF(ft, gnorm, GX, u1r  , u1i  , v1r, v1i, K))*phase*dx
    a_bc = avg(  (theta*icondition(kappa, u0r  , u1r  , v0r, v1r, u0i  , u1i  , v0i, v1i)))*abs(jump(phase))*dS;
    return a0+a1+a_bc


def NoTimeMatrices(u0r, u0i, v0r, v0i, u1r, u1i, v1r, v1i, K, GX, kappa, theta, phase):
    m0 = (u0r*v0r   + u0i*v0i)*(1-phase)*dx
    m1 = (u1r*v1r   + u1i*v1i)*phase*dx
    M = assemble(m0+m1);
        
    j0 = -GX*(u0i*v0r   - u0r*v0i)*(1-phase)*dx
    j1 = -GX*(u1i*v1r   - u1r*v1i)*phase*dx
    J = assemble(j0+j1);    
    s0 = K*( inner(grad(u0r), grad(v0r)) + inner(grad(u0i), grad(v0i)) )*(1-phase)*dx
    s1 = K*( inner(grad(u1r), grad(v1r)) + inner(grad(u1i), grad(v1i)) )*phase*dx
    S = assemble(s0+s1)

    im = avg(  icondition(kappa, u0r  , u1r  , v0r, v1r, u0i  , u1i  , v0i, v1i) )*abs(jump(phase))*dS;
    I = assemble(im)
    
    M.ident_zeros();
    
    return M, J, S, I

def ThetaMethod_A(ft, gnorm, theta, k, M, J, S, I):
    return 1./k*M + ft*gnorm*theta*J + theta*S + theta*I

# Theta method

stepcounter = 0;

theta = 0.5;
Nsteps = 100
k = T/Nsteps;

M, J, S, I = NoTimeMatrices(u0r, u0i, v0r, v0i, u1r, u1i, v1r, v1i, K, GX, kappa, theta, phase);

while t < T: # Time-stepping loop
    if stepcounter % nskip == 0:
        print 't=%f '%t, 'gnorm: %f'%gnorm;

    ft = FT(t, delta, Delta);

    L = ThetaMethod_L(ft, gnorm, GX, u0r, u0i, v0r, v0i, u1r, u1i, v1r, v1i, u0r_0, u0i_0, u1r_0, u1i_0,k, kappa, K, theta, phase);
    # a = ThetaMethod_a(ft, gnorm, GX, u0r, u0i, v0r, v0i, u1r, u1i, v1r, v1i, u0r_0, u0i_0, u1r_0, u1i_0,k,kappa, K, theta, phase);
    # a = a + 1e-6*(u0r*v0r+u0i*v0i)*theta*dx + 1e-6*(u1r*v1r+u1i*v1i)*(1-theta)*dx
    A = ThetaMethod_A(ft, gnorm, theta, k, M, J, S, I);
    # A = assemble(a);
    b = assemble(L);

    # print J.is_symmetric(1e-10);
    # stop
 
    u = Function(W)
    solve(A,u.vector(),b, "gmres", "ilu");
    
    u0r_0, u0i_0, u1r_0, u1i_0 = split(u)
    
    SubMeshSave(u0r_0, u0i_0, file_u0r, file_u0i, mesh0, nskip, stepcounter, dolfin_version);
    SubMeshSave(u1r_0, u1i_0, file_u1r, file_u1i, mesh1, nskip, stepcounter, dolfin_version);
    
    t += k;
    stepcounter += 1;

signal = assemble((phase*u1r_0+(1-phase)*u0r_0)*dx)/assemble(one*dx);
print 'b:',bvalue, 'Signal: ', signal


