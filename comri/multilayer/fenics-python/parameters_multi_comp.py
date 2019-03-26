""" This demo program solves the Bloch-Torrey equation applied to diffusion MRI """

# Copyright (C) 2017 Van-Dang Nguyen (vdnguyen@kth.se)

# This file is part of FEniCS
#
# FEniCS is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import sys
import math
from dolfin import *
from mshr import *

import logging; logging.getLogger('UFL').setLevel(logging.ERROR)

import numpy as np
import numpy.linalg as la

set_log_active(False)
parameters['allow_extrapolation'] = True

dolfin_version = dolfin.dolfin_version()
print 'dolfin version:', dolfin_version

def FuncF(ft, gnorm, GX, ur, ui, vr, vi, K):
    Fr = ft*gnorm*GX*ui*vr - K*inner(grad(ur), grad(vr))
    Fi = - ft*gnorm*GX*ur*vi - K*inner(grad(ui), grad(vi))
    return Fr + Fi

def icondition(kappa, u0rm, u1rm, v0r, v1r, u0im, u1im, v0i, v1i):
    F_bcr = kappa*(u0rm-u1rm)*(v0r-v1r)
    F_bci = kappa*(u0im-u1im)*(v0i-v1i)
    return F_bcr + F_bci

#####################
def SubMeshSave(ur, ui, file_ur, file_ui, mesh, n, stepcounter, dolfin_version):
  if dolfin_version=='1.6.0':
    V = FunctionSpace(mesh, "CG", porder)
  else:
    # For FEniCS 2016, 2017
    Ve = FiniteElement("CG", mesh.ufl_cell(), porder)
    V = FunctionSpace(mesh, Ve)
  if stepcounter % n == 0:
    ur_p = project(ur, V);
    ui_p = project(ui, V);
    ur_p.rename("Real", "label");
    ui_p.rename("Imag", "label");
    file_ur << ur_p;
    file_ui << ui_p;

def ieval(u,omega, phase):
  if omega==1:
    return u('+')*phase('+') + u('-')*phase('-');
  if omega==0:
    return u('+')*(1.-phase('+')) + u('-')*(1.-phase('-'))                              

#################################################################################
# GEOMETRY SETTINGS
xc, yc, zc = 0.0, 0.0, 0.0
mresolution=150

R0 = 5.0;
R = [R0, 3./2.*R0, 2*R0];

domains=[];
for r in R:
    domain0 = Circle(Point(xc, yc), r, int(r*10));
    domains.append(domain0);

for i in xrange(0, len(domains)-1):
    domains[-1].set_subdomain(i+1,domains[i]);

mesh = generate_mesh(domains[-1], mresolution);

V_DG = FunctionSpace(mesh, 'DG', 0)
dofmap_DG = V_DG.dofmap()
phase = Function(V_DG)
cellmarker = CellFunction("size_t", mesh)

for cell in cells(mesh):
    p = cell.midpoint();
    d = p.x()*p.x()+p.y()*p.y()+p.z()*p.z();
    for i in xrange(0, len(R)):
        if d<R[i]*R[i]:
            phase.vector()[dofmap_DG.cell_dofs(cell.index())] = i%2;
            cellmarker[cell.index()] = i%2;
            break;

plot(phase); interactive(); stop;

mesh0 = SubMesh(mesh, cellmarker, 0)
mesh1 = SubMesh(mesh, cellmarker, 1)
# V_DG = FunctionSpace(mesh, 'DG', 0)
# dofmap_DG = V_DG.dofmap()
# phase = Function(V_DG)
# vol = CellVolume(mesh)
# h = 0.5*CellSize(mesh);

# for cell in cells(mesh):
#     phase.vector()[dofmap_DG.cell_dofs(cell.index())] = cellmarker[cell.index()];

mesh_file = File("mesh.xml")
mesh_file << mesh
phase_file = File("Phi.xml")
phase_file << phase
########## END OF GEOMETRY SETTINGS
#################################################################################

# parameters
#####################################################
porder = 1;

delta, Delta = 1000, 50000
t, T = 0, Delta+delta;
K = 3e-3;
# k = 0.1*hmin*hmin/(4*K);
g0, g1 = 1, 0


kappa = 5e-5; # permeability

nskip = 2;

bvalue = 1000;
gnorm = sqrt(bvalue)/sqrt(delta*delta*(Delta-delta/3.0));
#####################################################

#################################################################################
# FUNCTION SPACES
if dolfin_version=='1.6.0':
    V = FunctionSpace(mesh , "CG", porder);
    W = MixedFunctionSpace([V, V, V, V])
else:
    # For FEniCS 2016, 2017
    Ve = FiniteElement("CG", mesh.ufl_cell(), porder)
    TH = MixedElement([Ve,Ve,Ve,Ve])
    V = FunctionSpace(mesh,Ve);
    W = FunctionSpace(mesh, TH)

v = TestFunction(W)
v0r, v0i, v1r, v1i = v[0], v[1], v[2], v[3]

w = TrialFunction(W);
u0r, u0i, u1r, u1i = w[0], w[1], w[2], w[3]
#################################################################################

#################################################################################
# Initial conditions
one = Function(V);
one.vector()[:] = 1;
u_0 = Function(W);
assign(u_0.sub(0), one)
assign(u_0.sub(2), one)
u0r_0, u0i_0, u1r_0, u1i_0 = u_0[0], u_0[1], u_0[2], u_0[3]
#################################################################################

GX=Expression("x[0]*g0+x[1]*g1", g0=g0, g1=g1,domain=mesh,degree=3);

#################################
## output files 
file_u0r = File("results/u0r.pvd")
file_u0i = File("results/u0i.pvd")
file_u1r = File("results/u1r.pvd")
file_u1i = File("results/u1i.pvd")
#################################

def FT(t, delta, Delta):
    ft1 = 1.0*(t>=0 and t<delta) 
    ft2 = -1.0*(t>=Delta and t<=Delta+delta);
    return ft1 + ft2;  
stepcounter = 0;
