# This is a collection of functions and classes used to solve the Bloch-Torrey equation 
# applied to computational diffusion MRI using the finite element method coupled with 
# the theta-method for the spatial discretization.

# Copyright (C) 2019 Van-Dang Nguyen (vdnguyen@kth.se or dang.1032170@gmail.com)

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
# Last changed: 2019-04-25

# This demo is maintained by Van-Dang Nguyen
# Please report possible problems to vdnguyen@kth.se or dang.1032170@gmail.com

from dolfin import *
import sympy as sp

import time, os, sys, shutil, mpi4py, numpy

def GdotX(gdir, mymesh):
  gdim = mymesh.geometry().dim()
  if (gdim==2):
    GX=Expression("x[0]*g0+x[1]*g1", g0=gdir.x(), g1=gdir.y(), domain=mymesh, degree=3);
  if (gdim==3):
    GX=Expression("x[0]*g0+x[1]*g1+x[2]*g2", g0=gdir.x(), g1=gdir.y(), g2=gdir.z(), domain=mymesh, degree=3);
  return GX;

def FuncF_wBC(ft, qvalue, gdir, ur, ui, vr, vi, D, mymesh, T2):
    GX=GdotX(gdir, mymesh)
    Fr =   ft*qvalue*GX*ui*vr - 1./T2*ur*vr - inner(D*grad(ur), grad(vr))
    Fi = - ft*qvalue*GX*ur*vi - 1./T2*ui*vi - inner(D*grad(ui), grad(vi))
    return Fr + Fi
  
def icondition_wBC(kappa, u0rm, u1rm, v0r, v1r, u0im, u1im, v0i, v1i):
    F_bcr = kappa*(u0rm-u1rm)*(v0r-v1r)
    F_bci = kappa*(u0im-u1im)*(v0i-v1i)
    return F_bcr + F_bci

def ieval(u,omega, phase):
  if omega==1:
    return u('+')*phase('+') + u('-')*phase('-');
  if omega==0:
    return u('+')*(1.-phase('+')) + u('-')*(1.-phase('-'))                      
  
def ThetaMethodL_wBC1c(ft, ift, mri_para , w , v, u_0, sp, mydomain, wpperiodic):
    D = mydomain.D
    T2 = mri_para.T2    
    mymesh = mydomain.mymesh
    k = sp.k
    theta = sp.theta
    qvalue = mri_para.qvalue
    gdir = mri_para.gdir
    u0r_0, u0i_0 = split(u_0)
    v0r, v0i = v[0], v[1]
    u0r, u0i = w[0], w[1]
    L0 = (u0r_0/k*v0r + u0i_0/k*v0i+(1-theta)*FuncF_wBC(ft, qvalue, gdir, u0r_0, u0i_0, v0r, v0i, D, mymesh, T2))*dx
    L_pbc = 0;
    if sum(mydomain.PeriodicDir)>0:
      # Start applying the weak pseudo-periodic BC
      wpperiodic.set_values(u_0, ift);
      u0r_bc, u0i_bc = wpperiodic[0], wpperiodic[1]
      L_pbc +=  (1-theta)*mydomain.kappa_e*(u0r_bc*v0r   + u0i_bc*v0i)*ds; # u0r_0, u0i_0                                                          
      # End of applying the weak pseudo-periodic BC
    return L0+L_pbc

def ThetaMethodF_wBC1c(ft, ift, mri_para, w , v, sp, mydomain):
    D = mydomain.D
    mymesh = mydomain.mymesh
    k = sp.k
    theta = sp.theta
    qvalue = mri_para.qvalue
    gdir = mri_para.gdir
    T2 = mri_para.T2
    v0r, v0i = v[0], v[1]
    u0r, u0i = w[0], w[1]
    a0 = (  -theta*FuncF_wBC(ft, qvalue, gdir, u0r, u0i, v0r, v0i, D, mymesh, T2))*dx
    a_pbc = 0;
    if sum(mydomain.PeriodicDir)>0:
      a_pbc = theta*mydomain.kappa_e*(u0r*v0r   + u0i*v0i)*ds;
    return a0+a_pbc
  
  
def ThetaMethodL_wBC2c(ft, ift, mri_para , w , v, u_0, sp, mydomain, wpperiodic):
    D = mydomain.D
    mymesh = mydomain.mymesh
    k = sp.k
    theta = sp.theta
    qvalue = mri_para.qvalue
    gdir = mri_para.gdir
    T2 = mri_para.T2
    kappa = mydomain.kappa
    phase = mydomain.phase
    u0r_0, u0i_0, u1r_0, u1i_0 = split(u_0)
    v0r, v0i, v1r, v1i = v[0], v[1], v[2], v[3]
    u0r, u0i, u1r, u1i = w[0], w[1], w[2], w[3]

    L0 = (u0r_0/k*v0r + u0i_0/k*v0i+(1-theta)*FuncF_wBC(ft, qvalue, gdir, u0r_0, u0i_0, v0r, v0i, D, mymesh, T2))*(1-phase)*dx
    L1 = (u1r_0/k*v1r +u1i_0/k*v1i+(1-theta)*FuncF_wBC(ft, qvalue, gdir, u1r_0, u1i_0, v1r, v1i, D, mymesh, T2))*phase*dx
    L_bc = avg((1-theta)*icondition_wBC(kappa, u0r_0, u1r_0, v0r, v1r, u0i_0, u1i_0, v0i, v1i))*abs(jump(phase))*dS;
    
    L_pbc = 0;
    if sum(mydomain.PeriodicDir)>0:
      # Start applying the weak pseudo-periodic BC
      wpperiodic.set_values(u_0, ift);
      u0r_bc, u0i_bc, u1r_bc, u1i_bc = wpperiodic[0], wpperiodic[1], wpperiodic[2], wpperiodic[3]
      L_pbc +=  (1-theta)*mydomain.kappa_e*(u1r_bc*v1r   + u1i_bc*v1i)*    phase*ds; # u1r_0, u1i_0                                                          
      L_pbc +=  (1-theta)*mydomain.kappa_e*(u0r_bc*v0r   + u0i_bc*v0i)*(1-phase)*ds; # u0r_0, u0i_0                                                          
      # End of applying the weak pseudo-periodic BC
    return L0+L1-L_bc+L_pbc
  
  
def ThetaMethodF_wBC2c(ft, ift, mri_para, w , v, sp, mydomain):
    D = mydomain.D
    mymesh = mydomain.mymesh
    k = sp.k
    theta = sp.theta  
    qvalue = mri_para.qvalue
    gdir = mri_para.gdir
    kappa = mydomain.kappa
    phase = mydomain.phase
    T2 = mri_para.T2
    v0r, v0i, v1r, v1i = v[0], v[1], v[2], v[3]
    u0r, u0i, u1r, u1i = w[0], w[1], w[2], w[3]
    a0 = (  -theta*FuncF_wBC(ft, qvalue, gdir, u0r  , u0i  , v0r, v0i, D, mymesh, T2))*(1-phase)*dx
    a1 = (  -theta*FuncF_wBC(ft, qvalue, gdir, u1r  , u1i  , v1r, v1i, D, mymesh, T2))*phase*dx
    a_bc = avg(  (theta*icondition_wBC(kappa, u0r  , u1r  , v0r, v1r, u0i  , u1i  , v0i, v1i)))*abs(jump(phase))*dS;
    
    a_pbc = 0;
    if sum(mydomain.PeriodicDir)>0:
      a_pbc = theta*mydomain.kappa_e*(u1r*v1r   + u1i*v1i)*phase*ds + theta*mydomain.kappa_e*(u0r*v0r   + u0i*v0i)*(1-phase)*ds;
      
    return a0+a1+a_bc+a_pbc


# Strong periodic boundary conditions  
def inner_interface(ift, kappa, qvalue, u0rm, u1rm, v0r, v1r, u0im, u1im, v0i, v1i, fn, g, D0, D1):
    F_bcr  = (-kappa*avg(u0rm-u1rm)-0.5*qvalue*ift*(avg(u0im)*inner(avg(D0)*avg(g),fn)+avg(u1im)*inner(avg(D1)*avg(g),fn)))*avg(v0r-v1r)                      
    F_bcr += -qvalue*ift*( avg(u0im)*inner(avg(D0)*avg(g),fn)-avg(u1im)*inner(avg(D1)*avg(g),fn) )*0.5*avg(v0r+v1r)                                                                                                                                                                                         
    F_bci  = (-kappa*avg(u0im-u1im)+0.5*qvalue*ift*(avg(u0rm)*inner(avg(D0)*avg(g),fn)+avg(u1rm)*inner(avg(D0)*avg(g),fn)))*avg(v0i-v1i)                      
    F_bci += qvalue*ift*(  avg(u0rm)*inner(avg(D0)*avg(g),fn)-avg(u1rm)*inner(avg(D1)*avg(g),fn) )*0.5*avg(v0i+v1i)                                         
    return -F_bcr - F_bci

def outer_interface(ift, qvalue, D, fn, ur, ui, vr, vi, g):
    F_bcr =  (ift*qvalue+1e-16)*inner(D*g, fn)*ui*vr
    F_bci = -(ift*qvalue+1e-16)*inner(D*g, fn)*ur*vi
    return F_bcr + F_bci

  
def FuncF_sBC(ift, qvalue, g, ur, ui, vr, vi, D, T2):
    Fr =   ift*qvalue*(inner(g,D*grad(ui))+inner(grad(ui),D*g))*vr - inner(g,D*g)*qvalue*qvalue*ift*ift*ur*vr-1./T2*ur*vr-inner(D*grad(ur), grad(vr))
    Fi = - ift*qvalue*(inner(g,D*grad(ur))+inner(grad(ur),D*g))*vi - inner(g,D*g)*qvalue*qvalue*ift*ift*ui*vi-1./T2*ui*vi-inner(D*grad(ui), grad(vi))
    return Fr + Fi 
  
def ThetaMethodF_sBC1c(ft, ift, mri_para, w, v, sp, mydomain):
    D = mydomain.D  
    T2 = mri_para.T2  
    k = sp.k
    theta = sp.theta  
    qvalue = mri_para.qvalue
    g = mri_para.g
    fn = mydomain.fn
    v0r, v0i = v[0], v[1]
    u0r, u0i = w[0], w[1]
    a0 = (  -theta*FuncF_sBC(ift, qvalue, g, u0r  , u0i  , v0r, v0i, D, T2))*dx
    a0_outer_bc = theta*outer_interface(ift, qvalue , D, fn, u0r, u0i, v0r, v0i, g)*ds
    return a0 + a0_outer_bc 

def ThetaMethodL_sBC1c(ft, ift, mri_para, w, v, u_0, sp, mydomain):
    D = mydomain.D 
    T2 = mri_para.T2      
    k = sp.k
    theta = sp.theta  
    qvalue = mri_para.qvalue
    g = mri_para.g
    fn = mydomain.fn
    v0r, v0i = v[0], v[1]
    u0r, u0i = w[0], w[1]
    u0r_0, u0i_0 = split(u_0)
    L0 = (u0r_0/k*v0r + u0i_0/k*v0i +theta*FuncF_sBC(ift, qvalue, g, u0r_0, u0i_0, v0r, v0i, D, T2))*dx
    L0_outer_bc = -theta*outer_interface(ift, qvalue, D, fn, u0r_0, u0i_0, v0r, v0i, g)*ds
    return L0 + L0_outer_bc

  
def ThetaMethodF_sBC2c(ft, ift, mri_para, w, v, sp, mydomain):
    D = mydomain.D 
    T2 = mri_para.T2      
    k = sp.k
    theta = sp.theta  
    qvalue = mri_para.qvalue
    g = mri_para.g
    fn = mydomain.fn
    fn0 = mydomain.fn0
    kappa = mydomain.kappa
    phase = mydomain.phase
    v0r, v0i, v1r, v1i = v[0], v[1], v[2], v[3]
    u0r, u0i, u1r, u1i = w[0], w[1], w[2], w[3]
    a0 = (  -theta*FuncF_sBC(ift, qvalue, g, u0r  , u0i  , v0r, v0i, D, T2))*(1-phase)*dx
    a1 = (  -theta*FuncF_sBC(ift, qvalue, g, u1r  , u1i  , v1r, v1i, D, T2))*phase*dx
    a_inner_bc  = (  (theta*inner_interface(ift, kappa, qvalue, u0r  , u1r  , v0r, v1r, u0i  , u1i  , v0i, v1i, fn0, g, D, D)))*abs(jump(phase))*dS;
    a0_outer_bc = theta*outer_interface(ift, qvalue , D, fn, u0r, u0i, v0r, v0i, g)*ds
    a1_outer_bc = theta*outer_interface(ift, qvalue , D, fn, u1r, u1i, v1r, v1i, g)*ds
    return a0+a1 +a_inner_bc + a0_outer_bc + a1_outer_bc

  
def ThetaMethodL_sBC2c(ft, ift, mri_para, w, v, u_0, sp, mydomain):
    D = mydomain.D 
    T2 = mri_para.T2      
    k = sp.k
    theta = sp.theta  
    qvalue = mri_para.qvalue
    g = mri_para.g
    fn = mydomain.fn
    fn0 = mydomain.fn0
    kappa = mydomain.kappa
    phase = mydomain.phase
    v0r, v0i, v1r, v1i = v[0], v[1], v[2], v[3]
    u0r, u0i, u1r, u1i = w[0], w[1], w[2], w[3]
    u0r_0, u0i_0, u1r_0, u1i_0 = split(u_0)

    L0 = (u0r_0/k*v0r + u0i_0/k*v0i +theta*FuncF_sBC(ift, qvalue, g, u0r_0, u0i_0, v0r, v0i, D, T2))*(1-phase)*dx
    L1 = (u1r_0/k*v1r + u1i_0/k*v1i +theta*FuncF_sBC(ift, qvalue, g, u1r_0, u1i_0, v1r, v1i, D, T2))*phase*dx
    L_inner_bc  = -(theta*inner_interface(ift, kappa, qvalue, u0r_0, u1r_0, v0r, v1r, u0i_0, u1i_0, v0i, v1i, fn0, g, D, D))*abs(jump(phase))*dS;
    L0_outer_bc = -theta*outer_interface(ift, qvalue, D, fn, u0r_0, u0i_0, v0r, v0i, g)*ds
    L1_outer_bc = -theta*outer_interface(ift, qvalue, D, fn, u1r_0, u1i_0, v1r, v1i, g)*ds
    return L0+L1+L_inner_bc + L0_outer_bc+L1_outer_bc
    
def MassMatrix2c(w, v, phase):
    v0r, v0i, v1r, v1i = v[0], v[1], v[2], v[3]
    u0r, u0i, u1r, u1i = w[0], w[1], w[2], w[3]
    m0 = (u0r*v0r   + u0i*v0i)*(1-phase)*dx
    m1 = (u1r*v1r   + u1i*v1i)*phase*dx
    M = assemble(m0+m1);
    M.ident_zeros()
    return M;
  
def MassMatrix1c(w, v):
    v1r, v1i = v[0], v[1]
    u1r, u1i = w[0], w[1]
    M = assemble( (u1r*v1r   + u1i*v1i)*dx );
    M.ident_zeros()
    return M;
  
class WeakPseudoPeriodic_2c(UserExpression):
    def __init__(self, mydomain, **kwargs):
        self.xmin, self.ymin, self.zmin, self.xmax, self.ymax, self.zmax = mydomain.xmin, mydomain.ymin, mydomain.zmin, mydomain.xmax, mydomain.ymax, mydomain.zmax 
        self.pdir = mydomain.PeriodicDir
        self.gdim = mydomain.gdim
        self.gdir = mydomain.gdir        
        self.qvalue = mydomain.qvalue        
        super().__init__(**kwargs)
        

    def set_values(self, u, ift):
        u.set_allow_extrapolation(True)
        self.u = u
        self.ift = ift        
    def eval(self, value, x):
        tol = 1E-7
        is_eval = False;
        xln, yln, zln = 0, 0, 0
        temp_u0r, temp_u0i, temp_u1r, temp_u1i = 0, 0, 0, 0
                
        if (self.pdir[0] == 1): # x-direction
            if abs(x[0]-self.xmin) <= tol:
                xln, yln, zln = self.xmax-x[0], 0, 0;
                if self.gdim ==2:
                  temp_u0r, temp_u0i, temp_u1r, temp_u1i = self.u([self.xmax, x[1]]);
                if self.gdim ==3:
                  temp_u0r, temp_u0i, temp_u1r, temp_u1i = self.u([self.xmax, x[1], x[2]]);

            if abs(x[0]-self.xmax) <= tol:
                xln, yln, zln = self.xmin-x[0], 0, 0;
                if (self.gdim==2):
                  temp_u0r, temp_u0i, temp_u1r, temp_u1i = self.u([self.xmin, x[1]]);
                if (self.gdim==3):
                  temp_u0r, temp_u0i, temp_u1r, temp_u1i = self.u([self.xmin, x[1], x[2]]);

        if (self.pdir[1] == 1): # y-direction
            if abs(x[1]-self.ymin) <= tol:
                xln, yln, zln = 0, self.ymax-x[1], 0;
                if self.gdim ==2:
                  temp_u0r, temp_u0i, temp_u1r, temp_u1i = self.u([x[0],self.ymax]);
                if self.gdim ==3:
                  temp_u0r, temp_u0i, temp_u1r, temp_u1i = self.u([x[0], self.ymax, x[2]]);

            if abs(x[1]-self.ymax) <= tol:
                xln, yln, zln = 0, self.ymin-x[1], 0;
                if (self.gdim==2):
                  temp_u0r, temp_u0i, temp_u1r, temp_u1i = self.u([x[0], self.ymin]);
                if (self.gdim==3):
                  temp_u0r, temp_u0i, temp_u1r, temp_u1i = self.u([x[0], self.ymin, x[2]]);

        if (self.pdir[2] == 1 and self.gdim==3): # z-direction
            if abs(x[2]-self.zmin) <= tol:
                xln, yln, zln = 0, 0, self.zmax-x[2];
                temp_u0r, temp_u0i, temp_u1r, temp_u1i = self.u([x[0], x[1], self.zmax]);

            if abs(x[2]-self.zmax) <= tol:
                xln, yln, zln = 0, 0, self.zmin-x[2];
                temp_u0r, temp_u0i, temp_u1r, temp_u1i = self.u([x[0], x[1], self.zmin]);

        # exp(i*theta)= cos(theta)+i*sin(theta); 
        theta_ln = self.qvalue*(self.gdir.x()*xln + self.gdir.y()*yln + self.gdir.z()*zln)*self.ift;
        value[0] = temp_u0r*cos(theta_ln)-temp_u0i*sin(theta_ln); # for Real part
        value[1] = temp_u0r*sin(theta_ln)+temp_u0i*cos(theta_ln); # for Imag part

        value[2] = temp_u1r*cos(theta_ln)-temp_u1i*sin(theta_ln); # for Real part
        value[3] = temp_u1r*sin(theta_ln)+temp_u1i*cos(theta_ln); # for Imag part

    def value_shape(self):
        return (4,)
      
        
class PeriodicBD(SubDomain):

    def __init__(self, mydomain, **kwargs):
        self.xmin, self.ymin, self.zmin = mydomain.xmin, mydomain.ymin, mydomain.zmin
        self.xmax, self.ymax, self.zmax = mydomain.xmax, mydomain.ymax, mydomain.zmax
        self.gdim = mydomain.gdim
        self.PeriodicDir = mydomain.PeriodicDir
        self.tol = 1e-2*mydomain.mymesh.hmin()
        super().__init__(**kwargs)
    
    # Left boundary is "target domain" G
    def inside(self, x, on_boundary):
        # return True if on left or bottom boundary AND NOT on one of the two corners (0, 1) and (1, 0) 
        bcx = abs(x[0] - self.xmin)<self.tol and self.PeriodicDir[0] == 1
        bcy = abs(x[1] - self.ymin)<self.tol and self.PeriodicDir[1] == 1
        bcz = 0
        if (self.gdim==3):
          bcz = abs(x[2] - self.zmin) < self.tol and self.PeriodicDir[2] == 1
          
        return bool((bcx or bcy or bcz) and on_boundary)

    def map(self, x, y):
        Lxx, Lyy, Lzz = 0, 0, 0
        if self.PeriodicDir[0] == 0:
          Lxx = 1e7;
        if self.PeriodicDir[1] == 0:
          Lyy = 1e7;
        if self.PeriodicDir[2] == 0:
          Lzz = 1e7;
          
        if abs(x[0] - self.xmax) < self.tol and self.PeriodicDir[0] == 1:
            y[0] = x[0] - (self.xmax - self.xmin) + Lxx
            y[1] = x[1]
            if (self.gdim==3):
              y[2] = x[2]
        elif abs(x[1] - self.ymax) < self.tol and self.PeriodicDir[1] == 1:
            y[0] = x[0]
            y[1] = x[1] - (self.ymax-self.ymin) + Lyy 
            if (self.gdim==3):
              y[2] = x[2]
        elif (self.gdim==3) and abs(x[2] - self.zmax)<self.tol and self.PeriodicDir[2] == 1:
            y[0] = x[0]
            y[1] = x[1]
            y[2] = x[2] - (self.zmax-self.zmin) + Lzz 
        else:
            y[0] = x[0] - (self.xmax-self.xmin) + Lxx
            y[1] = x[1] - (self.ymax-self.ymin) + Lyy
            if (self.gdim==3):
              y[2] = x[2] - (self.zmax-self.zmin) + Lzz


def MassMatrix(mydomain):
    # Wrapper for the mass matrix
    w = mydomain.w; v = mydomain.v;
    if mydomain.IsDomainMultiple==True:
        return MassMatrix2c(w, v, mydomain.phase)
    else:
        return MassMatrix1c(w, v)

class WeakPseudoPeriodic_1c(UserExpression):
    def __init__(self, mydomain, **kwargs):
        self.xmin, self.ymin, self.zmin, self.xmax, self.ymax, self.zmax = mydomain.xmin, mydomain.ymin, mydomain.zmin, mydomain.xmax, mydomain.ymax, mydomain.zmax 
        self.pdir = mydomain.PeriodicDir
        self.gdim = mydomain.gdim
        self.gdir = mydomain.gdir
        self.qvalue = mydomain.qvalue       
        super().__init__(**kwargs)
        

    def set_values(self, u, ift):
        u.set_allow_extrapolation(True)
        self.u   = u
        self.ift = ift
        
    def eval(self, value, x):
        tol = 1E-7
        
        xln, yln, zln = 0, 0, 0
        temp_u0r, temp_u0i = 0, 0
        if (self.pdir[0] == 1): # x-direction
            if abs(x[0]-self.xmin) <= tol:
                xln, yln, zln = self.xmax-x[0], 0, 0;
                if self.gdim ==2:
                  temp_u0r, temp_u0i = self.u([self.xmax, x[1]]);
                if self.gdim ==3:
                  temp_u0r, temp_u0i = self.u([self.xmax, x[1], x[2]]);

            if abs(x[0]-self.xmax) <= tol:
                xln, yln, zln = self.xmin-x[0], 0, 0;
                if (self.gdim==2):
                  temp_u0r, temp_u0i = self.u([self.xmin, x[1]]);
                if (self.gdim==3):
                  temp_u0r, temp_u0i = self.u([self.xmin, x[1], x[2]]);

        if (self.pdir[1] == 1): # y-direction
            if abs(x[1]-self.ymin) <= tol:
                xln, yln, zln = 0, self.ymax-x[1], 0;
                if self.gdim ==2:
                  temp_u0r, temp_u0i = self.u([x[0],self.ymax]);
                if self.gdim ==3:
                  temp_u0r, temp_u0i = self.u([x[0], self.ymax, x[2]]);

            if abs(x[1]-self.ymax) <= tol:
                xln, yln, zln = 0, self.ymin-x[1], 0;
                if (self.gdim==2):
                  temp_u0r, temp_u0i = self.u([x[0], self.ymin]);
                if (self.gdim==3):
                  temp_u0r, temp_u0i = self.u([x[0], self.ymin, x[2]]);

        if (self.pdir[2] == 1 and self.gdim==3): # z-direction
            if abs(x[2]-self.zmin) <= tol:
                xln, yln, zln = 0, 0, self.zmax-x[2];
                temp_u0r, temp_u0i = self.u([x[0], x[1], self.zmax]);

            if abs(x[2]-self.zmax) <= tol:
                xln, yln, zln = 0, 0, self.zmin-x[2];
                temp_u0r, temp_u0i = self.u([x[0], x[1], self.zmin]);

        # exp(i*theta)= cos(theta)+i*sin(theta); 
        theta_ln = self.qvalue*(self.gdir.x()*xln + self.gdir.y()*yln + self.gdir.z()*zln)*self.ift;
        value[0] = temp_u0r*cos(theta_ln)-temp_u0i*sin(theta_ln); # for Real part
        value[1] = temp_u0r*sin(theta_ln)+temp_u0i*cos(theta_ln); # for Imag part
    def value_shape(self):
        return (2,)


def MyFunctionSpaces(mydomain, periodicBD):  
  comm = MPI.comm_world
  rank = comm.Get_rank()
  
  Ve = FiniteElement("CG", mydomain.mymesh.ufl_cell(), mydomain.porder)
      
  if (mydomain.IsDomainMultiple==True):
        TH = MixedElement([Ve,Ve,Ve,Ve])
        if rank == 0:
            print("Function Space for Two-compartment Domains has 4 components");
            print("(ur0, ui0, ur1, ur1): r-real, i-imaginary")
  else:
        TH = MixedElement([Ve,Ve])   
        if rank == 0:
            print("Function Space for Single Domains has 2 components");
            print("(ur, ui): r-real, i-imaginary")
  if mydomain.IsDomainPeriodic==False or periodicBD==None:
        if rank == 0:
            print("Initialize a standard function space.")
        if sum(mydomain.PeriodicDir)>0 and rank==0:
            print("The pseudo-periodic BCS are weakly imposed.")
            print("The mesh does not need to be periodic.")
        V_DG = FunctionSpace(mydomain.mymesh, 'DG', 0)    
        V = FunctionSpace(mydomain.mymesh,Ve);
        W = FunctionSpace(mydomain.mymesh, TH)
  else:
        if rank == 0:
            print("Initialize peridodic function spaces.")
            print("The pseudo-periodic BCS are strongly imposed.")
            print("The mesh needs to be periodic.")
        V_DG = FunctionSpace(mydomain.mymesh, 'DG', 0, constrained_domain=periodicBD)
        V = FunctionSpace(mydomain.mymesh,Ve, constrained_domain=periodicBD)
        W = FunctionSpace(mydomain.mymesh, TH, constrained_domain=periodicBD)    
  return Ve, V, W, V_DG

def CheckAndCorrectPeriodicity(mesh, direction, tol):
    print("Check and correct the periodicity in direction "+str(direction))
    gdim = mesh.geometry().dim()
    if (direction>=gdim):
        print("  Direction "+str(direction)+" is invalid");
        return;
    bmesh  = BoundaryMesh(mesh, "exterior")   # surface boundary mesh
    mesh_points=mesh.coordinates()
    bmesh_points=bmesh.coordinates()
    
    vertmap = bmesh.entity_map(0)

    xmin = bmesh_points[:, 0].min()
    xmax = bmesh_points[:, 0].max()

    ymin = bmesh_points[:, 1].min()
    ymax = bmesh_points[:, 1].max()

    zmin, zmax = 0, 0;
    if (gdim==3):
        zmin = bmesh_points[:, 2].min()
        zmax = bmesh_points[:, 2].max()

    if direction==0:
        cmin, cmax = xmin, xmax
    if direction==1:
        cmin, cmax = ymin, ymax
    if direction==2:
        cmin, cmax = zmin, zmax

    master=[];slave=[];

    for v in vertices(bmesh):
        global_vindex = vertmap[v.index()]
        x = bmesh_points[v.index()];
        x2 = 0;
        if (gdim==3):
            x2 = x[2]
            
        if abs(x[direction]-cmin)<tol:
            master.append((x[0],x[1],x2, global_vindex))
        if abs(x[direction]-cmax)<tol:
            slave.append((x[0],x[1],x2, global_vindex))

    if not(len(master)==len(slave)):
          print("  The mesh is not periodic! The number of vertices on the opposite sides is not the same.");
          return
        
    sorter = lambda x: (x[0], x[1], x[2])

    sorted_master = sorted(master)
    sorted_slave = sorted(slave)

    x_str=["x-", "y-", "z-"];
    
    import numpy
    for i in range(0,3):
        if not(i==direction):
            masterX = list(zip(*sorted_master))[i]
            slaveX = list(zip(*sorted_slave))[i]
            error = abs(numpy.subtract(masterX, slaveX)).max()
            print("  The maximum error in "+x_str[i]+"component is: "+str(error))
            if (error>tol):
                print("  The mesh is not periodic with the given tol="+str(tol))
                return
              
    masterID = list(zip(*sorted_master))[3]
    slaveID = list(zip(*sorted_slave))[3] 

    for i in range(0,len(masterID)):
        for j in range(0,gdim):
          if not(j==direction):
              mesh_points[masterID[i]][j]=mesh_points[slaveID[i]][j]
          else:
              mesh_points[masterID[i]][j]=-mesh_points[slaveID[i]][j]
    print("  The mesh was successfully corrected for this direction.")

def GetGlobalDomainSize(mesh, mpi4py, numpy):
    gdim = mesh.geometry().dim()
    comm = mesh.mpi_comm()
    mcoors = mesh.coordinates()
    lxmin, lymin = mcoors[:,0].min(), mcoors[:,1].min()
    lxmax, lymax = mcoors[:,0].max(), mcoors[:,1].max()
    lzmin, lzmax = 0, 0
    if gdim==3:
        lzmin, lzmax = mcoors[:,2].min(), mcoors[:,2].max()

    Dsize = [lxmin, lymin, lzmin, lxmax, lymax, lzmax]
    local_size = numpy.array(Dsize, 'd')
    global_size = numpy.zeros(len(Dsize), dtype='d')

    comm.Allreduce(local_size[0:3], global_size[0:3], mpi4py.MPI.MIN)
    comm.Allreduce(local_size[3:6], global_size[3:6], mpi4py.MPI.MAX)
    return global_size
  
class MyDomain():  
    def __init__(self, mymesh, mri_para):
        comm = mymesh.mpi_comm()
        rank = comm.Get_rank()
        self.porder = 1                                  # order of basis functions of FEM
        self.hmin = MPI.min(comm, mymesh.hmin())
        self.hmax = MPI.max(comm, mymesh.hmax())
        self.tol = 1e-2*self.hmin
        self.gdim = mymesh.geometry().dim()
        self.tdim = mymesh.topology().dim()
        self.xmin, self.ymin, self.zmin, self.xmax, self.ymax, self.zmax=GetGlobalDomainSize(mymesh, mpi4py, numpy)
        if rank == 0:
            print("Domain size: xmin=%f, ymin=%f, zmin=%f, xmax=%f, ymax=%f, zmax=%f"%(self.xmin, self.ymin, self.zmin, self.xmax, self.ymax, self.zmax))
        self.mymesh = mymesh;            
        self.gdir = mri_para.gdir        
        self.qvalue = mri_para.qvalue 
        self.kappa_e_scalar = 3e-3/self.hmin

    def WeakPseudoPeridicMarker(self):        
        if self.gdim==2:
            pmk = self.kappa_e_scalar*Expression("(x[0]<xmin+eps || x[0]>xmax-eps)*p0 || (x[1]<ymin+eps || x[1]>ymax-eps)*p1", 
                             xmin=self.xmin, xmax=self.xmax, ymin=self.ymin, ymax=self.ymax, 
                             eps=self.tol, p0 = self.PeriodicDir[0], p1 = self.PeriodicDir[1], domain=self.mymesh, degree=1);
        if self.gdim==3:
            pmk = self.kappa_e_scalar*Expression("(x[0]<xmin+eps || x[0]>xmax-eps)*p0 || (x[1]<ymin+eps || x[1]>ymax-eps)*p1 || (x[2]<zmin+eps || x[2]>zmax-eps)*p2", 
                             xmin=self.xmin, xmax=self.xmax, ymin=self.ymin, ymax=self.ymax, zmin=self.zmin, zmax=self.zmax, 
                             eps=self.tol, p0 = self.PeriodicDir[0], p1 = self.PeriodicDir[1], p2 = self.PeriodicDir[2], domain=self.mymesh, degree=1);
        return pmk
    def ImposeDiffusionTensor(self, k00, k01, k02, k10, k11, k12, k20, k21, k22):
        print("Impose Diffusion Tensor ...")
        if self.gdim==2:
            self.D = as_matrix(((k00, k01), (k10, k11)))
        if self.gdim==3:
            self.D = as_matrix(((k00, k01, k02), (k10, k11, k12), (k20, k21, k22)))
        
    def Apply(self):      
        self.fn = FacetNormal(self.mymesh);
        
        if self.IsDomainMultiple == True: 
              self.fn0 = ieval(self.fn, 0, self.phase);
        
        self.kappa_e = self.WeakPseudoPeridicMarker()
        
        if (sum(self.PeriodicDir)>0):
                periodicBD = PeriodicBD(self) 
        else:
                periodicBD = None
      
        self.Ve, self.V, self.W, self.V_DG = MyFunctionSpaces(self, periodicBD)
        self.v = TestFunction(self.W); self.w = TrialFunction(self.W);
      
        if self.IsDomainPeriodic == True or sum(self.PeriodicDir)==0:
                self.wpperiodic = None
        else:
                self.wpperiodic = ComputeWeakBC(self)



def ThetaMethodF(ft, ift, mri_para, mri_simu, mydomain):
    w = mydomain.w
    v = mydomain.v
    if (mydomain.IsDomainMultiple==True):
      if (mydomain.IsDomainPeriodic==True) and sum(mydomain.PeriodicDir)>0:
          F = ThetaMethodF_sBC2c(ft, ift, mri_para, w, v, mri_simu, mydomain)
      else:
          F = ThetaMethodF_wBC2c(ft, ift, mri_para, w, v, mri_simu, mydomain)
    else:
      if (mydomain.IsDomainPeriodic==True) and sum(mydomain.PeriodicDir)>0:
          F = ThetaMethodF_sBC1c(ft, ift, mri_para, w, v, mri_simu, mydomain)
      else:
          F = ThetaMethodF_wBC1c(ft, ift, mri_para, w, v, mri_simu, mydomain)
    return F      

def ThetaMethodL(ft, ift, mri_para, mri_simu, mydomain):
    w = mydomain.w
    v = mydomain.v  
    u_0 = mri_simu.u_0
    wpperiodic = mydomain.wpperiodic
    if (mydomain.IsDomainMultiple==True):
      if (mydomain.IsDomainPeriodic==True) and sum(mydomain.PeriodicDir)>0:
          L = ThetaMethodL_sBC2c(ft, ift, mri_para, w, v ,u_0, mri_simu, mydomain);
      else:
          L = ThetaMethodL_wBC2c(ft, ift, mri_para, w, v ,u_0, mri_simu, mydomain, wpperiodic);
    else:
      if (mydomain.IsDomainPeriodic==True) and sum(mydomain.PeriodicDir)>0:
          L = ThetaMethodL_sBC1c(ft, ift, mri_para, w, v ,u_0, mri_simu, mydomain);
      else:
          L = ThetaMethodL_wBC1c(ft, ift, mri_para, w, v ,u_0, mri_simu, mydomain, wpperiodic);
    return L    

def ComputeWeakBC(mydomain):
      if mydomain.IsDomainMultiple == True:
          wpperiodic = WeakPseudoPeriodic_2c(mydomain, degree=1)
      else:
          wpperiodic = WeakPseudoPeriodic_1c(mydomain, degree=1)
      return wpperiodic
    
def convert_g2q(gvalue):
    g_ratio = 2.675e8;
    return gvalue*g_ratio*1e-12
  
def convert_q2g(qvalue):
    g_ratio = 2.675e8;
    return qvalue/g_ratio*1e12

def GetPartitionMarkers(mshfile, pfile=""):

  fext=mshfile[len(mshfile)-3:len(mshfile)]
  if not(fext =="msh"):
    print("Sorry!!! The program supports .msh only! Please check the input file format...");
    exit(0);

  print("Extracting cell markers from: "+mshfile+" ...");

  filewithoutextention = os.path.splitext(os.path.basename(mshfile))[0];

  f = open(mshfile, "r")
  lineList = f.readlines()
  f.close()

  lastline_index = lineList.index('$EndElements\n');

  lastline = lineList[lastline_index-1].split(" ")
  last_id = int(lastline[0])
    
  first_id = -100;

  for x in lineList:
    x = x.split(" ")
    if len(x) == len(lastline):
        first_id = int(x[0])
        break;

  if len(lastline)==8:
        dim = 2;
  elif len(lastline)==9:
        dim = 3;
  else:
        print("Invalid dimension. Please double check msh2xml!");

  numcells = last_id - first_id + 1

  if (pfile==""):
     outfile = 'pmk_'+filewithoutextention+'.xml'
  else:
     outfile = pfile

  filename = open(outfile,'w');
  filename.write('<?xml version="1.0"?>\n')
  filename.write('<dolfin xmlns:dolfin="http://fenicsproject.org">\n')
  filename.write('  <mesh_function>\n')
  filename.write('    <mesh_value_collection type="uint" dim="'+str(dim)+'" size="'+str(numcells)+'">\n');

  cmpts=[];
  for x in lineList:
    x = x.split(" ")
    if len(x) == len(lastline):
          filename.write('      <value cell_index="'+str(int(x[0])-first_id)+'" local_entity="0" value="'+x[3]+'" />\n');
          if (len(cmpts)==0 or not(x[3] in cmpts)):
              cmpts.append(x[3])

  filename.write('    </mesh_value_collection>\n');
  filename.write('  </mesh_function>\n');
  filename.write('</dolfin>\n');
  filename.close();

  print("Extracted successfully on: "+str(numcells)+" elements")
  print("Partition marker list: " + str(cmpts))
  print("Wrote to: "+outfile)
  return cmpts  
  
def CreatePhaseFunc(mymesh, evengroup, oddgroup, partition_marker):
    V_DG = FunctionSpace(mymesh, 'DG', 0)
    dofmap_DG = V_DG.dofmap()
    phase = Function(V_DG)
    if not(partition_marker==None): # partition_marker is given
        partion_list = [];
        for cell in cells(mymesh):
            cmk = partition_marker[cell.index()]
            if (len(partion_list)==0 or not(cmk in partion_list)):
                partion_list.append(cmk)
            phase.vector()[dofmap_DG.cell_dofs(cell.index())] = cmk % 2  
        return phase, partion_list
    else: # partition_marker is not given
        partition_marker = MeshFunction("size_t", mymesh, mymesh.topology().dim())
        partion_list = [-1]*(len(evengroup)+len(oddgroup))
        if len(evengroup)==0:
            phase.vector()[:] = 0
        elif len(oddgroup)==0:
            phase.vector()[:] = 1
        elif (len(evengroup)+len(oddgroup)==0):
            print("At least one of evengroup, oddgroup, partition_marker is not empty!")
            
        for cell in cells(mymesh):
            p = cell.midpoint();
            for submesh_id in range(0, len(evengroup)):
                is_inside = evengroup[submesh_id].bounding_box_tree().compute_first_entity_collision(p)<4294967295
                if is_inside==True:
                    cmk = 2*submesh_id + 2;
                    partion_list[submesh_id] = cmk;
                    partition_marker[cell.index()] = cmk;
                    phase.vector()[dofmap_DG.cell_dofs(cell.index())] = 0;
                    break;
            for submesh_id in range(0, len(oddgroup)):
                is_inside = oddgroup[submesh_id].bounding_box_tree().compute_first_entity_collision(p)<4294967295
                if is_inside==True:
                    cmk = 2*submesh_id + 1
                    partion_list[submesh_id+len(evengroup)] = cmk;
                    partition_marker[cell.index()] = cmk;
                    phase.vector()[dofmap_DG.cell_dofs(cell.index())] = 1;
                    break;
        for cell in cells(mymesh):
            cmk = partition_marker[cell.index()]
            if (len(partion_list)==0 or not(cmk in partion_list)):
                partion_list.insert(0,cmk)
        return phase, partion_list, partition_marker
      
class MRI_parameters():
    def __init__(self):
        # Initialize default parameters
        self.bvalue = None
        self.gvalue = None        
        self.gdir = [1, 0, 0];
        self.nperiod = 0; # number of period for OGSE sequences
        self.T2 = 1e16    # T2 relaxation time
        self.s = sp.Symbol('s')
    def set_gradient_dir(self, mymesh, g0, g1, g2):
        gdim = mymesh.geometry().dim()
        if gdim==2:
            self.gdir = Point(g0, g1)
            if abs(self.gdir.norm())>1e-10:
                  self.gdir /= self.gdir.norm()
            else:
                print("|g|=0! Please check again the gradient directions!"); sys.exit()
            self.g = Expression(("g0","g1"), g0=self.gdir.x(), g1=self.gdir.y(),domain=mymesh, degree=1);
        if gdim==3:
            self.gdir = Point(g0, g1, g2)
            if abs(self.gdir.norm())>1e-10:
                self.gdir /= self.gdir.norm()
            else:
                print("|g|=0! Please check again the gradient directions!"); sys.exit() 
            self.g = Expression(("g0","g1","g2"), g0=self.gdir.x(), g1=self.gdir.y(), g2=self.gdir.z(),domain=mymesh, degree=1);
            
    def integral_term_for_gb(self):
        self.int4gb = float(sp.integrate(self.ifs_sym*self.ifs_sym, (self.s, 0, self.T)))

    def itime_profile_sym(self):
        # Return symbolic int_0^s f(u) du, change variable to ift(s)
        u = sp.Symbol('u')
        self.ifs_sym = sp.integrate(self.fs_sym.subs(self.s, u), (u, 0, self.s))
        
    def time_profile(self, t):
        return (float(self.fs_sym.subs(self.s,t)))
      
    def itime_profile(self, t): 
        return (float(self.ifs_sym.subs(self.s, t)))

    def convert_b2q(self):
        self.qvalue = sqrt(self.bvalue)/sqrt(self.int4gb);
        return self.qvalue
    def convert_q2b(self):
        self.bvalue = self.qvalue*self.qvalue*self.int4gb;
        return self.bvalue
      
    def Apply(self):
        self.itime_profile_sym(); 
        self.integral_term_for_gb();
        if not(self.bvalue==None):
            self.qvalue = self.convert_b2q();
            self.gvalue = convert_q2g(self.qvalue);
        elif not(self.gvalue==None):
            self.qvalue = convert_g2q(self.gvalue);
            self.bvalue = self.convert_q2b();
        elif (self.bvalue==None and self.bvalue==None):
            print("bvalue or gvalue need to be specified.")
            sys.exit()      
          
class MRI_simulation():
    def __init__(self):
          self.nskip = 5;    # Output frequency (for visualization only)
          self.theta = 0.5;  # theta=0.5: midpoint method

    def InitialCondition(self, mydomain, Dirac_Delta):
          if Dirac_Delta==None:
              if mydomain.gdim==2:
                  Dirac_Delta = Expression("x[0]*x[0]+x[1]*x[1]<eps",eps=1e6, domain=mydomain.mymesh, degree=1);
              if mydomain.gdim==3:
                  Dirac_Delta = Expression("x[0]*x[0]+x[1]*x[1]+x[2]*x[2]<eps",eps=1e6, domain=mydomain.mymesh, degree=1);
              Dirac_Delta = interpolate(Dirac_Delta, mydomain.V);
          u_0 = Function(mydomain.W);
          assign(u_0.sub(0), Dirac_Delta)
          if (mydomain.IsDomainMultiple==True):
              assign(u_0.sub(2), Dirac_Delta)  
          return Dirac_Delta, u_0
        
    def solve(self, mydomain, mri_para, linsolver, ic=None): 

          self.Dirac_Delta, self.u_0 = self.InitialCondition(mydomain, ic)

          stepcounter = 0;

          M = MassMatrix(mydomain);

          tp = 0;
          self.t = tp;

          comm = MPI.comm_world
          rank = comm.Get_rank()

          ft_f, ift_f, ft_p_f, ift_p_f = Function(mydomain.V), Function(mydomain.V), Function(mydomain.V), Function(mydomain.V);
          F = ThetaMethodF(ft_f, ift_f, mri_para, self, mydomain)
          L = ThetaMethodL(ft_p_f, ift_p_f, mri_para, self, mydomain)

          start_time = time.time()
          while self.t < mri_para.T + self.k: # Time-stepping loop                                                                                                                                
              if stepcounter % self.nskip == 0 and rank==0:
                  print('t: %6.2f '%self.t, 'T: %6.2f'%mri_para.T, 'dt: %.1f'%self.k,'qvalue: %e'%mri_para.qvalue,'Completed %3.2f%%'%(float(self.t)/float(mri_para.T+self.k)*100.0));

              ft_f.vector()[:]   = mri_para.time_profile(self.t);   ift_f.vector()[:]   = mri_para.itime_profile(self.t);
              ft_p_f.vector()[:] = mri_para.time_profile(tp);       ift_p_f.vector()[:] = mri_para.itime_profile(tp);
                            
              A = 1/self.k*M + assemble(F);
              b = assemble(L);

              linsolver.solve(A, self.u_0.vector(),b);

              tp = self.t;
              self.t += self.k;
              stepcounter += 1;
 
          self.elapsed_time = time.time() - start_time
          if rank==0:
              print("Successfully Completed! Elapsed time: %f seconds"%self.elapsed_time)

def PostProcessing(mydomain, mri_para, mri_simu, plt, ms=''):
    comm = MPI.comm_world
    rank = comm.Get_rank()    

    one = Function(mydomain.V)
    one.vector()[:] = 1
    whole_vol = assemble(one*dx)
    voi = assemble(mri_simu.Dirac_Delta*dx)
    if mydomain.IsDomainMultiple == True:
        u0r_0, u0i_0, u1r_0, u1i_0 = split(mri_simu.u_0)
        initial0 = assemble((1-mydomain.phase)*mri_simu.Dirac_Delta*dx);
        signal0 = assemble(((1-mydomain.phase)*u0r_0)*dx);
        initial1 = assemble(mydomain.phase*mri_simu.Dirac_Delta*dx)
        signal1 = assemble((mydomain.phase*u1r_0)*dx);
        signal = assemble((mydomain.phase*u1r_0+(1-mydomain.phase)*u0r_0)*dx);
        if numpy.isscalar(mydomain.kappa)==True:
          out_text = 'b: %.3f, g: %.3f, q: %.3e, Signal: %.3e, Normalized signal: %.6e, kappa: %.3e, dt: %.3f, hmin: %.3e, hmax: %.3e, whole_vol: %.3f, vol_of_interest: %.3f, elasped time %.3f (s)\n'%(mri_para.bvalue, mri_para.gvalue, mri_para.qvalue, signal, signal/voi, mydomain.kappa, mri_simu.k, mydomain.hmin, mydomain.hmax, whole_vol, voi, mri_simu.elapsed_time)
        else:
          out_text = 'b: %.3f, g: %.3f, q: %.3e, Signal: %.3e, Normalized signal: %.6e, dt: %.3f, hmin: %.3e, hmax: %.3e, whole_vol: %.3f, vol_of_interest: %.3f, elasped time %.3f (s)\n'%(mri_para.bvalue, mri_para.gvalue, mri_para.qvalue, signal, signal/voi, mri_simu.k, mydomain.hmin, mydomain.hmax, whole_vol, voi, mri_simu.elapsed_time)            
        if rank==0:
            print('Signal on each compartment')
            print('Sum initial0: %.3e, Signal0: %.3e'%(initial0, signal0))
            print('Sum initial1: %.3e, Signal1: %.3e'%(initial1, signal1))
            print(out_text)
        try:               
            dm = mydomain.V_DG.dofmap()
            mydomain.mphase = MeshFunction("size_t", mydomain.V_DG.mesh(), mydomain.V_DG.mesh().topology().dim())
            for cell in cells(mydomain.V_DG.mesh()):
                  mydomain.mphase[cell] = mydomain.phase.vector()[dm.cell_dofs(cell.index())]
                
            mydomain.mesh0 = SubMesh(mydomain.mymesh, mydomain.mphase, 0)
            mydomain.mesh1 = SubMesh(mydomain.mymesh, mydomain.mphase, 1)
            
            V0 = FunctionSpace(mydomain.mesh0, mydomain.Ve);
            V1 = FunctionSpace(mydomain.mesh1, mydomain.Ve);
            u0r_0p = project(u0r_0,V0)
            u1r_0p = project(u1r_0,V1)
            u0i_0p = project(u0i_0,V0)
            u1i_0p = project(u1i_0,V1)
            
            if mydomain.tdim==mydomain.gdim and not(plt==None):
                plt.figure(10000);
                plot(u0r_0p, cmap="coolwarm")
                plt.figure(10001);            
                plot(u1r_0p, cmap="coolwarm")  
            File("u0r.pvd")<<u0r_0p
            File("u1r.pvd")<<u1r_0p
            File("u0i.pvd")<<u0i_0p
            File("u1i.pvd")<<u1i_0p
        except:
            if rank==0:
                print("Could not post-process the solutions for the visualization purposes due to some reasons.")
    else:
        ur, ui = split(mri_simu.u_0)
        signal = assemble(ur*dx);
        out_text = 'b: %.3f, g: %.3f, q: %.3e, Signal: %.3e, Normalized signal: %.6e, dt: %.3f, hmin: %.3e, hmax: %.3e, whole_vol: %.3f, vol_of_interest: %.3f, elasped time %.3f (s)\n'%(mri_para.bvalue, mri_para.gvalue, mri_para.qvalue, signal, signal/voi, mri_simu.k, mydomain.hmin, mydomain.hmax, whole_vol, voi, mri_simu.elapsed_time)
        if rank==0:
            print(out_text)
        V = FunctionSpace(mydomain.mymesh,mydomain.Ve);
        ur_p = project(ur,V)
        if mydomain.tdim==mydomain.gdim and not(plt==None): 
            plt.figure(10000);
            plot(ur_p, cmap="coolwarm")
        File("ur.pvd")<<ur_p
        
    if int(rank) == 0:
        print("save to log.txt")
        outfile = open('log.txt', 'a')
        if not(ms == ''):
            outfile.write('%'+ms+'\n')
        outfile.write(out_text)
        outfile.close()
