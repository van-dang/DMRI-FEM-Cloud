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

import time, os, sys, shutil

def GdotX(gdir, mymesh):
  gdim = mymesh.geometry().dim()
  if (gdim==2):
    GX=Expression("x[0]*g0+x[1]*g1", g0=gdir.x(), g1=gdir.y(), domain=mymesh, degree=3);
  if (gdim==3):
    GX=Expression("x[0]*g0+x[1]*g1+x[2]*g2", g0=gdir.x(), g1=gdir.y(), g2=gdir.z(), domain=mymesh, degree=3);
  return GX;

def FuncF_wBC(ft, gnorm, gdir, ur, ui, vr, vi, D, mymesh, T2):
    GX=GdotX(gdir, mymesh)
    Fr =   ft*gnorm*GX*ui*vr - 1./T2*ur*vr - inner(D*grad(ur), grad(vr))
    Fi = - ft*gnorm*GX*ur*vi - 1./T2*ui*vi - inner(D*grad(ui), grad(vi))
    return Fr + Fi
  
def icondition_wBC(kappa, u0rm, u1rm, v0r, v1r, u0im, u1im, v0i, v1i):
    F_bcr = kappa*(u0rm-u1rm)*(v0r-v1r)
    F_bci = kappa*(u0im-u1im)*(v0i-v1i)
    return F_bcr + F_bci

def SubMeshSave(ur, ui, file_ur, file_ui, mymesh, n, stepcounter, dolfin_version):
  if dolfin_version=='1.6.0':
    V = FunctionSpace(mymesh, "CG", porder)
  else:
    # For FEniCS 2016, 2017
    Ve = FiniteElement("CG", mymesh.ufl_cell(), porder)
    V = FunctionSpace(mymesh, Ve)
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
  
def ThetaMethodL_wBC1c(ft, ift, mri_para , w , v, u_0, sp, mydomain, wpperiodic):
    D = mydomain.D
    T2 = mri_para.T2    
    mymesh = mydomain.mymesh
    k = sp.k
    theta = sp.theta
    gnorm = mri_para.gnorm
    gdir = mri_para.gdir
    u0r_0, u0i_0 = split(u_0)
    v0r, v0i = v[0], v[1]
    u0r, u0i = w[0], w[1]
    L0 = (u0r_0/k*v0r + u0i_0/k*v0i+(1-theta)*FuncF_wBC(ft, gnorm, gdir, u0r_0, u0i_0, v0r, v0i, D, mymesh, T2))*dx
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
    gnorm = mri_para.gnorm
    gdir = mri_para.gdir
    T2 = mri_para.T2
    kappa = mydomain.kappa
    v0r, v0i = v[0], v[1]
    u0r, u0i = w[0], w[1]
    a0 = (  -theta*FuncF_wBC(ft, gnorm, gdir, u0r, u0i, v0r, v0i, D, mymesh, T2))*dx
    a_pbc = 0;
    if sum(mydomain.PeriodicDir)>0:
      a_pbc = theta*mydomain.kappa_e*(u0r*v0r   + u0i*v0i)*ds;
    return a0+a_pbc
  
  
def ThetaMethodL_wBC2c(ft, ift, mri_para , w , v, u_0, sp, mydomain, wpperiodic):
    D = mydomain.D
    mymesh = mydomain.mymesh
    k = sp.k
    theta = sp.theta
    gnorm = mri_para.gnorm
    gdir = mri_para.gdir
    T2 = mri_para.T2
    kappa = mydomain.kappa
    phase = mydomain.phase
    u0r_0, u0i_0, u1r_0, u1i_0 = split(u_0)
    v0r, v0i, v1r, v1i = v[0], v[1], v[2], v[3]
    u0r, u0i, u1r, u1i = w[0], w[1], w[2], w[3]

    L0 = (u0r_0/k*v0r + u0i_0/k*v0i+(1-theta)*FuncF_wBC(ft, gnorm, gdir, u0r_0, u0i_0, v0r, v0i, D, mymesh, T2))*(1-phase)*dx
    L1 = (u1r_0/k*v1r +u1i_0/k*v1i+(1-theta)*FuncF_wBC(ft, gnorm, gdir, u1r_0, u1i_0, v1r, v1i, D, mymesh, T2))*phase*dx
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
    gnorm = mri_para.gnorm
    gdir = mri_para.gdir
    kappa = mydomain.kappa
    phase = mydomain.phase
    T2 = mri_para.T2
    v0r, v0i, v1r, v1i = v[0], v[1], v[2], v[3]
    u0r, u0i, u1r, u1i = w[0], w[1], w[2], w[3]
    a0 = (  -theta*FuncF_wBC(ft, gnorm, gdir, u0r  , u0i  , v0r, v0i, D, mymesh, T2))*(1-phase)*dx
    a1 = (  -theta*FuncF_wBC(ft, gnorm, gdir, u1r  , u1i  , v1r, v1i, D, mymesh, T2))*phase*dx
    a_bc = avg(  (theta*icondition_wBC(kappa, u0r  , u1r  , v0r, v1r, u0i  , u1i  , v0i, v1i)))*abs(jump(phase))*dS;
    
    a_pbc = 0;
    if sum(mydomain.PeriodicDir)>0:
      a_pbc = theta*mydomain.kappa_e*(u1r*v1r   + u1i*v1i)*phase*ds + theta*mydomain.kappa_e*(u0r*v0r   + u0i*v0i)*(1-phase)*ds;
      
    return a0+a1+a_bc+a_pbc


# Strong periodic boundary conditions  
def inner_interface(ift, kappa, gnorm, u0rm, u1rm, v0r, v1r, u0im, u1im, v0i, v1i, fn, g, D0, D1):
    F_bcr  = (-kappa*avg(u0rm-u1rm)-0.5*gnorm*ift*(avg(u0im)*inner(avg(D0)*avg(g),fn)+avg(u1im)*inner(avg(D1)*avg(g),fn)))*avg(v0r-v1r)                      
    F_bcr += -gnorm*ift*( avg(u0im)*inner(avg(D0)*avg(g),fn)-avg(u1im)*inner(avg(D1)*avg(g),fn) )*0.5*avg(v0r+v1r)                                                                                                                                                                                         
    F_bci  = (-kappa*avg(u0im-u1im)+0.5*gnorm*ift*(avg(u0rm)*inner(avg(D0)*avg(g),fn)+avg(u1rm)*inner(avg(D0)*avg(g),fn)))*avg(v0i-v1i)                      
    F_bci += gnorm*ift*(  avg(u0rm)*inner(avg(D0)*avg(g),fn)-avg(u1rm)*inner(avg(D1)*avg(g),fn) )*0.5*avg(v0i+v1i)                                         
    return -F_bcr - F_bci

def outer_interface(ift, gnorm, D, fn, ur, ui, vr, vi, g):
    F_bcr =  (ift*gnorm+1e-16)*inner(D*g, fn)*ui*vr
    F_bci = -(ift*gnorm+1e-16)*inner(D*g, fn)*ur*vi
    return F_bcr + F_bci

  
def FuncF_sBC(ift, gnorm, g, ur, ui, vr, vi, D, T2):
    Fr =   ift*gnorm*(inner(g,D*grad(ui))+inner(grad(ui),D*g))*vr - inner(g,D*g)*gnorm*gnorm*ift*ift*ur*vr-1./T2*ur*vr-inner(D*grad(ur), grad(vr))
    Fi = - ift*gnorm*(inner(g,D*grad(ur))+inner(grad(ur),D*g))*vi - inner(g,D*g)*gnorm*gnorm*ift*ift*ui*vi-1./T2*ui*vi-inner(D*grad(ui), grad(vi))
    return Fr + Fi


def ThetaMethodF_sBC1c(ft, ift, mri_para, w, v, sp, mydomain):
    D = mydomain.D  
    T2 = mri_para.T2  
    k = sp.k
    theta = sp.theta  
    gnorm = mri_para.gnorm
    g = mri_para.g
    fn = mydomain.fn
    kappa = mydomain.kappa
    v0r, v0i = v[0], v[1]
    u0r, u0i = w[0], w[1]
    a0 = (  -theta*FuncF_sBC(ift, gnorm, g, u0r  , u0i  , v0r, v0i, D, T2))*dx
    a0_outer_bc = theta*outer_interface(ift, gnorm , D, fn, u0r, u0i, v0r, v0i, g)*ds
    return a0 + a0_outer_bc 

  

def ThetaMethodL_sBC1c(ft, ift, mri_para, w, v, u_0, sp, mydomain):
    D = mydomain.D 
    T2 = mri_para.T2      
    k = sp.k
    theta = sp.theta  
    gnorm = mri_para.gnorm
    g = mri_para.g
    fn = mydomain.fn
    v0r, v0i = v[0], v[1]
    u0r, u0i = w[0], w[1]
    u0r_0, u0i_0 = split(u_0)
    L0 = (u0r_0/k*v0r + u0i_0/k*v0i +theta*FuncF_sBC(ift, gnorm, g, u0r_0, u0i_0, v0r, v0i, D, T2))*dx
    L0_outer_bc = -theta*outer_interface(ift, gnorm, D, fn, u0r_0, u0i_0, v0r, v0i, g)*ds
    return L0 + L0_outer_bc

  
def ThetaMethodF_sBC2c(ft, ift, mri_para, w, v, sp, mydomain):
    D = mydomain.D 
    T2 = mri_para.T2      
    k = sp.k
    theta = sp.theta  
    gnorm = mri_para.gnorm
    g = mri_para.g
    fn = mydomain.fn
    fn0 = mydomain.fn0
    kappa = mydomain.kappa
    phase = mydomain.phase
    v0r, v0i, v1r, v1i = v[0], v[1], v[2], v[3]
    u0r, u0i, u1r, u1i = w[0], w[1], w[2], w[3]
    a0 = (  -theta*FuncF_sBC(ift, gnorm, g, u0r  , u0i  , v0r, v0i, D, T2))*(1-phase)*dx
    a1 = (  -theta*FuncF_sBC(ift, gnorm, g, u1r  , u1i  , v1r, v1i, D, T2))*phase*dx
    a_inner_bc  = (  (theta*inner_interface(ift, kappa, gnorm, u0r  , u1r  , v0r, v1r, u0i  , u1i  , v0i, v1i, fn0, g, D, D)))*abs(jump(phase))*dS;
    a0_outer_bc = theta*outer_interface(ift, gnorm , D, fn, u0r, u0i, v0r, v0i, g)*ds
    a1_outer_bc = theta*outer_interface(ift, gnorm , D, fn, u1r, u1i, v1r, v1i, g)*ds
    return a0+a1 +a_inner_bc + a0_outer_bc + a1_outer_bc

  
def ThetaMethodL_sBC2c(ft, ift, mri_para, w, v, u_0, sp, mydomain):
    D = mydomain.D 
    T2 = mri_para.T2      
    k = sp.k
    theta = sp.theta  
    gnorm = mri_para.gnorm
    g = mri_para.g
    fn = mydomain.fn
    fn0 = mydomain.fn0
    kappa = mydomain.kappa
    phase = mydomain.phase
    v0r, v0i, v1r, v1i = v[0], v[1], v[2], v[3]
    u0r, u0i, u1r, u1i = w[0], w[1], w[2], w[3]
    u0r_0, u0i_0, u1r_0, u1i_0 = split(u_0)

    L0 = (u0r_0/k*v0r + u0i_0/k*v0i +theta*FuncF_sBC(ift, gnorm, g, u0r_0, u0i_0, v0r, v0i, D, T2))*(1-phase)*dx
    L1 = (u1r_0/k*v1r + u1i_0/k*v1i +theta*FuncF_sBC(ift, gnorm, g, u1r_0, u1i_0, v1r, v1i, D, T2))*phase*dx
    L_inner_bc  = -(theta*inner_interface(ift, kappa, gnorm, u0r_0, u1r_0, v0r, v1r, u0i_0, u1i_0, v0i, v1i, fn0, g, D, D))*abs(jump(phase))*dS;
    L0_outer_bc = -theta*outer_interface(ift, gnorm, D, fn, u0r_0, u0i_0, v0r, v0i, g)*ds
    L1_outer_bc = -theta*outer_interface(ift, gnorm, D, fn, u1r_0, u1i_0, v1r, v1i, g)*ds
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
        self.gnorm = mydomain.gnorm        
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
        theta_ln = self.gnorm*(self.gdir.x()*xln + self.gdir.y()*yln + self.gdir.z()*zln)*self.ift;
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
        self.gnorm = mydomain.gnorm       
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
        theta_ln = self.gnorm*(self.gdir.x()*xln + self.gdir.y()*yln + self.gdir.z()*zln)*self.ift;
        value[0] = temp_u0r*cos(theta_ln)-temp_u0i*sin(theta_ln); # for Real part
        value[1] = temp_u0r*sin(theta_ln)+temp_u0i*cos(theta_ln); # for Imag part
    def value_shape(self):
        return (2,)


def MyFunctionSpaces(mydomain, periodicBD):  
  Ve = FiniteElement("CG", mydomain.mymesh.ufl_cell(), mydomain.porder)
      
  if (mydomain.IsDomainMultiple==True):
        TH = MixedElement([Ve,Ve,Ve,Ve])
        print("Function Space for Two-compartment Domains has 4 components");
        print("(ur0, ui0, ur1, ur1): r-real, i-imaginary")
  else:
        TH = MixedElement([Ve,Ve])    
        print("Function Space for Single Domains has 2 components");
        print("(ur, ui): r-real, i-imaginary")
  if mydomain.IsDomainPeriodic==False or periodicBD==None:
        print("Initialize a standard function space.")
        if sum(mydomain.PeriodicDir)>0:
            print("The pseudo-periodic BCS are weakly imposed.")
            print("The mesh does not need to be periodic.")
        V_DG = FunctionSpace(mydomain.mymesh, 'DG', 0)    
        V = FunctionSpace(mydomain.mymesh,Ve);
        W = FunctionSpace(mydomain.mymesh, TH)
  else:
        print("Initialize peridodic function spaces.")
        print("The pseudo-periodic BCS are strongly imposed.")
        print("The mesh needs to be periodic.")
        V_DG = FunctionSpace(mydomain.mymesh, 'DG', 0, constrained_domain=periodicBD)
        V = FunctionSpace(mydomain.mymesh,Ve, constrained_domain=periodicBD)
        W = FunctionSpace(mydomain.mymesh, TH, constrained_domain=periodicBD)    
  return Ve, V, W, V_DG


class MyDomain():
    def __init__(self, mymesh, mri_para):
        self.mymesh = mymesh;
        self.porder = 1                                  # order of basis functions of FEM
        self.tol = 1e-6*mymesh.hmin()
        self.gdim = mymesh.geometry().dim()
        self.tdim = mymesh.topology().dim()
        self.hmin = mymesh.hmin()
        self.hmax = mymesh.hmax()      
        self.xmin = mymesh.coordinates()[:, 0].min()
        self.xmax = mymesh.coordinates()[:, 0].max()
        self.ymin = mymesh.coordinates()[:, 1].min()
        self.ymax = mymesh.coordinates()[:, 1].max()
        self.zmin, self.zmax = 0, 0 
        if (self.gdim==3):
            self.zmin = mymesh.coordinates()[:, 2].min()
            self.zmax = mymesh.coordinates()[:, 2].max()        

        self.gdir = mri_para.gdir        
        self.gnorm = mri_para.gnorm 
            
    def WeakPseudoPeridicMarker(self):        
        if self.gdim==2:
            pmk = 3e-3/self.hmin*Expression("(x[0]<xmin+eps || x[0]>xmax-eps)*p0 || (x[1]<ymin+eps || x[1]>ymax-eps)*p1", 
                             xmin=self.xmin, xmax=self.xmax, ymin=self.ymin, ymax=self.ymax, 
                             eps=1e-10, p0 = self.PeriodicDir[0], p1 = self.PeriodicDir[1], domain=self.mymesh, degree=1);
        if self.gdim==3:
            pmk = 3e-3/self.hmin*Expression("(x[0]<xmin+eps || x[0]>xmax-eps)*p0 || (x[1]<ymin+eps || x[1]>ymax-eps)*p1 || (x[2]<zmin+eps || x[2]>zmax-eps)*p2", 
                             xmin=self.xmin, xmax=self.xmax, ymin=self.ymin, ymax=self.ymax, zmin=self.zmin, zmax=self.zmax, 
                             eps=1e-10, p0 = self.PeriodicDir[0], p1 = self.PeriodicDir[1], p2 = self.PeriodicDir[2], domain=self.mymesh, degree=1);
        return pmk
    def ImposeDiffusionTensor(self, k00, k01, k02, k10, k11, k12, k20, k21, k22):
        print("Impose Diffusion Tensor")
        if self.gdim==2:
            self.D = as_matrix(((k00, k01), (k10, k11)))
        if self.gdim==3:
            self.D = as_matrix(((k00, k01, k02), (k10, k11, k12), (k20, k21, k22)))
        print(self.D)
        
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
    
def convert_q2g(qvalue):
    g_ratio = 2.675e8;
    return qvalue*g_ratio*1e-12
  
def convert_g2q(gnorm):
    g_ratio = 2.675e8;
    return gnorm/g_ratio*1e12

def Create_phase_func(mymesh, cmpt_mesh , pmk):
    V_DG = FunctionSpace(mymesh, 'DG', 0)
    dofmap_DG = V_DG.dofmap()
    phase = Function(V_DG)
    cellmarker = MeshFunction("size_t", mymesh, mymesh.topology().dim())
    partion_list = [];
    for cell in cells(mymesh):
        cmk = 0
        if not(pmk==None):
            cmk0 = pmk[cell.index()]
            cmk = cmk0 % 2    
        if not(cmpt_mesh==None):
            p = cell.midpoint();
            cmk0 = cmpt_mesh.bounding_box_tree().compute_first_entity_collision(p)<4294967295
            cmk = cmk0
        phase.vector()[dofmap_DG.cell_dofs(cell.index())] = cmk;
        cellmarker[cell.index()] = cmk; 
        
        if (len(partion_list)==0 or not(cmk0 in partion_list)):
            partion_list.append(cmk0)
            
    if pmk==None:
        pmk = cellmarker
    return cellmarker, phase, pmk, partion_list 

class MRI_parameters():
    def __init__(self):
        # Initialize default parameters
        self.bvalue = None
        self.qvalue = None        
        self.gdir = [1, 0, 0];
        self.nperiod = 0; # number of period for OGSE sequences
        self.T2 = 1e16    # T2 relaxation time
        self.s = sp.Symbol('s')
    def set_gradient_dir(self, mymesh, g0, g1, g2):
        gdim = mymesh.geometry().dim()
        if gdim==2:
            self.gdir = Point(g0, g1)
            self.gdir /= self.gdir.norm()
            self.g = Expression(("g0","g1"), g0=self.gdir.x(), g1=self.gdir.y(),domain=mymesh, degree=1);
        if gdim==3:
            self.gdir = Point(g0, g1, g2)
            self.gdir /= self.gdir.norm()        
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

    def convert_b2g(self):
        self.gnorm = 0;
        self.gnorm = sqrt(self.bvalue)/sqrt(self.int4gb);
        return self.gnorm
    def convert_g2b(self):
        self.bvalue = self.gnorm*self.gnorm*self.int4gb;
        return self.bvalue
      
    def Apply(self):
        self.itime_profile_sym(); 
        self.integral_term_for_gb();
        if not(self.bvalue==None):
            self.gnorm = self.convert_b2g();
            self.qvalue = convert_g2q(self.gnorm);
        elif not(self.qvalue==None):
            self.gnorm = convert_q2g(self.qvalue);
            self.bvalue = self.convert_g2b();
        elif (self.bvalue==None and self.bvalue==None):
            print("bvalue or qvalue need to be specified.")
            sys.exit()      
          
class MRI_simulation():
    def __init__(self):
          self.nskip = 5;    # Output frequency (for visualization only)
          self.theta = 0.5;  # theta=0.5: midpoint method

    def InitialCondition(self, mydomain):
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
        
    def solve(self, mydomain, mri_para, linsolver): 
      
          self.Dirac_Delta, self.u_0 = self.InitialCondition(mydomain)
          
          stepcounter = 0;

          M = MassMatrix(mydomain);

          self.t = 0;
          ft_prev  =  mri_para.time_profile(self.t);
          ift_prev = mri_para.itime_profile(self.t);

          start_time = time.time()
          
          while self.t < mri_para.T + self.k: # Time-stepping loop
              if stepcounter % self.nskip == 0:
                  print('t: %6.2f '%self.t, 'T: %6.2f'%mri_para.T, 'dt: %.1f'%self.k,'gnorm: %e'%mri_para.gnorm,'Completed %3.2f%%'%(float(self.t)/float(mri_para.T+self.k)*100.0));

              ft = mri_para.time_profile(self.t);
              ift = mri_para.itime_profile(self.t);
                            
              L = ThetaMethodL(ft_prev, ift_prev, mri_para, self, mydomain);
              A = 1/self.k*M + assemble(ThetaMethodF(ft, ift, mri_para, self, mydomain))

              b = assemble(L);
              
              linsolver.solve(A, self.u_0.vector(),b);

              ft_prev  = ft;
              ift_prev = ift;

              self.t += self.k;
              stepcounter += 1;

          self.elapsed_time = time.time() - start_time
          print("Successfully Completed! Elapsed time: %f seconds"%self.elapsed_time)
          
def Post_processing(mydomain, mri_para, mri_simu, ms=''):
    one = Function(mydomain.V)
    one.vector()[:] = 1
    whole_vol = assemble(one*dx)
    voi = assemble(mri_simu.Dirac_Delta*dx)
    if mydomain.IsDomainMultiple == True:
        u0r_0, u0i_0, u1r_0, u1i_0 = split(mri_simu.u_0)
        signal0 = assemble(((1-mydomain.phase)*u0r_0)*dx)/assemble((1-mydomain.phase)*mri_simu.Dirac_Delta*dx);
        signal1 = assemble((mydomain.phase*u1r_0)*dx)/assemble(mydomain.phase*mri_simu.Dirac_Delta*dx);
        signal = assemble((mydomain.phase*u1r_0+(1-mydomain.phase)*u0r_0)*dx);
        print('Signal on each compartment')
        print('Signal0: %.3e'%signal0)
        print('Signal1: %.3e'%signal1)
        out_text = 'b: %.3f, Signal: %.3e, Normalized signal: %.3e, kappa: %.3e, dt: %.3f, hmin: %.3f, whole_vol: %.3f, vol_of_interest: %.3f, Free signal: %.3e, elasped time %.3f (s)\n'%(mri_para.bvalue, signal, signal/voi, mydomain.kappa, mri_simu.k, mydomain.hmin, whole_vol, voi, exp(-mri_para.bvalue*D0), mri_sim.elapsed_time)
        
        V0 = FunctionSpace(mesh0, mydomain.Ve);
        V1 = FunctionSpace(mesh1, mydomain.Ve);
        u0r_0p = project(u0r_0,V0)
        u1r_0p = project(u1r_0,V1)
        if mydomain.tdim==mydomain.gdim:
            plot(u0r_0p, cmap="coolwarm")
            plot(u1r_0p, cmap="coolwarm")  
        File("u0r.pvd")<<u0r_0p
        File("u1r.pvd")<<u1r_0p
    else:
        ur, ui = split(mri_simu.u_0)
        signal = assemble(ur*dx);
        out_text = 'b: %.3f, Signal: %.3e, Normalized signal: %.3e, dt: %.3f, hmin: %.3e, whole_vol: %.3f, vol_of_interest: %.3f, Free signal: %.3e, elasped time %.3f (s)\n'%(mri_para.bvalue, signal, signal/voi, mri_simu.k, mydomain.hmin, whole_vol, voi, exp(-mri_para.bvalue*mydomain.D0), mri_sim.elapsed_time)
        V = FunctionSpace(mymesh,mydomain.Ve);
        ur_p = project(ur,V)
        if mydomain.tdim==mydomain.gdim:        
            plot(ur_p, cmap="coolwarm")
        File("ur.pvd")<<ur_p

    print(out_text)
    print("save to log.txt")
    outfile = open('log.txt', 'a')
    if not(ms == ''):
        outfile.write('%'+ms+'\n')
    outfile.write(out_text)
    outfile.close()
