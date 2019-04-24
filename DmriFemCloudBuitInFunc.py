from dolfin import *
def GdotX(gdir, mesh):
  gdim = mesh.geometry().dim()
  if (gdim==2):
    GX=Expression("x[0]*g0+x[1]*g1", g0=gdir.x(), g1=gdir.y(), domain=mesh, degree=3);
  if (gdim==3):
    GX=Expression("x[0]*g0+x[1]*g1+x[2]*g2", g0=gdir.x(), g1=gdir.y(), g2=gdir.z(), domain=mesh, degree=3);
  return GX;

def FuncF_wBC(ft, gnorm, gdir, ur, ui, vr, vi, D, mesh):
    GX=GdotX(gdir, mesh)
    Fr = ft*gnorm*GX*ui*vr - inner(D*grad(ur), grad(vr))
    Fi = - ft*gnorm*GX*ur*vi - inner(D*grad(ui), grad(vi))
    return Fr + Fi
  
def icondition_wBC(kappa, u0rm, u1rm, v0r, v1r, u0im, u1im, v0i, v1i):
    F_bcr = kappa*(u0rm-u1rm)*(v0r-v1r)
    F_bci = kappa*(u0im-u1im)*(v0i-v1i)
    return F_bcr + F_bci

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
  
def ThetaMethodL_wBC1c(ft, ift, mri_para , w , v, u_0, sp, mydomain, wpperiodic):
    D = mydomain.D
    mesh = mydomain.mesh
    k = sp.k
    theta = sp.theta
    gnorm = mri_para.gnorm
    gdir = mri_para.gdir
    u0r_0, u0i_0 = split(u_0)
    v0r, v0i = v[0], v[1]
    u0r, u0i = w[0], w[1]
    L0 = (u0r_0/k*v0r + u0i_0/k*v0i+(1-theta)*FuncF_wBC(ft, gnorm, gdir, u0r_0, u0i_0, v0r, v0i, D, mesh))*dx
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
    mesh = mydomain.mesh
    k = sp.k
    theta = sp.theta
    gnorm = mri_para.gnorm
    gdir = mri_para.gdir
    kappa = mydomain.kappa
    v0r, v0i = v[0], v[1]
    u0r, u0i = w[0], w[1]
    a0 = (  -theta*FuncF_wBC(ft, gnorm, gdir, u0r, u0i, v0r, v0i, D, mesh))*dx
    a_pbc = 0;
    if sum(mydomain.PeriodicDir)>0:
      a_pbc = theta*mydomain.kappa_e*(u0r*v0r   + u0i*v0i)*ds;
    return a0+a_pbc
  
  
def ThetaMethodL_wBC2c(ft, ift, mri_para , w , v, u_0, sp, mydomain, wpperiodic):
    D = mydomain.D
    mesh = mydomain.mesh
    k = sp.k
    theta = sp.theta
    gnorm = mri_para.gnorm
    gdir = mri_para.gdir
    kappa = mydomain.kappa
    phase = mydomain.phase
    u0r_0, u0i_0, u1r_0, u1i_0 = split(u_0)
    v0r, v0i, v1r, v1i = v[0], v[1], v[2], v[3]
    u0r, u0i, u1r, u1i = w[0], w[1], w[2], w[3]

    L0 = (u0r_0/k*v0r + u0i_0/k*v0i+(1-theta)*FuncF_wBC(ft, gnorm, gdir, u0r_0, u0i_0, v0r, v0i, D, mesh))*(1-phase)*dx
    L1 = (u1r_0/k*v1r +u1i_0/k*v1i+(1-theta)*FuncF_wBC(ft, gnorm, gdir, u1r_0, u1i_0, v1r, v1i, D, mesh))*phase*dx
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
    mesh = mydomain.mesh
    k = sp.k
    theta = sp.theta  
    gnorm = mri_para.gnorm
    gdir = mri_para.gdir
    kappa = mydomain.kappa
    phase = mydomain.phase

    v0r, v0i, v1r, v1i = v[0], v[1], v[2], v[3]
    u0r, u0i, u1r, u1i = w[0], w[1], w[2], w[3]
    a0 = (  -theta*FuncF_wBC(ft, gnorm, gdir, u0r  , u0i  , v0r, v0i, D, mesh))*(1-phase)*dx
    a1 = (  -theta*FuncF_wBC(ft, gnorm, gdir, u1r  , u1i  , v1r, v1i, D, mesh))*phase*dx
    a_bc = avg(  (theta*icondition_wBC(kappa, u0r  , u1r  , v0r, v1r, u0i  , u1i  , v0i, v1i)))*abs(jump(phase))*dS;
    
    a_pbc = 0;
    if sum(mydomain.PeriodicDir)>0:
      a_pbc = theta*mydomain.kappa_e*(u1r*v1r   + u1i*v1i)*phase*ds + theta*mydomain.kappa_e*(u0r*v0r   + u0i*v0i)*(1-phase)*ds;
      
    return a0+a1+a_bc+a_pbc


# Strong periodic boundary conditions  
def inner_interface(kappa, gnorm, u0rm, u1rm, v0r, v1r, u0im, u1im, v0i, v1i, fn, g, D0, D1):
    F_bcr  = (-kappa*avg(u0rm-u1rm)-0.5*gnorm*ift*(avg(u0im)*inner(D0*avg(g),fn)+avg(u1im)*inner(D1*avg(g),fn)))*avg(v0r-v1r)                      
    F_bcr += -gnorm*ift*( avg(u0im)*inner(D0*avg(g),fn)-avg(u1im)*inner(D1*avg(g),fn) )*0.5*avg(v0r+v1r)                                                                                                                                                                                         
    F_bci  = (-kappa*avg(u0im-u1im)+0.5*gnorm*ift*(avg(u0rm)*inner(D0*avg(g),fn)+avg(u1rm)*inner(D0*avg(g),fn)))*avg(v0i-v1i)                      
    F_bci += gnorm*ift*(  avg(u0rm)*inner(D0*avg(g),fn)-avg(u1rm)*inner(D1*avg(g),fn) )*0.5*avg(v0i+v1i)                                         
    return -F_bcr - F_bci

def outer_interface(ift, gnorm, D, fn, ur, ui, vr, vi, g):
    F_bcr =  (ift*gnorm+1e-16)*inner(D*g, fn)*ui*vr
    F_bci = -(ift*gnorm+1e-16)*inner(D*g, fn)*ur*vi
    return F_bcr + F_bci

  
def FuncF_sBC(ift, gnorm, g, ur, ui, vr, vi, D):
    Fr =   ift*gnorm*(inner(g,D*grad(ui))+inner(grad(ui),D*g))*vr - inner(g,D*g)*gnorm*gnorm*ift*ift*ur*vr-D*inner(grad(ur), grad(vr))
    Fi = - ift*gnorm*(inner(g,D*grad(ur))+inner(grad(ur),D*g))*vi - inner(g,D*g)*gnorm*gnorm*ift*ift*ui*vi-D*inner(grad(ui), grad(vi))
    return Fr + Fi


def ThetaMethodF_sBC1c(ft, ift, mri_para, w, v, sp, mydomain):
    D = mydomain.D  
    k = sp.k
    theta = sp.theta  
    gnorm = mri_para.gnorm
    g = mri_para.g
    fn = mydomain.fn
    fn0 = mydomain.fn0
    kappa = mydomain.kappa
    v0r, v0i = v[0], v[1]
    u0r, u0i = w[0], w[1]
    a0 = (  -theta*FuncF_sBC(ift, gnorm, g, u0r  , u0i  , v0r, v0i, D))*dx
    a0_outer_bc = theta*outer_interface(ift, gnorm , D, fn, u0r, u0i, v0r, v0i, g)*ds
    return a0 + a0_outer_bc 

  

def ThetaMethodL_sBC1c(ft, ift, mri_para, w, v, u_0, sp, mydomain):
    D = mydomain.D  
    k = sp.k
    theta = sp.theta  
    gnorm = mri_para.gnorm
    g = mri_para.g
    fn = mydomain.fn
    v0r, v0i = v[0], v[1]
    u0r, u0i = w[0], w[1]
    u0r_0, u0i_0 = split(u_0)
    L0 = (u0r_0/k*v0r + u0i_0/k*v0i +theta*FuncF_sBC(ift, gnorm, g, u0r_0, u0i_0, v0r, v0i, D))*dx
    L0_outer_bc = -theta*outer_interface(ift, gnorm, D, fn, u0r_0, u0i_0, v0r, v0i, g)*ds
    return L0 + L0_outer_bc

  
def ThetaMethodF_sBC2c(ft, ift, mri_para, w, v, sp, mydomain):
    D = mydomain.D  
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
    a0 = (  -theta*FuncF_sBC(ift, gnorm, g, u0r  , u0i  , v0r, v0i, D))*(1-phase)*dx
    a1 = (  -theta*FuncF_sBC(ift, gnorm, g, u1r  , u1i  , v1r, v1i, D))*phase*dx
    a_inner_bc  = (  (theta*inner_interface(kappa, gnorm, u0r  , u1r  , v0r, v1r, u0i  , u1i  , v0i, v1i, fn0, g, D, D)))*abs(jump(phase))*dS;
    a0_outer_bc = theta*outer_interface(ift, gnorm , D, fn, u0r, u0i, v0r, v0i, g)*ds
    a1_outer_bc = theta*outer_interface(ift, gnorm , D, fn, u1r, u1i, v1r, v1i, g)*ds
    return a0+a1 +a_inner_bc + a0_outer_bc + a1_outer_bc

  
def ThetaMethodL_sBC2c(ft, ift, mri_para, w, v, u_0, sp, mydomain):
    D = mydomain.D  
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

    L0 = (u0r_0/k*v0r + u0i_0/k*v0i +theta*FuncF_sBC(ift, gnorm, g, u0r_0, u0i_0, v0r, v0i, D))*(1-phase)*dx
    L1 = (u1r_0/k*v1r + u1i_0/k*v1i +theta*FuncF_sBC(ift, gnorm, g, u1r_0, u1i_0, v1r, v1i, D))*phase*dx
    L_inner_bc  = -(theta*inner_interface(kappa, gnorm, u0r_0, u1r_0, v0r, v1r, u0i_0, u1i_0, v0i, v1i, fn0, g, D, D))*abs(jump(phase))*dS;
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
    def __init__(self, mydomain, mri_para, **kwargs):
        self.xmin, self.ymin, self.zmin, self.xmax, self.ymax, self.zmax = mydomain.xmin, mydomain.ymin, mydomain.zmin, mydomain.xmax, mydomain.ymax, mydomain.zmax 
        self.pdir = mydomain.PeriodicDir
        self.gdim = mydomain.gdim
        self.gdir = mri_para.gdir        
        self.gnorm = mri_para.gnorm        
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
        self.tol = 1e-2*mydomain.mesh.hmin()
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


def MassMatrix(w, v, mydomain):
    # Wrapper for the mass matrix
    if mydomain.IsDomainMultiple==True:
        return MassMatrix2c(w, v, mydomain.phase)
    else:
        return MassMatrix1c(w, v)

class WeakPseudoPeriodic_1c(UserExpression):
    def __init__(self, mydomain, mri_para, **kwargs):
        self.xmin, self.ymin, self.zmin, self.xmax, self.ymax, self.zmax = mydomain.xmin, mydomain.ymin, mydomain.zmin, mydomain.xmax, mydomain.ymax, mydomain.zmax 
        self.pdir = mydomain.PeriodicDir
        self.gdim = mydomain.gdim
        self.gdir = mri_para.gdir
        self.gnorm = mri_para.gnorm       
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


def MyFunctionSpaces(mydomain, porder):  
  Ve = FiniteElement("CG", mydomain.mesh.ufl_cell(), porder)
  
  if (sum(mydomain.PeriodicDir)>0):
        periodicBD = PeriodicBD(mydomain) 
  else:
        periodicBD = None
      
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
        V_DG = FunctionSpace(mydomain.mesh, 'DG', 0)    
        V = FunctionSpace(mydomain.mesh,Ve);
        W = FunctionSpace(mydomain.mesh, TH)
  else:
        print("Initialize peridodic function spaces.")
        print("The pseudo-periodic BCS are strongly imposed.")
        print("The mesh needs to be periodic.")
        V_DG = FunctionSpace(mydomain.mesh, 'DG', 0, constrained_domain=periodicBD)
        V = FunctionSpace(mydomain.mesh,Ve, constrained_domain=periodicBD)
        W = FunctionSpace(mydomain.mesh, TH, constrained_domain=periodicBD)    
  return Ve, V, W, V_DG


class MyDomain():
    def __init__(self, mesh):
        self.mesh = mesh;
        self.tol = 1e-6*mesh.hmin()
        self.gdim = mesh.geometry().dim()
        self.hmin = mesh.hmin()
        self.hmax = mesh.hmax()      
        self.xmin = mesh.coordinates()[:, 0].min()
        self.xmax = mesh.coordinates()[:, 0].max()
        self.ymin = mesh.coordinates()[:, 1].min()
        self.ymax = mesh.coordinates()[:, 1].max()
        self.zmin, self.zmax = 0, 0 
        if (self.gdim==3):
            self.zmin = mesh.coordinates()[:, 2].min()
            self.zmax = mesh.coordinates()[:, 2].max()        
        
    def WeakPseudoPeridicMarker(self):        
        if self.gdim==2:
            pmk = 3e-3/self.hmin*Expression("(x[0]<xmin+eps || x[0]>xmax-eps)*p0 || (x[1]<ymin+eps || x[1]>ymax-eps)*p1", 
                             xmin=self.xmin, xmax=self.xmax, ymin=self.ymin, ymax=self.ymax, 
                             eps=1e-10, p0 = self.PeriodicDir[0], p1 = self.PeriodicDir[1], domain=self.mesh, degree=1);
        if self.gdim==3:
            pmk = 3e-3/self.hmin*Expression("(x[0]<xmin+eps || x[0]>xmax-eps)*p0 || (x[1]<ymin+eps || x[1]>ymax-eps)*p1 || (x[2]<zmin+eps || x[2]>zmax-eps)*p2", 
                             xmin=self.xmin, xmax=self.xmax, ymin=self.ymin, ymax=self.ymax, zmin=self.zmin, zmax=self.zmax, 
                             eps=1e-10, p0 = self.PeriodicDir[0], p1 = self.PeriodicDir[1], p2 = self.PeriodicDir[2], domain=self.mesh, degree=1);
        return pmk
    def ImposeDiffusionTensor(self, k00, k01, k02, k10, k11, k12, k20, k21, k22):
        print("Impose Diffusion Tensor")
        if self.gdim==2:
            self.D = as_matrix(((k00, k01), (k10, k11)))
        if self.gdim==3:
            self.D = as_matrix(((k00, k01, k02), (k10, k11, k12), (k20, k21, k22)))
        print(self.D)
        
    def Apply(self):      
        self.fn = FacetNormal(self.mesh);
        self.fn0 = ieval(self.fn, 0, self.phase);
        self.kappa_e = self.WeakPseudoPeridicMarker()


def ThetaMethodF(ft, ift, mri_para, mri_simu, mydomain):
    w = mri_simu.w
    v = mri_simu.v
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
    w = mri_simu.w
    v = mri_simu.v  
    u_0 = mri_simu.u_0
    wpperiodic = mri_simu.wpperiodic
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

def ComputeWeakBC(mydomain, mri_para):
      if mydomain.IsDomainMultiple == True:
          wpperiodic = WeakPseudoPeriodic_2c(mydomain, mri_para, degree=1)
      else:
          wpperiodic = WeakPseudoPeriodic_1c(mydomain, mri_para, degree=1)
      return wpperiodic
    
def convert_q2g(qvalue):
    g_ratio = 2.675e8;
    return qvalue*g_ratio*1e-12
  
def convert_g2q(gnorm):
    g_ratio = 2.675e8;
    return gnorm/g_ratio*1e12
