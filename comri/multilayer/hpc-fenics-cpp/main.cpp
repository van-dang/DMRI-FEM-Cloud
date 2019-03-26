// This demo solves the Bloch-Torrey equation applied to computational diffusion MRI using
// the finite element method coupled with the theta-method for the spatial discretization.                                                                
// The scope of usage:                                                                                                                                    
//   (1) multilayered structures with permeability between intermediate interfaces
//   (2) pure homogeneous Neumann on boundaries
// Run on beskow: 
//   aprun -n 16 ./demo -m multi_layer_torus.xml -d 10000 -D 10000 -N 1000 -b 1000 -r 1

// Copyright (C) 2017 Van-Dang Nguyen

// This file is part of DOLFIN.
//
// DOLFIN is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// DOLFIN is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
//
// First added:  2017-10-10
// Last changed: 2017-10-10


#include <dolfin/config/dolfin_config.h>
#include <dolfin/mesh/RivaraRefinement.h>
#include <dolfin/main/init.h>
#include <dolfin/parameter/parameters.h>
#include <dolfin/common/constants.h>
#include <dolfin/common/TimeDependent.h>
#include <dolfin/fem/DirichletBC.h>
#include <dolfin/la/KrylovSolver.h>
#include <dolfin/la/LUSolver.h>
#include <dolfin/function/Function.h>
#include <dolfin/function/SpecialFunctions.h>
#include <dolfin/fem/UFC.h>
#include <dolfin/fem/Assembler.h>
#include <dolfin/mesh/Vertex.h>

#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/MeshData.h>
#include <dolfin/la/Matrix.h>
#include <dolfin/la/Vector.h>

#include <iomanip>
#include <iostream>
#include <sstream>


#include "./ufc/Bloch_Torrey_DG3D.h"
#include "./ufc/Bloch_Torrey3D.h"
#include "./ufc/Bloch_Torrey_NoTime3D.h"
#include "./ufc/Comp_Sig3D.h"
#include "./ufc/InitialCondition3D.h"

using namespace dolfin;

Point gdir(0, 0, 1); 

class GdotX : public Function
{
  public:
    
    GdotX(Mesh& mesh) : Function(mesh) {}
    
    void eval(double * value, const double* x) const
    {
      value[0] = x[0]*gdir.x()+x[1]*gdir.y()+x[2]*gdir.z();
    }

    uint rank() const
    {
      return 0;
    }

    uint dim(uint i) const
    {
      return 1;
    }

};


double FT(double t, double delta, double Delta)
{
    double ft1 = 1.0*(t>=0 && t<delta); 
    double ft2 = -1.0*(t>=Delta && t<=Delta+delta);
    return ft1 + ft2;  
}

int Nsteps = 100;
double bvalue = 1000;
double delta = 1000;
double Delta = 1000;
int nskip = 5;
std::string dir="results";
int nrefine = 0;

Form *L, *a_DG;

std::string fmesh="mesh.xml", fphase="Phi.xml";

int ncomps = 3;
double R0 = 5.0;
double r[] = {R0, 1.5*R0,2*R0};
double R = 20;

void MarkPhase(Mesh*mesh, Vector*ei, MeshFunction<real>& meshfun )
{
  meshfun.init(*mesh, mesh->topology().dim());

  for (CellIterator cell(*mesh); !cell.end(); ++cell)
    {
      int id = (*cell).index();
      uint ci = id;
      dolfin::Point p = cell->midpoint();
      if(dolfin::MPI::numProcesses() > 1)
        ci = mesh->distdata().get_cell_global(ci);

      double d = pow(R-sqrt(p.x()*p.x()+p.y()*p.y()),2)+p.z()*p.z();
      for (int i=0; i<ncomps; i++)
	{
	  if (d<r[i]*r[i])
	    {
	      real value = i%2;
	      uint number = 1;
	      ei->set(&value, number,&ci);
	      meshfun.set(*cell,value);
	      break;
	    }
	}
    }
}

int main(int argc, char *argv[])
{

  MPI_Init(&argc, &argv);
  dolfin_set("output destination","silent");
  if(dolfin::MPI::processNumber() == 0)
    dolfin_set("output destination","terminal");

  double start, end;
  MPI_Barrier(MPI_COMM_WORLD); 
  start = MPI_Wtime();

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  

  message("%d arguments\n",argc);
  message("Command: %s", argv[0]);
  for (int optind = 1; optind < argc; optind++)
    {
	message(" %s ",argv[optind]);
	if (argv[optind][0] == '-')
	  {
	    switch (argv[optind][1])
	      {
	      case 'N':
		Nsteps = atoi(argv[optind+1]);
		break;
	      case 'b':
		bvalue = atof(argv[optind+1]);
		break;
	      case 'd':
		delta = atof(argv[optind+1]);
		break;
	      case 'D':
		Delta = atof(argv[optind+1]);
		break;
	      case 'j':
		nskip = atoi(argv[optind+1]);
		break;
	      case 'm':
		fmesh = argv[optind+1];
		break;
              case 'f':
                dir = argv[optind+1];
                break;
              case 'r':
                nrefine = atoi(argv[optind+1]);
                break;
	      } 
	  } 
    }

  // Read mesh
  message("\nReading mesh...");
  Mesh mesh (fmesh);
  message("done\n");

  for (int i=0; i<nrefine; i++)
    {
      message("Refining %d time ...", i+1);
      mesh.refine();
    }

  printf("Generating phase function...");
  a_DG = new Bloch_Torrey_DG3DBilinearForm();
  Function Phi; Vector Phi_vec;
  Phi.init(mesh, Phi_vec, *a_DG, 0);

  MeshFunction<real> meshfun;
  MarkPhase(&mesh, &Phi_vec, meshfun);

  File file_phi((dir+"/Phi.bin").c_str());
  file_phi << Phi;
  printf("done\n");

  File ff((dir+"/Phi_viz.bin").c_str());
  ff << meshfun;
 
  Assembler assembler(mesh);
    
  int gdim = mesh.topology().dim();
    
  // # parameters
  // #####################################################
  double t = 0; double T = Delta+delta;
  double ft = 0; double kcoeff = 3e-3; double kappa = 5e-5;  
  double dt = T/Nsteps; double theta = 0.5;
  double gnorm = sqrt(bvalue)/sqrt(delta*delta*(Delta-delta/3.0));
  // #####################################################
  Function ft_f; PETScVector ft_fx;
  Function gnorm_f(mesh, gnorm), kcoeff_f(mesh, kcoeff), dt_f(mesh, dt), theta_f(mesh, theta), kappa_f(mesh, kappa);
    
  GdotX GX(mesh);

  KrylovSolver sol(bicgstab, none);
 
  Function u;
  Vector ux;
  double s0;

  L = new Bloch_Torrey3DLinearForm(u, Phi, GX, ft_f, gnorm_f, kcoeff_f, theta_f, dt_f, kappa_f);

  // Initial conditions
  u.init(mesh, ux, *L, 0);
  ft_f.init(mesh, ft_fx, *L, 0);
  ft_fx = ft;

  message("Setting initial conditions ...");
  InitialCondition3DBilinearForm ai(Phi);
  InitialCondition3DLinearForm Li(Phi);
  PETScVector bi; PETScMatrix Ai;
  message("Assembling matrices ...");
  assembler.assemble(bi, Li);
  assembler.assemble(Ai, ai);
  sol.solve(Ai, ux , bi);
  message("done");

  Comp_Sig3DFunctional S0(Phi,u);
  s0 = assembler.assemble(S0);
  message("s0=%f\n",s0);
    
  message("Preparing no-time matrices ...");
  // Construct time-indepent matrices
    
  PETScMatrix M, SI, J, MSI;
    
  if (gdim==3)
    {
      Function mmk0(mesh,1./dt),smk0(mesh, 0.0),imk0(mesh, 0.0),jmk0(mesh, 0.0);
      Bloch_Torrey_NoTime3DBilinearForm a_notime3D0(Phi, GX, kcoeff_f, kappa_f, mmk0, smk0, jmk0, imk0);
      assembler.assemble(M, a_notime3D0, true); 
	 
      Function mmk1(mesh,0.0),smk1(mesh, theta),imk1(mesh, theta),jmk1(mesh, 0.0);                                                                              
      Bloch_Torrey_NoTime3DBilinearForm a_notime3D1(Phi, GX, kcoeff_f, kappa_f, mmk1, smk1, jmk1, imk1);
      assembler.assemble(SI, a_notime3D1, true); 
	
      MSI.dup(M);
      MSI += SI; // return 1./dt*M + theta*S
      MSI.apply();
        
      Function mmk2(mesh,0.0),smk2(mesh, 0.0),imk2(mesh, 0.0),jmk2(mesh, theta);                                                                                
      Bloch_Torrey_NoTime3DBilinearForm a_notime3D2(Phi, GX, kcoeff_f, kappa_f, mmk2, smk2, jmk2, imk2);
      assembler.assemble(J, a_notime3D2, true); // return theta*J
    }  
    
  if (rank==0)
    {
      std::string com = "mkdir -p "+dir;
      system (com.c_str());
      message("Marking directory: %s", com.c_str());
    }
  
  File m_file((dir+"/mesh.bin").c_str());
  m_file<<mesh;
  
  File u_file((dir+"/u.bin").c_str());
  
  int step_counter = 0;
  while (t < T)
    {
      message("t=%f, dt=%f, gnorm=%f, step_counter=%d, Completed %.1f%%\n",t, dt, gnorm, step_counter, t/T*100);
      ft =  FT(t, delta, Delta);
      ft_fx = ft;
      
      PETScMatrix A, Jt;
      A.dup(MSI);
      Jt.dup(J);
      Jt *= ft*gnorm;
      Jt.apply();
      
      A += Jt;
      A.apply();
      
      message("u.vector().norm(l2): %f\n", u.vector().norm());
      message("Reassembling ...");
      Vector b;
      assembler.assemble(b, *L, true);
      message("done.\n");
      
      sol.solve(A, ux, b);
      if (step_counter%nskip==0)
	{
	  message("Write to %s", dir.c_str());
	  u_file<<u;
	}
      
      t += dt;
      b = 0;
      step_counter += 1;
    }
  message("done\n");
  Comp_Sig3DFunctional S(Phi,u);
  double s = assembler.assemble(S);
  message("s=%f\n",s/s0);
  
  MPI_Barrier(MPI_COMM_WORLD); 
  end = MPI_Wtime();
  
  if (rank == 0) { 
    message("Runtime = %f\n", end-start);
  }
}
