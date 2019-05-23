// This demo solves the Bloch-Torrey equation applied to computational diffusion MRI using
// the finite element method coupled with the theta-method for the spatial discretization.                                                                
// The scope of usage:                                                                                                                                    
//   (1) multilayered structures with permeability between intermediate interfaces
//   (2) pure homogeneous Neumann on boundaries
//   (3) Allow marking the phase function from a give submesh (mesh for cell).
//   (4) Allow water exchange at the external boundaries

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
// Last changed: 2019-03-25
//

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
#include <dolfin/mesh/IntersectionDetector.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/MeshData.h>
#include <dolfin/la/Matrix.h>
#include <dolfin/la/Vector.h>


#include <iomanip>
#include <iostream>
#include <sstream>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "./ufc/Bloch_Torrey_DG3D.h"
#include "./ufc/Bloch_Torrey3D.h"
#include "./ufc/Bloch_Torrey_NoTime3D.h"
#include "./ufc/Comp_Sig3D.h"
#include "./ufc/InitialCondition3D.h"


using namespace std;
using namespace dolfin;

int mpi_rank, mpi_nprocs;

void GatherBoundaryVertices(std::vector<int>&bdispl, Mesh*mesh, std::vector<double> &xcoors,std::vector<double> &ycoors,std::vector<double> &zcoors, std::vector<int>& boundarymeshsizes, Form* aM)
{
  MeshFunction<bool> used_vertex;
  used_vertex.init(*mesh,0);
  used_vertex = false;

  BoundaryMesh boundarymesh(*mesh);
  MeshFunction<uint>* cell_map = boundarymesh.data().meshFunction("cell map");

  boundarymeshsizes.resize(mpi_nprocs);
  std::vector<double> std_xcoors_local, std_ycoors_local, std_zcoors_local;

  int d = 4;
  UFC ufc(aM->form(), *mesh, aM->dofMaps());
  Cell c(*mesh, 0);
  uint local_dim = c.numEntities(0);
  uint *idx  = new uint[d * local_dim];
  double *XX_block = new double[d * local_dim];

  for(CellIterator bf(boundarymesh); !bf.end(); ++bf)
    {
      Facet f(*mesh, cell_map->get(*bf));

      for (CellIterator cell(f); !cell.end(); ++cell)
      //for (CellIterator cell(*mesh); !cell.end(); ++cell)
        {
          ufc.update(*cell, mesh->distdata());
          (aM->dofMaps())[0].tabulate_dofs(idx, ufc.cell, cell->index());

          for (VertexIterator v(*cell); !v.end(); ++v)
            {
              //if (!mesh->distdata().is_ghost(v->index(), 0) && !used_vertex.get(*v))
	      //if (!mesh->distdata().is_ghost(v->index(), 0))
                {
                  std_xcoors_local.push_back(v->x()[0]);
                  std_ycoors_local.push_back(v->x()[1]);
                  std_zcoors_local.push_back(v->x()[2]);
		  used_vertex.set(*v, true);
                }
            }
        }
    }

  double*xcoors_local = &std_xcoors_local[0];
  double*ycoors_local = &std_ycoors_local[0];
  double*zcoors_local = &std_zcoors_local[0];

  int lvertices=std_xcoors_local.size();
  MPI_Allgather(&lvertices, 1, MPI_INT, &boundarymeshsizes[0], 1, MPI_INT, dolfin::MPI::DOLFIN_COMM);
  int bnumvertices;
  MPI_Allreduce(&lvertices, &bnumvertices, 1, MPI_INT, MPI_SUM, dolfin::MPI::DOLFIN_COMM);
  xcoors.resize(bnumvertices);
  ycoors.resize(bnumvertices);
  zcoors.resize(bnumvertices);

  bdispl.resize(mpi_nprocs);
  bdispl[0]=0;
  for (int iproc=1; iproc<mpi_nprocs; iproc++)
  {
    bdispl[iproc] =bdispl[iproc-1] + boundarymeshsizes[iproc-1];
  }

  MPI_Allgatherv(xcoors_local,lvertices,MPI_DOUBLE,&xcoors[0],&boundarymeshsizes[0],&bdispl[0],MPI_DOUBLE,dolfin::MPI::DOLFIN_COMM); 
  MPI_Allgatherv(ycoors_local,lvertices,MPI_DOUBLE,&ycoors[0],&boundarymeshsizes[0],&bdispl[0],MPI_DOUBLE,dolfin::MPI::DOLFIN_COMM); 
  MPI_Allgatherv(zcoors_local,lvertices,MPI_DOUBLE,&zcoors[0],&boundarymeshsizes[0],&bdispl[0],MPI_DOUBLE,dolfin::MPI::DOLFIN_COMM);
} 



void get_diameter(Mesh*mesh, double&hmax, double&hmin)
{
  double hmax_local = 0;
  double hmin_local = 1e15;
 
  for (CellIterator cell(*mesh); !cell.end(); ++cell)
    {
      double d = cell->diameter();		
      hmax_local = max(d, hmax_local);
      hmin_local = min(d, hmin_local);
    }

  MPI_Allreduce(&hmax_local, &hmax, 1, MPI_DOUBLE,MPI_MAX,dolfin::MPI::DOLFIN_COMM);
  MPI_Allreduce(&hmin_local, &hmin, 1, MPI_DOUBLE,MPI_MIN,dolfin::MPI::DOLFIN_COMM);
}

Point gdir(1, 0, 0); 
Function u;
double xmin = -250., xmax=250., ymin=-500., ymax=200., zmin = -150., zmax=100.;
// double xmin = -10., xmax=10., ymin=-10., ymax=10., zmin = -10., zmax=10.;
// double xmin = -5., xmax=5., ymin=-5., ymax=5., zmin = -5., zmax=5.;

double ift = 0;


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
double IFT(double t, double delta, double Delta)
{
  double ft1 = t*(t>=0 && t<delta); 
  double ft2 = delta*(t>=delta && t<Delta); 
  double ft3 = (delta - t + Delta)*(t>=Delta && t<=Delta+delta); 
  return ft1 + ft2 + ft3;  
}

int Nsteps = 100;
double bvalue = 4000;
double delta = 40000;
double Delta = 40000;
double dt = 100;
double qvalue = sqrt(bvalue)/sqrt(delta*delta*(Delta-delta/3.0));
double g_ratio = 2.675e8;
double gvalue = qvalue/g_ratio*1e12;
double kappa = 5e-5;
double kcoeff = 3e-3;

int nskip = 5;
std::string dir="results";
int nrefine = 0;
bool is_dt = 0;
bool is_save = 0;

Form *L, *a_DG;

bool is_input_b = 0, is_input_q = 0;

Mesh *mesh, *cell_mesh;
bool existRegionFile = false;

std::string fmesh="mesh.xml", fcell="cell.xml";

void MarkPhase(Mesh*mesh, Mesh*cell_mesh, Vector*ei, MeshFunction<double>& meshfun)
{
  // MeshFunction<bool> solid_cells(mesh, mesh.topology().dim());
  meshfun.init(*mesh, mesh->topology().dim());

  bool existRegionFile = cell_mesh;
   
  IntersectionDetector *idetector;
   
  if (existRegionFile) 
  {
    message("Region file exists");
    idetector = new IntersectionDetector(*cell_mesh);

    for (CellIterator c(*mesh); !c.end(); ++c)
      {
	Cell& cell = *c;
	Point mp = cell.midpoint();
	int id = (*c).index();
	uint ci = id;
	if(dolfin::MPI::numProcesses() > 1)
	  ci = mesh->distdata().get_cell_global(ci);

	Array<unsigned int> overlap_cells;
	overlap_cells.clear();
	idetector->overlap(mp, overlap_cells);
     
	bool bfnd = false;
	for(int i=0; i < overlap_cells.size(); i++)
	  {	
	    Cell testcell(*cell_mesh, overlap_cells[i]);
	    if (cell_mesh->type().intersects(testcell,mp))
	      {
		bfnd = true;

		double value = 1.0;
		uint number = 1;
		ei->set(&value, number,&ci);
		meshfun.set(cell,value);

		break;
	      }			
	  }
      }
  }
  else
  {
    std::cout << "Region file doesn't exist" << std::endl;
    exit(0);
  }
}

class InitialVals: public Function
{
public:

  InitialVals(Mesh& mesh) : Function(mesh) {}

  void eval(double * value, const double* x) const
  {
    value[0] = 1.0;
    value[1] = 0.0;
    value[2] = 1.0;
    value[3] = 0.0;
  }

  uint rank() const
  {
    return 1;
  }

  uint dim(uint i) const
  {
    return 4;
  }

};

uint NumVertices(Mesh*mesh)
{
  uint numvertices_local = mesh->numVertices();
  uint numvertices;
  MPI_Allreduce(&numvertices_local, &numvertices, 1, MPI_UNSIGNED, MPI_SUM, dolfin::MPI::DOLFIN_COMM);
  return(numvertices);
}

int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);

  
  dolfin_set("output destination","silent");
  if(dolfin::MPI::processNumber() == 0)
    dolfin_set("output destination","terminal");
  

  double start, end;
  MPI_Barrier(dolfin::MPI::DOLFIN_COMM);
  start = MPI_Wtime();

  MPI_Comm_rank( dolfin::MPI::DOLFIN_COMM, &mpi_rank);
  MPI_Comm_size (dolfin::MPI::DOLFIN_COMM, &mpi_nprocs);


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
	      is_input_b = 1;
	      bvalue = atof(argv[optind+1]);
	      break;
	    case 'q':
	      is_input_q = 1;
	      gvalue = atof(argv[optind+1]);
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
	    case 'c':
	      {                                           
		fcell = argv[optind+1];   
		existRegionFile = true;
	      }
	      
	      break; 
	    case 's':
		is_save = atoi(argv[optind+1]);
		break; 
            case 'k':
                is_dt = 1;
                dt = atof(argv[optind+1]);
                break;
	    case 'v':
	      {
		double gx = atof(argv[optind+1]);
		double gy = atof(argv[optind+2]);
		double gz = atof(argv[optind+3]);
		gdir = Point(gx, gy, gz);
		gdir /= gdir.norm();
	      } 
	      break;
	    case 'K':
	      kcoeff = atof(argv[optind+1]);
	      break;
	    case 'p':
	      kappa = atof(argv[optind+1]);
	      break;
	    }
	}
    }

  if (is_input_b)
    {
      qvalue = sqrt(bvalue)/sqrt(delta*delta*(Delta-delta/3.0));
      gvalue = qvalue/g_ratio*1e12;
    }
  else if (is_input_q)
    {
      qvalue = gvalue*g_ratio*1e-12;
      bvalue = qvalue*qvalue*delta*delta*(Delta-delta/3.0);
    }

  if (mpi_rank==0)
    {
      std::string com = "mkdir -p "+dir;
      int check = system (com.c_str());
      message("Marking directory: %s", com.c_str());
    }

  // Read mesh
  message("\nReading mesh...");
  mesh = new Mesh(fmesh);

  if (existRegionFile)
    {
      message("Reading a given submesh");
      bool readInSerial = dolfin_get("Mesh read in serial");
      dolfin_set("Mesh read in serial",true);
      cell_mesh = new Mesh(fcell);
      dolfin_set("Mesh read in serial",readInSerial);
    }
  else
    {    
      message("Submesh is not given");
    }

  message("done\n");

  for (int i=0; i<nrefine; i++)
    {
      message("Refining %d time ...", i+1);
      mesh->refine();
    }


  message("Generating phase function...");
  a_DG = new Bloch_Torrey_DG3DBilinearForm();
  Function Phi; Vector Phi_vec;
  Phi.init(*mesh, Phi_vec, *a_DG, 0);

  MeshFunction<double> meshfun;

  if (existRegionFile)
    {
      MarkPhase(mesh, cell_mesh, &Phi_vec, meshfun);
    }

  File file_phi((dir+"/Phi.bin").c_str());
  file_phi << Phi;
  message("done");

  File ff((dir+"/Phi_viz.bin").c_str());
  ff << meshfun;
 
  File m_file((dir+"/mesh.bin").c_str());
  m_file<<*mesh;

  message("done");
  
  Assembler assembler(*mesh);
    
  int gdim = mesh->topology().dim();
    
  // # parameters
  // #####################################################
  double t = 0; double T = Delta+delta;
  double ft = 0; 

  if (!is_dt)
    dt = T/Nsteps;
  double theta = 0.5;

  message("kcoeff: %e, perm: %e, Nsteps: %d, dt: %f, delta: %f, Delta: %f, q: %e\n",kcoeff,kappa, Nsteps, dt, delta, Delta, qvalue);
  message("Gradient direction: %f %f %f",gdir.x(), gdir.y(), gdir.z());

  Function ft_f; PETScVector ft_fx;
  Function qvalue_f(*mesh, qvalue), kcoeff_f(*mesh, kcoeff), dt_f(*mesh, dt), theta_f(*mesh, theta), kappa_f(*mesh, kappa);
    
  GdotX GX(*mesh);


  dolfin_set("Krylov relative tolerance", 1e-4);
  dolfin_set("Krylov absolute tolerance", 1e-10);
  dolfin_set("Krylov maximum iterations", 100000);


  KrylovSolver sol(bicgstab, jacobi);
 

  Vector ux;
  double s0;

  Function up;   Vector upx;

  double hmax, hmin;
  get_diameter(mesh, hmax, hmin);
  message("hmax: %e, hmin: %e\n", hmax, hmin);

  MeshSize hf(*mesh);

  L = new Bloch_Torrey3DLinearForm(u, Phi, GX, ft_f, qvalue_f, kcoeff_f, theta_f, dt_f, kappa_f);

  // Initial conditions
  up.init(*mesh, upx, *L, 0);
  u.init(*mesh, ux, *L, 0);
  ft_f.init(*mesh, ft_fx, *L, 0);
  ft_fx = ft;

  message("Setting initial conditions ...");
  InitialVals initialvals(*mesh);

  InitialCondition3DBilinearForm ai(Phi);
  InitialCondition3DLinearForm Li(Phi, initialvals);

  PETScVector bi; PETScMatrix Ai;
  message("Assembling matrices ...");
  assembler.assemble(bi, Li);
  assembler.assemble(Ai, ai);
  sol.solve(Ai, ux , bi);
  message("done");

  File ui_file((dir+"/initial.bin").c_str());
  ui_file<<u;
    
  Comp_Sig3DFunctional S0(Phi,u);
  s0 = assembler.assemble(S0);
  message("s0=%f\n",s0);
    
  message("Preparing no-time matrices ...");
  // Construct time-indepent matrices
    
  PETScMatrix M, SI, J, MSI;
    
  Function mmk0(*mesh,1./dt),smk0(*mesh, 0.0),imk0(*mesh, 0.0),jmk0(*mesh, 0.0), zero(*mesh, 0.0);
  Bloch_Torrey_NoTime3DBilinearForm a_notime3D0(Phi, GX, kcoeff_f, kappa_f, mmk0, smk0, jmk0, imk0, zero, hf);
  assembler.assemble(M, a_notime3D0, true); 
  
  Function mmk1(*mesh,0.0),smk1(*mesh, theta),imk1(*mesh, theta),jmk1(*mesh, 0.0);                                                                              
  Bloch_Torrey_NoTime3DBilinearForm a_notime3D1(Phi, GX, kcoeff_f, kappa_f, mmk1, smk1, jmk1, imk1, zero, hf);
  assembler.assemble(SI, a_notime3D1, true); 
	
  MSI.dup(M);
  MSI += SI; // return 1./dt*M + theta*S
  MSI.apply();
  
  Function mmk2(*mesh,0.0),smk2(*mesh, 0.0),imk2(*mesh, 0.0),jmk2(*mesh, theta);                                                                                
  Bloch_Torrey_NoTime3DBilinearForm a_notime3D2(Phi, GX, kcoeff_f, kappa_f, mmk2, smk2, jmk2, imk2, zero, hf);
  assembler.assemble(J, a_notime3D2, true); // return theta*J
     
  File u_file((dir+"/u.bin").c_str());
  
  int step_counter = 0;
  while (t < T + dt)
    {
      message("t=%f, dt=%f, qvalue=%e, step_counter=%d, Completed %.1f%%\n",t, dt, qvalue, step_counter, t/T*100);
      ft =  FT(t, delta, Delta);
      ift =  IFT(t, delta, Delta);

      ft_fx = ft;
      
      PETScMatrix A, Jt;
      A.dup(MSI);
      Jt.dup(J);
      Jt *= ft*qvalue;
      Jt.apply();
      
      A += Jt;
      
      A.apply();

      if (up.vector().norm()>1e5)
	exit(0);
      message("Reassembling ...");
      Vector b;
      assembler.assemble(b, *L, true);
      message("done.\n");
      
      sol.solve(A, ux, b);
      if (step_counter%nskip==0 && is_save)
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
  
  message("b: %f, q: %e, g: %e, perm: %e, gdir: (%f, %f, %f), s: %e\n", bvalue, qvalue, gvalue, kappa, gdir.x(), gdir.y(), gdir.z(), s/s0);
  
  MPI_Barrier(dolfin::MPI::DOLFIN_COMM); 
  end = MPI_Wtime();
  
  if (mpi_rank == 0) { 
    message("Runtime = %f\n", end-start);
  }
  delete mesh, cell_mesh;
}
