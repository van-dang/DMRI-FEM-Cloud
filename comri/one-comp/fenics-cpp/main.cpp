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
// Last changed: 2019-03-16
//

#include <iomanip>
#include <iostream>
#include <sstream>

#include <dolfin.h>

#include "./ufc/Bloch_Torrey3D.h"
#include "./ufc/Bloch_Torrey_NoTime3D.h"
#include "./ufc/Comp_Sig3D.h"
#include "./ufc/InitialCondition3D.h"

using namespace dolfin;

Point gdir(0, 1, 0);

class GdotX : public Expression
{
  void eval(Array<double>& values, const Array<double>& x) const
  {
    Point px(x[0], x[1], x[2]);
    values[0] = px.dot(gdir);
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
double delta = 40000;
double Delta = 40000;
double dt = 100;
double gnorm = sqrt(bvalue)/sqrt(delta*delta*(Delta-delta/3.0));
double g_ratio = 2.675e8;
double qvalue = gnorm/g_ratio*1e12;
int nskip = 5;
bool is_save = 0;
bool is_dt = 0;
std::string dir="results";
int nrefine = 0;

std::string fmesh="mesh.xml";

bool is_input_b = 0, is_input_q = 0;
int main(int argc, char *argv[])
{
  printf("%d arguments\n",argc);
  printf("Command: %s", argv[0]);
  for (int optind = 1; optind < argc; optind++)
    {
	printf(" %s ",argv[optind]);
	if (argv[optind][0] == '-')
	  {
	    switch (argv[optind][1])
	      {
	      case 'N':
		Nsteps = atoi(argv[optind+1]);
		break;
	      case 'b':
		bvalue = atof(argv[optind+1]);
		is_input_b = 1;
		break;
	      case 'q':
		qvalue = atof(argv[optind+1]);
		is_input_q = 1;
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
	      case 's':
		is_save = atoi(argv[optind+1]);
		break;
	      case 'k':
		is_dt = 1;
		dt = atof(argv[optind+1]);
		break;
	      case 'v':
                double gx = atof(argv[optind+1]);
                double gy = atof(argv[optind+2]);
                double gz = atof(argv[optind+3]);
		gdir = Point(gx, gy, gz);
		gdir /= gdir.norm();
                break;
	      } 
	  } 
    }

  if (is_input_b)
    {
      gnorm = sqrt(bvalue)/sqrt(delta*delta*(Delta-delta/3.0));
      qvalue = gnorm/g_ratio*1e12;
    }
  else if (is_input_q)
    {
      gnorm = qvalue*g_ratio*1e-12;
      bvalue = gnorm*gnorm*delta*delta*(Delta-delta/3.0);
    }
  // Read mesh
  printf("\nReading mesh...");
  auto mesh = std::make_shared<Mesh>(fmesh);
  //Mesh mesh (fmesh);
  printf("done\n");

  for (int i=0; i<nrefine; i++)
    {
      printf("Refining %d time ...", i+1);
      // *mesh = refine(*mesh);
    }
    
  int gdim = mesh->topology().dim();
    
  // # parameters
  // #####################################################
  double t = 0; double T = Delta+delta;
  double ft = 0; double kcoeff = 2.4e-3;
  
  if (!is_dt)
    dt = T/Nsteps;
  double theta = 0.5;
  printf("\ndiffusion coefficient: %e, Nsteps: %d, dt: %f, delta: %f, Delta: %f, gnorm: %f\n",kcoeff, Nsteps, dt, delta, Delta, gnorm);
  printf("\nGradient direction: %f %f %f",gdir.x(), gdir.y(), gdir.z());

  // #####################################################
  // Bloch_Torrey3D::FunctionSpace V(mesh);
  auto W = std::make_shared<Bloch_Torrey3D::FunctionSpace>(mesh);
  
  
  // Function ft_f(W);
  // ft_f.vector() = ft;

  auto ft_f = std::make_shared<Function>(W);
  *ft_f->vector() = ft;

  auto gnorm_f = std::make_shared<Constant>(gnorm);
  auto kcoeff_f = std::make_shared<Constant>(kcoeff);
  auto dt_f = std::make_shared<Constant>(dt);
  auto theta_f = std::make_shared<Constant>(theta);
  auto GX = std::make_shared<GdotX>();
  
  KrylovSolver sol("gmres","ilu");
  sol.parameters["absolute_tolerance"] = 1e-6;
  sol.parameters["relative_tolerance"] = 1e-4;
  sol.parameters["maximum_iterations"] = 1000;

  double s0;

  auto u = std::make_shared<Function>(W);
  Bloch_Torrey3D::LinearForm L(W);
  L.u_0 = u;
  L.GX = GX;
  L.ft_a = ft_f;
  L.gnorm = gnorm_f;
  L.K = kcoeff_f;
  L.theta = theta_f;
  L.dt = dt_f;
  
  printf("\nSetting initial conditions ...");
  InitialCondition3D::BilinearForm ai(W, W);
  InitialCondition3D::LinearForm Li(W);
  Vector bi; Matrix Ai;
  assemble(bi, Li);
  assemble(Ai, ai);
  sol.solve(Ai, *u->vector() , bi);
  printf("done\n");
  
  Comp_Sig3D::Functional S0(mesh, u);
  s0 = assemble(S0);
  printf("\ns0=%f\n",s0);
    
  printf("Preparing no-time matrices ...");
  // Construct time-indepent matrices
    
  // PETScMatrix M, SI, J;
  std::shared_ptr<Matrix> M(new Matrix), SI(new Matrix), J(new Matrix);

  auto zero = std::make_shared<Constant>(0.0);
  
  if (gdim==3)
    {
      Bloch_Torrey_NoTime3D::BilinearForm a_notime3D(W,W);
      auto mmk0 = std::make_shared<Constant> (1./dt);
      a_notime3D.GX = GX;
      a_notime3D.K = kcoeff_f;
      a_notime3D.mmk = mmk0;
      a_notime3D.smk = zero;
      a_notime3D.jmk = zero;
      assemble(*M, a_notime3D); 

      a_notime3D.mmk = zero;
      a_notime3D.smk = theta_f;
      a_notime3D.jmk = zero;
      assemble(*SI, a_notime3D);
      
      a_notime3D.mmk = zero;
      a_notime3D.smk = zero;
      a_notime3D.jmk = theta_f;
      assemble(*J, a_notime3D);
    }

  std::shared_ptr<Matrix> MSI(M);
  *MSI = *M;
  *MSI += *SI; // return 1./dt*M + theta*S + theta*I
  MSI->apply("add");
      
  
  std::string com = "mkdir -p "+dir;
  int check = system (com.c_str());
  printf("\nMarking directory: %s\n", com.c_str());
  
  File u_file((dir+"/u.pvd").c_str());
  
  int step_counter = 0;
  while (t < T + dt)
    {
      if (step_counter%nskip==0)
	printf("t=%f, dt=%f, gnorm=%f, step_counter=%d, Completed %.1f%%\n",t, dt, gnorm, step_counter, t/T*100);
      ft =  FT(t, delta, Delta);
      *ft_f->vector() = ft;
      
      Matrix A = *MSI;
      Matrix Jt = *J;
      Jt *= ft*gnorm;
      Jt.apply("add");
      
      A += Jt;
      A.apply("add");
      
      Vector b;
      assemble(b, L);
      
      sol.solve(A, *u->vector(), b);
      if (step_counter%nskip==0 && is_save)
	{
	  printf("Write to %s", dir.c_str());
	  u_file<<*u;
	}
      
      t += dt;
      b = 0;
      step_counter += 1;
    }
  printf("done\n");
  Comp_Sig3D::Functional S(mesh, u);
  double s = assemble(S);
  printf("b: %f, gnorm: %f, q: %f, gdir: (%f, %f, %f), s: %f\n", bvalue, gnorm, qvalue, gdir.x(), gdir.y(), gdir.z(), s/s0);
}
