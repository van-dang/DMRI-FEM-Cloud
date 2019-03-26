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
//

#include <dolfin.h>
#include <mshr.h> 
#include "./ufc/Bloch_Torrey.h"
#include "./ufc/Bloch_Torrey_NoTime.h"
#include "./ufc/Comp_Sig.h"
#include "./ufc/Bloch_Torrey3D.h"
#include "./ufc/Bloch_Torrey_NoTime3D.h"
#include "./ufc/Comp_Sig3D.h"

using namespace dolfin;

Point gdir(0, 0, 1); 

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

// Matrix ThetaMethod_A(double ft, double gnorm, double theta, double dt, Matrix M, Matix J, Matix S, Matrix I):
//    return 1./dt*M + ft*gnorm*theta*J + theta*S + theta*I
int Nsteps = 100;
double bvalue = 1000;
double delta = 1000;
double Delta = 1000;
int nskip = 2;

// Bloch_Torrey::LinearForm*L;
Form *L;
std::string fmesh="mesh.xml", fphase="Phi.xml";

std::shared_ptr<FunctionSpace> V_DG, V, W;

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
	      case 'p':
		fphase = argv[optind+1];
		break;		
	      } 
	  } 
      }

    /*
    printf("\n Generating the mesh...\n");
    int is_run_geom = system("python parameters_multi_comp_3D.py");
    if (is_run_geom != 0)
    {
        printf("Error with geometry!\n");
        exit(0);
    } */
    
    // Read mesh
    printf("\nReading mesh...");
    auto mesh = std::make_shared<Mesh>(fmesh);
    printf("done\n");

    int gdim = mesh->topology().dim();
    
    // Create velocity FunctionSpace
    if (gdim==2)
	V_DG = std::make_shared<Bloch_Torrey::CoefficientSpace_phase>(mesh);
     
   if (gdim==3)
	V_DG = std::make_shared<Bloch_Torrey3D::CoefficientSpace_phase>(mesh);

   // Reading
   printf("Reading phase function...");
   auto Phi = std::make_shared<Function>(V_DG);
   File file_phi(fphase);
   file_phi >> *Phi;
   printf("done\n");
   // plot(Phi); interactive(); exit(0);

   // Create function space
   if (gdim==2)
     {
       W = std::make_shared<Bloch_Torrey::FunctionSpace>(mesh);
       V = std::make_shared<Bloch_Torrey::CoefficientSpace_u0r_0>(mesh);
     }
   if (gdim==3)
      {
	W = std::make_shared<Bloch_Torrey3D::FunctionSpace>(mesh);
	V = std::make_shared<Bloch_Torrey3D::CoefficientSpace_u0r_0>(mesh);
      } 
    // #################################################################################
    // # Initial conditions
    auto one = std::make_shared<Constant>(1.0);
    auto zero = std::make_shared<Constant>(0.0);
    auto gx = std::make_shared<GdotX>();
    
    // # parameters
    // #####################################################
    double t=0; double T = Delta+delta;
    double ft = 0; double kcoeff = 3e-3;  
    double kappa = 5e-5;
    double dt = T/Nsteps; double theta = 0.5;
    double gnorm = sqrt(bvalue)/sqrt(delta*delta*(Delta-delta/3.0));
    // #####################################################
    
    auto ft_f = std::make_shared<Constant>(ft);
    auto gnorm_f = std::make_shared<Constant>(gnorm);
    auto kcoeff_f = std::make_shared<Constant>(kcoeff);
    auto kappa_f = std::make_shared<Constant>(kappa);
    auto dt_f = std::make_shared<Constant>(dt);
    auto theta_f = std::make_shared<Constant>(theta);

    auto u0r_0 = std::make_shared<Function>(V);
    auto u0i_0 = std::make_shared<Function>(V);
    auto u1r_0 = std::make_shared<Function>(V);
    auto u1i_0 = std::make_shared<Function>(V);
    
    *u0r_0->vector() = 1;
    *u0i_0->vector() = 0;
    *u1r_0->vector() = 1;
    *u1i_0->vector() = 0;
    Bloch_Torrey::LinearForm*L2D;
    Bloch_Torrey3D::LinearForm*L3D;
    printf("Generate linear forms...");
    if (gdim==2)
      {
	L2D = new Bloch_Torrey::LinearForm(W);
	L2D->phase = Phi;
	L2D->u0r_0 = u0r_0;
	L2D->u0i_0 = u0i_0;
	L2D->u1r_0 = u1r_0;
	L2D->u1i_0 = u1i_0;
	L2D->GX = gx;
	L2D->ft = ft_f;
	L2D->gnorm = gnorm_f;
	L2D->K = kcoeff_f;
	L2D->kappa = kappa_f;
	L2D->theta = theta_f;
	L2D->dt = dt_f;
	L = L2D;
      }
    if (gdim==3)
      {
	L3D = new Bloch_Torrey3D::LinearForm(W);
	L3D->phase = Phi;
	L3D->u0r_0 = u0r_0;
	L3D->u0i_0 = u0i_0;
	L3D->u1r_0 = u1r_0;
	L3D->u1i_0 = u1i_0;
	L3D->GX = gx;
	L3D->ft = ft_f;
	L3D->gnorm = gnorm_f;
	L3D->K = kcoeff_f;
	L3D->kappa = kappa_f;
	L3D->theta = theta_f;
	L3D->dt = dt_f;
	L = L3D;
      }
    printf("done\nPreparing no-time matrices ...");
    // Construct time-indepent matrices
    std::shared_ptr<Matrix> M(new Matrix), SI(new Matrix), J(new Matrix);
    std::shared_ptr<Matrix> MSI(M);
    double s0;

    if (gdim==2)
      {
	Bloch_Torrey_NoTime::BilinearForm a_notime2D(W, W);
	a_notime2D.phase = Phi;
	a_notime2D.K = kcoeff_f;
	a_notime2D.kappa = kappa_f;
	a_notime2D.GX = gx;
    
	a_notime2D.mmk = std::make_shared<Constant>(1./dt);
	a_notime2D.smk = zero;
	a_notime2D.imk = zero;
	a_notime2D.jmk = zero;
	assemble(*M, a_notime2D); 
	M->ident_zeros();
    
	a_notime2D.mmk = zero;
	a_notime2D.smk = theta_f;
	a_notime2D.imk = theta_f;
	a_notime2D.jmk = zero;
	assemble(*SI, a_notime2D);
    
    
	*MSI = *M;
	*MSI += *SI; // return 1./dt*M + theta*S + theta*I
	MSI->apply("add");
          
	a_notime2D.mmk = zero;
	a_notime2D.smk = zero;
	a_notime2D.imk = zero;
	a_notime2D.jmk = theta_f;
	assemble(*J, a_notime2D); // return theta*J

	Comp_Sig::Functional S0(mesh, u0r_0, u1r_0, Phi);
	s0 = assemble(S0);
      }

    if (gdim==3)
      {
	Bloch_Torrey_NoTime3D::BilinearForm a_notime3D(W, W);
	a_notime3D.phase = Phi;
	a_notime3D.K = kcoeff_f;
	a_notime3D.kappa = kappa_f;
	a_notime3D.GX = gx;
    
	a_notime3D.mmk = std::make_shared<Constant>(1./dt);
	a_notime3D.smk = zero;
	a_notime3D.imk = zero;
	a_notime3D.jmk = zero;
	assemble(*M, a_notime3D); 
	M->ident_zeros();
    
	a_notime3D.mmk = zero;
	a_notime3D.smk = theta_f;
	a_notime3D.imk = theta_f;
	a_notime3D.jmk = zero;
	assemble(*SI, a_notime3D);
    
    
	*MSI = *M;
	*MSI += *SI; // return 1./dt*M + theta*S + theta*I
	MSI->apply("add");
          
	a_notime3D.mmk = zero;
	a_notime3D.smk = zero;
	a_notime3D.imk = zero;
	a_notime3D.jmk = theta_f;
	assemble(*J, a_notime3D); // return theta*J

	Comp_Sig3D::Functional S0(mesh, u0r_0, u1r_0, Phi);
	s0 = assemble(S0);
      }
    printf("done\n");

    printf("Intial signal: %f\n",s0);

    auto u = std::make_shared<Function>(W);
    int step_count = 0;
    while (t < T)
    {
        if (step_count % nskip == 0)
	  printf("t=%f, dt=%f, gnorm: %f\n",t, dt,gnorm);
	
        ft =  FT(t, delta, Delta);
        *ft_f = Constant(ft);
        
        Matrix A = *MSI;
        Matrix Jt = *J;
        Jt *= ft*gnorm;
        A += Jt;
	
	Vector b;
        assemble(b, *L);
        solve(A,*u->vector(), b,"gmres");        
     
        u0r_0->interpolate((*u)[0]);
        u0i_0->interpolate((*u)[1]);
        u1r_0->interpolate((*u)[2]);
        u1i_0->interpolate((*u)[3]); 
        
        t += dt;
	step_count += 1;
    }

    u0r_0->interpolate((*u)[0]);
    u0i_0->interpolate((*u)[1]);
    u1r_0->interpolate((*u)[2]);
    u1i_0->interpolate((*u)[3]);
    printf("delta: %f, Delta: %f, Nsteps: %d, dt: %f\n",delta, Delta, Nsteps, dt);
    printf("bvalue: %f, gnorm: %f, kcoeff: %e, kappa: %e\n", bvalue, gnorm, kcoeff, kappa);

    double st;
    if (gdim==2)
      {
	Comp_Sig::Functional St(mesh, u0r_0, u1r_0, Phi);
	st = assemble(St);
      }
    if (gdim==3)
      {
	Comp_Sig3D::Functional St(mesh, u0r_0, u1r_0, Phi);
	st = assemble(St);	
      }
    printf("bvalue: %f,signal %f\n", bvalue, st/s0);

    // Memory release
    if (gdim==2)
      {
	delete L2D;
      }
    if (gdim==3)
      {
	delete L3D;
      }
    //delete a_notime;
}
