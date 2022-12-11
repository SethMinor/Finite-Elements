/*
--------------------------------------------------------------------------------
Author : Georges SADAKA - LAMFA - AMIENS
http://www.lamfa.u-picardie.fr/sadaka/FreeFem++/GTA3/GTA3FreeFem++2D.html

Ionut Danaila:
http://ionut.danaila.perso.math.cnrs.fr/zdownload/Tutorial_2019_Singapore/Singapore_Day_01.pdf
http://ionut.danaila.perso.math.cnrs.fr/zdownload/Tutorial_2019_Singapore/Singapore_Day_02.pdf

Periodic BC's
https://www.um.es/freefem/ff++/pmwiki.php?n=Main.PeriodicBoundaryConditions

Metric tensor
http://ionut.danaila.perso.math.cnrs.fr/zdownload/Tutorial_2016_Fields/Tutorial_2016_Fields_Course_05_to_08.pdf
--------------------------------------------------------------------------------
*/

// Define macros
// Eventually replace these guys with Laplace-Beltrami operators
macro grad(u) [dx(u), dy(u)]//EOM

// Discretization and end-time
// (Note dx is pre-defined)
verbosity = 0;
real Dx = .01, dt = 0.001, idt = 1./dt, T = 0.5;

// Sign of non-linearity
real sigma = -1.;

// Build the border of the mesh
// Parametric borders must have continuous arrows
border bottom(t = 0, 1){x = 2*pi*t; y = 0; label = 1;};
border top(t = 0, 1){x = -2*pi*t + 2*pi; y = 2*pi; label = 3;};
border left(t = 0, 1){x = 0; y = -2*pi*t + 2*pi; label = 4;};
border right(t = 0, 1){x = 2*pi; y = 2*pi*t; label = 2;};

// Enforcing periodic BCs
func perio = [[4,y],[2,y],[3,x],[1,x]];

// Discretize the domain
// Adpative mesh refinement (number = constant mesh size, nbvx = max number of vertices)
// Can adapt the mesh in a loop over time, based on value of u:
// https://doc.freefem.org/documentation/mesh-generation.html#the-command-adaptmesh
mesh Th = buildmesh(bottom(1/Dx) + top(1/Dx) + left(1/Dx) + right(1/Dx), fixedborder = true);
real mymeshsize = 0.05;
Th = adaptmesh(Th, mymeshsize, IsMetric = 1, nbvx = 100000, periodic = perio); // periodic square
plot(Th, wait = true, cmm = "Refined mesh");

// Initialize finite element space type
fespace Vh(Th, P2, periodic = perio);
Vh<complex> u, v;

// Define mu
real mu = 20;

// Define IC
//func u0 = -exp(-(x - pi)^2-(y - pi)^2-5); // Gaussian
//func u0 = sqrt(mu)*tanh(sqrt(mu)*sqrt((x - pi)^2 + (y - pi)^2)); // one vortex
func u0 = sqrt(mu)*tanh(sqrt(mu)*sqrt((x - 3*pi/2)^2 + (y - 3*pi/2)^2)) *
					tanh(sqrt(mu)*sqrt((x - pi/2)^2 + (y - pi/2)^2)); // two vortices
//func u0 = sqrt(mu); // flat background
cout << "Initial mass is " << int2d(Th)((u0*conj(u0))^2) << endl;

// Initialize IC
Vh<complex> uold = u0, nonlin = u0*conj(u0);

// Adaptive mesh refinement based on IC
Vh fh = abs(u0);
plot(fh, wait = true, cmm = "Initial condition", fill = true, dim = 3);
Th = adaptmesh(Th, fh, periodic = perio); // periodic square
plot(Th, wait = true, cmm = "Adaptive mesh refinement");

// Weak form of NLS
problem NLS(u,v) = int2d(Th)(1i*idt*u*v) - int2d(Th)(1i*idt*uold*v)
									+ int2d(Th)(sigma*nonlin*u*v)
									+ int2d(Th)((grad(u)'*grad(v))/2);

// Loop over time
real t = 0;
while (t <= T){
	t += dt;
	NLS;
	nonlin = u*conj(u);
	uold = u;
	plot(nonlin, dim = 3,
		cmm = "t = " + t + "  ; \nMass = " + int2d(Th)(nonlin), fill = 1);

	// Adaptive mesh refinement?
	Vh fh = abs(u);
	Th = adaptmesh(Th, fh, periodic = perio); // periodic square
	// Plot as a movie?
	//plot(Th, cmm = "Adaptive mesh refinement");
}
