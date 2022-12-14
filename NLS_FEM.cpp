// Radius of torus
real r = 2;

// Semi-major axis of ellipse
real a = 3;

// Semi-minor axis of ellipse (a>=b)
real b = 3;

// Toroidal curvature quantity
// phi -> x
// theta -> y
func gamma = sqrt( (a + r*cos(y))^2*sin(x)^2 + (b + r*cos(y))^2*cos(x)^2 );

// Define macros
// Gradient on the surface of the torus
macro grad(u) [gamma*dx(u), r*dy(u)]//EOM


// Discretization and end-time
// (Note dx is pre-defined)
verbosity = 0;
real Dx = .01, dt = 0.0005, idt = 1./dt, T = 0.08;

// Sign of non-linearity
real sigma = -1;

// Build the border of the mesh
border bottom(t = 0, 1){x = 2*pi*t; y = 0; label = 1;};
border top(t = 0, 1){x = -2*pi*t + 2*pi; y = 2*pi; label = 3;};
border left(t = 0, 1){x = 0; y = -2*pi*t + 2*pi; label = 4;};
border right(t = 0, 1){x = 2*pi; y = 2*pi*t; label = 2;};

// Enforcing periodic BCs
func periodicity = [[4,y],[2,y],[3,x],[1,x]];

// Discretize the domain
mesh Th = buildmesh(bottom(1/Dx) + top(1/Dx) + left(1/Dx) + right(1/Dx), fixedborder = true);
real mymeshsize = 0.0001;
Th = adaptmesh(Th, mymeshsize, IsMetric = 1, nbvx = 100000, periodic = periodicity); // periodic square
plot(Th, wait = true, cmm = "Refined mesh");

// Initialize finite element space type
fespace Vh(Th, P1, periodic = periodicity);
Vh<complex> u, v;

// Define mu
real mu = 20;

// Define IC
//func u0 = exp(-(x - pi)^2 - (y - pi)^2); // Gaussian
//func u0 = exp(-(x - pi)^2 - (y - pi)^2*1i) + exp(-(x - pi/1.2)^2*1i - (y - pi/1.2)^2); // Two Gaussians
//func u0 = 2 + sqrt(mu)*(tanh(sqrt(mu)*(x - pi)))^2 + 1i*sqrt(mu); // Dark soliton stripe
//func u0 = 2 + sqrt(mu)*tanh(sqrt(mu)*(sin(x)/4 - y + pi))^2; // A snaking dark soliton stripe
//func u0 = 2 + sqrt(mu)*tanh(sqrt(mu)*(sin(y)/4 - x + pi))^2; // Snaking stripe (other direction) GOOD FOR MESH MOVIE
func u0 = sqrt(mu)*tanh(sqrt(mu)*(sin(y)/4 - x + pi/2))^2 + sqrt(mu)*tanh(sqrt(mu)*(sin(y)/4 - x + 3*pi/2))^2; // two poloidal stripes
//func u0 = sqrt(mu)*tanh(sqrt(mu)*sqrt((x - 3*pi/2)^2 + (y - 3*pi/2)^2)); // one vortex
//func u0 = sqrt(mu)*tanh(sqrt(mu)*sqrt((x - 3*pi/2)^2 + (y - 3*pi/2)^2))*tanh(sqrt(mu)*sqrt((x - pi/2)^2 + (y - pi/2)^2)); // two vortices
//func u0 = sqrt(mu); // flat background
cout << "Initial mass is " << int2d(Th)(u0*conj(u0)) << endl;

// Initialize IC
Vh<complex> uold = u0, nonlin = u0*conj(u0);

// Adaptive mesh refinement based on IC
Vh fh = abs(u0), fhplot = fh^2;
plot(fhplot, wait = true, cmm = "Initial condition", fill = true, dim = 3, value = true);
Th = adaptmesh(Th, fh, periodic = periodicity); // periodic square
plot(Th, wait = true, cmm = "Adaptive mesh refinement");

// Weak form of NLS
problem NLS(u,v) = int2d(Th)(1i*idt*u*v) - int2d(Th)(1i*idt*uold*v)
									- int2d(Th)((grad(u)'*grad(v))/2)
									+ int2d(Th)(sigma*nonlin*u*v);

// Loop over time
real t = 0;
while (t <= T){
	t += dt;
	NLS;
	nonlin = u*conj(u);
	uold = u;
	plot(nonlin, dim = 2,
		cmm = "Mass = " + int2d(Th)(nonlin) + ", t = " + t, fill = 1, value = true);

	// Adaptive mesh refinement?
	Vh fh = abs(u);
	Th = adaptmesh(Th, fh, periodic = periodicity); // periodic square
	// Plot as a movie?
	//plot(Th, cmm = "Adaptive mesh refinement");
}
