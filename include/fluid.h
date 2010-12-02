#ifndef FLUID_H
	#define FLUID_H
	using namespace std;
	
	typedef double Real;
	
	// the weight for the equilibrium distribution according with Hecht
	Real w[19] = {(2./36.),(2./36.),(2./36.),(2./36.),(2./36.),(2./36.),
	      (1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),
	      (1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),
	      (12./36.)};

	float e_x[19] = {1.,-1.,0.,0.,0.,0.,1.,1.,1.,1.,-1.,-1.,-1.,-1.,0.,0.,0.,0.,0.};
	float e_y[19] = {0.,0.,1.,-1.,0.,0.,1.,-1.,0.,0.,1.,-1.,0.,0.,1.,1.,-1.,-1.,0.};
	float e_z[19] = {0.,0.,0.,0.,1.,-1.,0.,0.,1.,-1.,0.,0.,1.,-1.,1.,-1.,1.,-1.,0.};
	int dfInv[19] = {1,0,3,2,5,4,11,10,13,12,7,6,9,8,17,16,15,14,18};
	
	// different types of cells
	const int FLUIDO  = 3, TOP = 1, BOTTOM  = 2, NOSLIP = 0;
	// velocity for the moving wall
	const Real V = 0.00, W = 0.00, cs=1.0/sqrt(3.0);
	
	// Declare some constants and globals ...
	const Real   omega = 1.00;     // viscosity of the fluid, 0..2
	int current = 0, other = 1; // which grid is the current and the old one?
	
	class fluid{

	private:

		int ts;
		int X, Y, Z;
		Real *****cells;
		Real ***flags;	
		Real ****vel;
		Real ***rho;
		Real ****fuerza;
		Real U;
		
	public:
		
		void inicializar(int x, int y, int z);
		void stream();
		void collide();
		int guardar(int s);
		float darVelocidad(int x, int y, int z, int f);
		void setFuerza(int x, int y, int z, float f[3]);
		void addFuerza(int x, int y, int z, float f[3]);
		void calcularMacro();
		void setVelocidad(float u);
	};
#endif
