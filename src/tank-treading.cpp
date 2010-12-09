// --------------------------------------------------------------
// 3D Sphere suspended in shear flow
// Terminado Octubre 09 - 2010
// --------------------------------------------------------------

typedef double Real;

// Librerias estandar
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Librerias del método
#include "../include/lbm.h"
#include "../include/ibm.h"
#include "../include/fem.h"

// Librerias para manejar objetos
#include "../include/fluid.h"
#include "../include/mesh.h"

using namespace std;

// ----------------------------------------------------------------
// 0. Conversion a cantidades adimensionales y parametros de flujo
// ----------------------------------------------------------------
float rho = 1.0;
float Re = 0.02;
float G = 0.01;
float tau = 1.0;
float nu = 1./6.;
float dt = 1.0;
int X = 25;
int Y = 25;
int Z = 25;
float H = Z;
float R = (X/10.0);
float gamma_dot = (Re*nu)/(R*R);
float u = (gamma_dot*H)/2.0;
float ks=(gamma_dot*nu*R)/G;
float k = gamma_dot/G;
float dx=H/X;
int STEPS = 200.0/k;
int VTK=1;
float fNeta[3] = {0,0,0};
float fMomento[3] = {0,0,0};

// -----------------------------------------------------------------------------
// Propiedades y variables para el fluido
// -----------------------------------------------------------------------------
fluid fluido;

// -----------------------------------------------------------------------------
// Propiedades y variables para la membrana
// -----------------------------------------------------------------------------
mesh membrana, referencia;

int main(int argc, char *argv[])
{
	// -----------------------------------------------------------------------//
	// 0. Inicializar
	// -----------------------------------------------------------------------//
	fluido.inicializar(X,Y,Z);
	fluido.setVelocidad(u);
	membrana.mesh_refine_tri4();
	membrana.mesh_refine_tri4();
	membrana.mesh_refine_tri4();
	membrana.mesh_refine_tri4();
	membrana.proyectarEsfera(R);
	membrana.moverCentro((X-1.0)/2.0, (Y-1.0)/2.0, (Z-1.0)/2.0);
	membrana.iniciarVelocidad();
	membrana.iniciarFuerzas();		
	
	referencia.mesh_refine_tri4();
	referencia.mesh_refine_tri4();
	referencia.mesh_refine_tri4();
	referencia.mesh_refine_tri4();
	referencia.proyectarEsfera(R);
	referencia.moverCentro((X-1.0)/2.0, (Y-1.0)/2.0, (Z-1.0)/2.0);
	referencia.iniciarVelocidad();
	referencia.iniciarFuerzas();		
	
	for(int ts = 0 ; ts < STEPS ; ts++)
	{
	// -----------------------------------------------------------------------//
	// 1. Interpolation
	// -----------------------------------------------------------------------//
	interpolation(fluido, membrana, X, Y, Z);

	// -----------------------------------------------------------------------//
	// 2. Encontrar nuevas posiciones de la membrana
	// -----------------------------------------------------------------------//
	//referencia.actualizarNodos(membrana.darNodos());
	membrana.moverNodos(dt, dx);

	// -----------------------------------------------------------------------//
	// 3. Calcular fuerzas en los nodos de la membrana
	// -----------------------------------------------------------------------//
	membrana.calcularFuerzas(referencia);

	// -----------------------------------------------------------------------//
	// 4. Propagar la densidad de fuerza hacia el fluido
	// -----------------------------------------------------------------------//
	spread(fluido, membrana, X, Y, Z);

	// -----------------------------------------------------------------------//
	// 5. Solucionar la dinámica del fluido
	// -----------------------------------------------------------------------//
	fluido.stream();
	fluido.collide();
	
	// -----------------------------------------------------------------------//
	// 6. Calcular propiedades macro del fluido
	// -----------------------------------------------------------------------//
	fluido.calcularMacro();

	// -----------------------------------------------------------------------//
	// 7. Calcular propiedades macro de la membrana
	// -----------------------------------------------------------------------//
	membrana.calcularCambioArea(referencia);
	
	// -----------------------------------------------------------------------//
	// 8. Visualización
	// -----------------------------------------------------------------------//
	if(ts%VTK==0)
	{
		fluido.guardar(ts);
		membrana.guardarVTU(ts);
		membrana.calcularVolumen();
		membrana.calcularFuerzaNeta(fNeta);
		membrana.calcularMomentoNeto(fMomento);
		membrana.calcularEnergia();
		printf("Equilibrio %f %f %f\n", fNeta[0], fNeta[1], fNeta[2]);
		printf("%d\n",ts);
	}
}//Ciclo principal
	return 0;
}
