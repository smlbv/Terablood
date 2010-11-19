// --------------------------------------------------------------
// Flujo cortante en 3D usando Método de Lattice Boltzmann
// Terminado Octubre 09 - 2010
// --------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

using namespace std;

// global parameters
//const int SIZE = 50;
const int X = 50;
const int Y = 150;
const int Z = 50;

// single / double precision?
typedef double Real;

// different types of cells
const int FLUIDO  = 3, TOP = 1, BOTTOM  = 2, NOSLIP = 0;
// velocity for the moving wall
const Real U = 1.0, V = 0.00, W = 0.00;

// Declare some constants and globals ...
const Real   omega = 1.00;     // viscosity of the fluid, 0..2
const int    STEPS = 1;     // max no. of steps to simulate
const int    VTK = 10;  // write image every ? steps
int current = 0, other = 1; // which grid is the current and the old one?

// main arrays (note - for simplicity these are allocated here, for larger 
//              simulations, better use malloc and a manual offset calculation)
Real   cells[2][X][Y][Z][19]; // the two grids of LBM cells with 19 distribution functions
char   flags[X][Y][Z];        // flags for each cell either FLUID, NOSLIP or VELOCITY

// the weight for the equilibrium distribution according with Hecht
Real w[19] = {(2./36.),(2./36.),(2./36.),(2./36.),(2./36.),(2./36.),
	      (1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),
	      (1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),
	      (12./36.)};

// arrays with cartesian components of the 19 velocity vectors
int e_x[19] = {1,-1,0,0,0,0,1,1,1,1,-1,-1,-1,-1,0,0,0,0,0};
int e_y[19] = {0,0,1,-1,0,0,1,-1,0,0,1,-1,0,0,1,1,-1,-1,0};
int e_z[19] = {0,0,0,0,1,-1,0,0,1,-1,0,0,1,-1,1,-1,1,-1,0};
int dfInv[19] = {1,0,3,2,5,4,11,10,13,12,7,6,9,8,17,16,15,14,18};


//---------------------------------------------------------------------------//
// Declare auxiliary functions
//---------------------------------------------------------------------------//
void nodoSuperior(Real g[19], Real f[19]);
void nodoInferior(Real g[19], Real f[19]);
int guardarFluido(int s);

// run the solver
int main(int argc, char *argv[])
{
	printf("Init\n");
	
	// initialize grid
	for (int i=0;i<X;i++)
		for (int j=0;j<Y;j++)
			for (int k=0;k<Z;k++) {
				for (int l=0;l<19;l++)  {
					cells[0][i][j][k][l] = cells[1][i][j][k][l] = w[l];
				}
				if(k==0){flags[i][j][k]==BOTTOM;}
				else if(k==Z-1){flags[i][j][k]==TOP;}
				else{flags[i][j][k]==FLUIDO;}
			}
	
	printf("Starting simulation...\n");
	for(int s=0 ; s<STEPS ; s++)
	{
		printf("Step %d\n",s); 
	
		// Stream from the other to the current grid...
		int a=0,b=0,c=0;
		for (int i=0;i<X;i++)
			for (int j=0;j<Y;j++)
				for (int k=1;k<Z-1;k++)
					for (int l=0;l<19;l++) {
						int inv = dfInv[l];
						a = i + e_x[inv];
					        b = j + e_y[inv];
					        c = k + e_z[inv];

					        // Periodico en x
					        if(a<0){a=X-1;}
					        if(a>(X-1)){a=0;}
					    
					        // Periodico en y
					        if(b<0){b=Y-1;}
					        if(b>(Y-1)){b=0;}
					        if(flags[a][b][c]==TOP){
							// Bounce - back
							cells[current][i][j][k][l] = cells[other][i][j][k][inv];}
						else{
							// Streaming - normal
							cells[current][i][j][k][l] = cells[other][a][b][c][l];}
					}//Stream
		for (int i=0;i<X;i++)
			for (int j=0;j<Y;j++){
			nodoInferior(cells[current][i][j][0],cells[other][i][j][0]);
			nodoSuperior(cells[current][i][j][Z-1],cells[other][i][j][Z-1]);}
	
		// collision step
		for (int i=0;i<X;i++)
			for (int j=0;j<Y;j++)
				for (int k=0;k<Z;k++) {

					// standard collide without turbulence model
					Real rho = 0.0, u_x=0.0, u_y=0.0, u_z=0.0;
					// normal fluid cell
					for (int l=0;l<19;l++) {
						const Real fi = cells[current][i][j][k][l];
						rho += fi;
						u_x += e_x[l]*fi;
						u_y += e_y[l]*fi;
						u_z += e_z[l]*fi;
					}
					if(k==0){u_x=-U;}
					if(k==Z-1){u_x=U;}
					for (int l=0;l<19;l++) {
						const Real tmp = (e_x[l]*u_x + e_y[l]*u_y + e_z[l]*u_z);
						Real feq = w[l] * ( rho - 
							(3.0/2.0 * (u_x*u_x + u_y*u_y + u_z*u_z)) +
							(3.0 *     tmp) +
							(9.0/2.0 * tmp*tmp ) );
						cells[current][i][j][k][l] = 
							(1.0-omega) * cells[current][i][j][k][l] +
							omega * feq; 
					}
				} // ijk
		// We're done for one time step, switch the grid... 
		other = current;
		current = (current+1)%2;
		
		// should we write an image for this step?
		if(s%VTK==0){
			guardarFluido(s);
		}
	}// End simulation
	printf("LBM-simulation done!\n");
	return 0;
}


void nodoSuperior(Real g[19], Real f[19])
{
	// Calculate new distributions functions
	g=f;
	Real A=0.0, B=0.0, rho=0.0, Nx=0.0, Ny=0.0;
    	A=f[0]+f[1]+f[2]+f[3]+f[6]+f[10]+f[11]+f[7]+f[18];
    	B=f[4]+f[8]+f[12]+f[14]+f[16];
    	rho = (A+2*B)/(W+1);
    
   	Nx=(1./2.)*(f[0]+f[6]+f[7]-(f[1]+f[10]+f[11]))-(1./3.)*rho*U;
    	Ny=(1./2.)*(f[2]+f[6]+f[10]-(f[3]+f[7]+f[11]))-(1./3.)*rho*V;
    
    	g[5]=f[4]-(1./3.)*rho*W;
    	g[9]=f[12]+(rho/6)*(-W+U)-Nx;
    	g[13]=f[8]+(rho/6)*(-W-U)+Nx;
    	g[15]=f[16]+(rho/6)*(-W+V)-Ny;
    	g[17]=f[14]+(rho/6)*(-W-V)+Ny;
}

void nodoInferior(Real g[19], Real f[19])
{
	// Calculate new distributions functions
	g=f;
	Real A=0.0, B=0.0, rho=0.0, Nx=0.0, Ny=0.0;
	A=f[0]+f[1]+f[2]+f[3]+f[6]+f[7]+f[10]+f[11]+f[18];
	B=f[5]+f[9]+f[13]+f[15]+f[17];
    	rho = (A+2*B)/(1-W);
    
    	Nx=(1./2.)*(f[0]+f[6]+f[7]-(f[1]+f[10]+f[11]))-(1./3.)*rho*-U;
    	Ny=(1./2.)*(f[2]+f[6]+f[10]-(f[3]+f[7]+f[11]))-(1./3.)*rho*V;
    
    	g[4]=f[5]+(1./3.)*rho*W;
    	g[8]=f[13]+(rho/6)*(W-U)-Nx;
    	g[12]=f[9]+(rho/6)*(W+U)+Nx;
    	g[14]=f[17]+(rho/6)*(W+V)-Ny;
    	g[16]=f[15]+(rho/6)*(W-V)+Ny;
}


// Save fluid in structured grid format .vts
int guardarFluido(int s)
{
	FILE *archivo;/*El manejador de archivo*/
	char ruta[80];
	char numero[4];
	strcpy(ruta, "temp/fluido-");
	sprintf(numero,"%d",s);
	strcat(ruta,numero);
	strcat(ruta,".vtk");
        archivo=fopen(ruta, "w");
        if(archivo==NULL){/*Si no lo logramos abrir, salimos*/
		printf("No se puede guardar archivo");
		return 1;}
	else{
	// Escribir datos al archivo 
	// 1. Escribir cabecera.
	fprintf(archivo, "# vtk DataFile Version 3.0\n");
	fprintf(archivo, "vtk output\n");
	fprintf(archivo, "ASCII\n");
	fprintf(archivo, "DATASET STRUCTURED_GRID\n");
	fprintf(archivo, "DIMENSIONS %i %i %i\n", X, Y, Z);
	fprintf(archivo, "POINTS %i double\n",X*Y*Z);
	// 2. Escribir coordenadas de cada punto
	for(int i = 0 ;i<X;i++)
		for(int j = 0 ;j<Y;j++)
			for(int k = 0 ;k<Z;k++){
				fprintf(archivo,"%d %d %d\n", i,j,k);
			}			
	// 3. Escribir datos sobre puntos 
	
	// Densidad
	Real rho=0.0, u_x=0.0, u_y=0.0, u_z=0.0;
	fprintf(archivo, "POINT_DATA %d\n", X*Y*Z);
	fprintf(archivo, "SCALARS Densidad double\n");
	fprintf(archivo, "LOOKUP_TABLE default\n");
	for(int i = 0 ;i<X;i++)
		for(int j = 0 ;j<Y;j++)
			for(int k = 0 ;k<Z;k++){
				rho=0.0;			
				for(int l = 0 ;l<19;l++){
					const Real fi = cells[current][i][j][k][l];
					rho+= fi;
				}
				fprintf(archivo, "%f \n", rho);
			}
			
	// Velocidad
	fprintf(archivo, "\nVECTORS Velocidad double\n");
	for(int i = 0 ;i<X;i++)
		for(int j = 0 ;j<Y;j++)
			for(int k = 0 ;k<Z;k++){
				rho=0.0;
				u_x=0.0;
				u_y=0.0;
				u_z=0.0;
				for(int l = 0 ;l<19;l++){
					const Real fi = cells[current][i][j][k][l];
					rho+= fi;
					u_x+=fi*e_x[l];
					u_y+=fi*e_y[l];
					u_z+=fi*e_z[l];
				}
				u_x=u_x/rho;
				u_y=u_y/rho;
				u_z=u_z/rho;
				fprintf(archivo, "%f %f %f\n", u_x,u_y,u_z);
			}
        fclose(archivo);/*Cerramos el archivo*/
        return 0;}
}
