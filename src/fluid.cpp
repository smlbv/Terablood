// Librerias estandar
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../include/fluid.h"
#include "../include/fronteras.h"

void fluid::inicializar(int x, int y, int z)
{
	X = x;
	Y = y;
	Z = z;
	
	// Construye las estructuras de información macroscopicas
	vel = new Real***[x];
	fuerza = new Real***[x];
	rho = new Real**[x];
	flags = new Real**[x];
	for(int i=0;i<x;i++)
	{
		vel[i] = new Real**[y];
		fuerza[i] = new Real**[y];
		rho[i] = new Real*[y];
		flags[i] = new Real*[y];
		for(int j = 0; j<y;j++)
		{
			vel[i][j] = new Real*[z];
			fuerza[i][j] = new Real*[z];
			rho[i][j] = new Real[z];
			flags[i][j] = new Real[z];
			for(int k=0; k<z ; k++)
			{
				vel[i][j][k] = new Real[3];
				fuerza[i][j][k] = new Real[3];
				rho[i][j][k] = 0;
				flags[i][j][k]=0;
			}
		}
	}
	
	// Construye las estructuras de información microscopicas
	cells = new Real****[2];
	for(int s=0;s<2;s++)
	{
		cells[s]=new Real***[x];
		for(int i=0;i<x;i++)
		{
			cells[s][i]=new Real**[y];
			for(int j = 0; j<y;j++)
			{
				cells[s][i][j]=new Real*[z];
				for(int k=0; k<z ; k++)
				{
					cells[s][i][j][k]=new Real[19];
					for(int l=0;l<19;l++)
					{
						cells[s][i][j][k][l] = w[l];
					}
				}
			}
		}
	}
	
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
}


void fluid::stream()
{

		int a=0,b=0,c=0;
		omp_set_num_threads(2);
#		pragma omp parallel
		{
#		pragma omp for  schedule(static)
		for (int i=0;i<X;i++)
			for (int j=0;j<Y;j++)
				for (int k=1;k<Z-1;k++)
					for (int l=0;l<19;l++) {
						int inv = dfInv[l];
						int a = i + e_x[inv];
					    int b = j + e_y[inv];
					    int c = k + e_z[inv];

					        // Periodico en 
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
			velNodoInferior(cells[current][i][j][0],cells[other][i][j][0], U, V, W);
			velNodoSuperior(cells[current][i][j][Z-1],cells[other][i][j][Z-1], U, V, W);
			}
		}
}


void fluid::collide()
{
	omp_set_num_threads(2);
#		pragma omp parallel
		{

#		pragma omp for  schedule(static)
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

					// Fronteras
					if(k==0){u_x=-U;}
					if(k==Z-1){u_x=U;}
					
					for (int l=0;l<19;l++) {
						const Real tmp = (e_x[l]*u_x + e_y[l]*u_y + e_z[l]*u_z);
						// Función de equilibrio
						Real feq = w[l] * rho * ( 1.0 - 
							(3.0/2.0 * (u_x*u_x + u_y*u_y + u_z*u_z)) +
							(3.0 *     tmp) +
							(9.0/2.0 * tmp*tmp ) );
//						Real feq = w[l]*rho*(1.0+3.0*tmp);
						// Fuerza por cada dirección i
						Real v1[3]={0.0,0.0,0.0};
						v1[0]=(e_x[l]-u_x)/(cs*cs);
						v1[1]=(e_y[l]-u_y)/(cs*cs);
						v1[2]=(e_z[l]-u_z)/(cs*cs);
						
						v1[0]+=(tmp*e_x[l])/(cs*cs*cs*cs);
						v1[1]+=(tmp*e_y[l])/(cs*cs*cs*cs);
						v1[2]+=(tmp*e_z[l])/(cs*cs*cs*cs);
						
						Real Fi=0.0;
						Fi = (v1[0]*fuerza[i][j][k][0] + v1[1]*fuerza[i][j][k][1] + v1[2]*fuerza[i][j][k][2]);
						Fi = (1.0-(1.0/(2.0*omega)))*w[l]*Fi;
						
						cells[current][i][j][k][l] = cells[current][i][j][k][l] - omega*(cells[current][i][j][k][l] - feq) + Fi;
					}
					fuerza[i][j][k][0]=0.0;
					fuerza[i][j][k][1]=0.0;
					fuerza[i][j][k][2]=0.0;					
				} // ijk
			}
		// We're done for one time step, switch the grid... 
		other = current;
		current = (current+1)%2;
}

void fluid::calcularMacro()
{
	for(int i=0;i<X;i++)
		for(int j=0;j<Y;j++)
			for(int k=0;k<Z;k++){
				float rho=0.0, ux=0.0, uy=0.0, uz=0.0;
				for(int l=0;l<19;l++)
				{
					const Real fi = cells[current][i][j][k][l];
					rho+=fi;
					ux+=fi*e_x[l];
					uy+=fi*e_y[l];
					uz+=fi*e_z[l];
				}
				vel[i][j][k][0] = (ux+(fuerza[i][j][k][0]/2.0))/rho;
				vel[i][j][k][1] = (uy+(fuerza[i][j][k][1]/2.0))/rho;
				vel[i][j][k][2] = (uz+(fuerza[i][j][k][2]/2.0))/rho;
				}
}

// Save fluid in structured grid format .vts
int fluid::guardar(int s)
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
		// Escribir vectores fuerza en el algoritmo
		fprintf(archivo, "VECTORS fuerza double\n");
		for(int i = 0 ;i<X;i++)
			for(int j = 0 ;j<Y;j++)
				for(int k = 0 ;k<Z;k++)
				{
					fprintf(archivo, "%f %f %f\n", fuerza[i][j][k][0], fuerza[i][j][k][0], fuerza[i][j][k][0]);
				}
        fclose(archivo);/*Cerramos el archivo*/
        return 0;
	}
}


float fluid::darVelocidad(int x, int y, int z, int f)
{
	return vel[x][y][z][f];
}

void fluid::setFuerza(int x, int y, int z, float f[3])
{
	fuerza[x][y][z][0] = f[0];
	fuerza[x][y][z][1] = f[1];
	fuerza[x][y][z][2] = f[2];
}

void fluid::addFuerza(int x, int y, int z, float f[3])
{
	fuerza[x][y][z][0] += f[0];
	fuerza[x][y][z][1] += f[1];
	fuerza[x][y][z][2] += f[2];
}


void fluid::setVelocidad(float u)
{
	U = u;
}
