#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/ibm.h"
#include "../include/mesh.h"
#include "../include/elemento.h"

int main(int argc, char *argv[])
{
	mesh mesh;
	float** nodos = mesh.darNodos();
	int** celdas = mesh.darCeldas();
	int nNodos = mesh.darNumeroNodos();
	int nCeldas = mesh.darNumeroCeldas();
	
	printf("Generando malla..\n");
	mesh.mesh_refine_tri4();
	mesh.proyectarEsfera(10);
	
	// Prueba de producto punto 
	float a[3]={1,2,3};
	float b[3]={1,2,3};
	float c[3]={2,4,6};
	cross(a,b,c[0], c[1], c[2]);
	printf("Resultado: %f %f %f\n", c[0], c[1], c[2]);

	// Prueba de funcion norma
	printf("Norma de a: %f\n", norm(a));


	// Prueba de fuerzas
	float referencia[3][3]={{0,0,0},{2.0,0,0},{1.0,1.0,0}};
	float deformado[3][3]={{0,0,0},{2.0,0,0},{1.0,1.0,0}};	
	float **fuerza;
	fuerza = fuerzas(referencia, deformado);
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			printf("%f ", fuerza[i][j]);
		}
		printf("\n");
	}	

	// Probar resultado Matriz dot Vector
	float v[3];
	MdotV(referencia, a, v[0], v[1], v[2]);	
	printf("Dot: %f %f %f\n", v[0], v[1], v[2]);

	// Prueba de fuerzas + rotacion 
	float ref[3][3]={{0,0,0},{2.0,0,0},{1.0,1.0,0}};
	float def[3][3]={{0,0,0},{2.0,0,0},{1.0,0.0,2.0}};	
	float **fuerza3;
	fuerza3 = rotacion(ref, def);
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			printf("%f ", fuerza3[i][j]);
		}
		printf("\n");
	}
	return 0;
}
