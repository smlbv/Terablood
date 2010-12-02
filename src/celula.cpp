#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/ibm.h"
#include "../include/mesh.h"

int main(int argc, char *argv[])
{
	mesh mesh;
	float** nodos = mesh.darNodos();
	int** celdas = mesh.darCeldas();
	int nNodos = mesh.darNumeroNodos();
	int nCeldas = mesh.darNumeroCeldas();
	
	printf("Generando malla..\n");
	mesh.mesh_refine_tri4();
	mesh.mesh_refine_tri4();
	mesh.mesh_refine_tri4();
	mesh.proyectarEsfera(10);
	mesh.proyectarRBC(10);
	mesh.guardarVTU(1);
	return 0;
}
