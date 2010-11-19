#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

/**
*	Función phi_2 para calcular función delta de Dirac con soporte de dos nodos
*	@param float r, Distancia entre nodos (Lagrangiana - Euleriana)
*   @return d, ponderación para distancia
*/
float phi_2(float r)
{
	if (0 <= abs(r) <= 1)
        return (1.0-abs(r));
    if ( 1 <= abs(r))
        return 0.0;
}

/**
*	Función phi_3 para calcular función delta de Dirac con soporte cuatro
*	@param float r, Distancia entre nodos (Lagrangiana - Euleriana)
*   @return float d, ponderación para distancia
*/
float phi_3(float r)
{
	if(0 <= abs(r) <= (1./2.))
        return ((1./3.)*(1+sqrt(1-3*r*r)));
    if((1./2.) <= abs(r) <= (3./2.))
        return ((1./6.)*(5-3+abs(r)-sqrt(-2+6*abs(r)-3*r*r)));
    if((3./2.) <= abs(r))
        return 0.0;
}


/**
*	Función phi_4 para calcular función delta de Dirac con soporte cuatro
*	@param float r, Distancia entre nodos (Lagrangiana - Euleriana)
*   @return d, ponderación para distancia
*/
float phi_4(float r)
{
	if(0 <= abs(r) <= 1)
        return ((1./8.)*(3-2*abs(r)+sqrt(1+4*abs(r)-4*r*r)));
    if(1 <= abs(r) <= 2)
        return ((1./8.)*(5-2*abs(r)-sqrt(-7+12*abs(r)-4*r*r)));
    if(2 <= abs(r))
        return 0.0;
}


/**
*	Función dirac_4 para calcular función delta de Dirac con soporte cuatro
*	@param x, Vector distancia entre nodos (Lagrangiana - Euleriana)
*   @return d, ponderación por distancia
*/
float dirac_2(float *x)
{
    float d = phi_2(x[0])*phi_2(x[1])*phi_2(x[2]);
    return d;
}

/**
*	Función dirac_3 para calcular función delta de Dirac con soporte cuatro
*	@param x, Vector distancia entre nodos (Lagrangiana - Euleriana)
*   @return d, ponderación por distancia
*/
float dirac_3(float *x)
{
    float d = phi_3(x[0])*phi_3(x[1])*phi_3(x[2]);
    return d;
}

/**
*	Función dirac_4 para calcular función delta de Dirac con soporte cuatro
*	@param x, Vector distancia entre nodos (Lagrangiana - Euleriana)
*   @return d, ponderación por distancia
*/
float dirac_4(float *x)
{
    float d = phi_4(x[0])*phi_4(x[1])*phi_4(x[2]);
    return d;
}

/**
* TODO: Implementar las funcion spread
*/

void spread()
{
	printf("\nSpread function\n");
	printf("Segunda compilacion\n");
} 

/**
* TODO: Implementar las funcion interpolation
*/
void interpolation()
{


}
