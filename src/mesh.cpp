#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../include/mesh.h"
#include "../include/elemento.h"

// Constructor
mesh::mesh()
{
	float t = (1.+sqrt(5.))/2.;
	float tau = t/sqrt(1.+t*t);
        float one = 1./sqrt(1.+t*t); // Unit sphere

        // Twelve vertices of icosahedron on unit sphere
        // Creacion dinamica de arreglos al apuntador
        // Apuntadores a filas
        nNodos = 12;       
        vertex = new float*[12];
        
        // Apuntadores a columnas
        for(int i = 0; i<12 ; i++){
        	vertex[i] = new float[3];
        } 

	float nodos[12][3] ={{  tau,  one,    0},
                     { -tau,  one,    0 },
                     { -tau, -one,    0 },
                     {  tau, -one,    0 },
                     {  one,   0 ,  tau },
                     {  one,   0 , -tau },
                     { -one,   0 , -tau },
                     { -one,   0 ,  tau },
                     {   0 ,  tau,  one },
                     {   0 , -tau,  one },
                     {   0 , -tau, -one },
                     {   0 ,  tau, -one }};

	for(int i = 0; i<nNodos ; i++)
		for(int j = 0 ; j<3 ; j++)
		{
			vertex[i][j]=nodos[i][j];
		}

        		
        // Structure for unit icosahedron
        nCeldas = 20;
	faces = new int*[20];
	for(int i=0;i<nCeldas;i++)
	{
		faces[i]=new int[3];
	}
        int prueba[3][2]={(1,2,3),
        		  (4,5,6)};
        int f[20][3]= {{4,  7,  8},
                        {4,  9,  7},
                        {5, 11,  6},
                        {5,  6, 10},
                        {0,  3,  4},
                        {0,  5,  3},
                        {2,  1,  7},
                        {2,  6,  1},
                        {8, 11,  0},
                        {8,  1, 11},
                        {9,  3, 10},
                        {9, 10,  2},
                        {8,  0,  4},
                        {11,  5, 0},
                        {4,  3, 9},
                        {5, 10,  3},
                        {7,  1,  8},
                        {6, 11,  1},
                        {7, 9,  2},
                        {6,  2, 10}};
	for(int i = 0; i<nCeldas ; i++)
		for(int j = 0 ; j<3 ; j++)
		{
			faces[i][j]=f[i][j];
		}
}

//Destructor: Destruye los elementos de la malla nodos y celdas
mesh::~mesh()
{
	
}

// Trasladar en el espacio el centro de la esfera
void mesh::moverCentro(float x, float y, float z)
{
	cX = x;
	cY = y;
	cZ = z;
	for(int i=0;i<nNodos;i++)
	{
		vertex[i][0] = vertex[i][0]  + x;
		vertex[i][1] = vertex[i][1]  + y;
		vertex[i][2] = vertex[i][2]  + z;
	}
}

// Rotar la estructura completa dados tres angulos 
// theta, phi, psi
void mesh::rotarEstructura(float alpha, float phi, float theta)
{
	float Rx[3][3] = {{1,          0,          0},
			  {0, cos(alpha), -sin(alpha)},
			  {0, sin(alpha),  cos(alpha)}};	

	float Ry[3][3] = {{ cos(phi), 0, sin(phi)},
			  {        0, 1,        0},
			  {-sin(phi), 0, cos(phi)}};	
			  
        float Rz[3][3] = {{cos(theta), -sin(theta),0},
			  {sin(theta), cos(theta), 0},
			  {0, 0, 1}};		
	float Rxy[3][3];
	float Rxyz[3][3];
	
	// Realizar composición de rotaciones
	/* Realiza el producto de matrices y guarda el resultado en una tercera matriz*/
	for(int i=0;i<3;i++){ 
	      for(int j=0;j<3;j++){
		  Rxy[i][j]=0;
		  for(int k=0;k<3;k++){
		      Rxy[i][j]=Rxy[i][j]+(Rx[i][k]*Ry[k][j]);
		  }
	      }
	  }
	  
	  for(int i=0;i<3;i++){ 
	      for(int j=0;j<3;j++){
		  Rxyz[i][j]=0;
		  for(int k=0;k<3;k++){
		      Rxyz[i][j]=Rxyz[i][j]+(Rxy[i][k]*Rz[k][j]);
		  }
	      }
	  }
	  
	  // Rotar cada uno de los nodos por la matriz
	  for(int i = 0; i<nNodos ; i++)
	  {
	  	float x = Rxyz[0][0]*vertex[i][0] + Rxyz[0][1]*vertex[i][1] + Rxyz[0][2]*vertex[i][2];
	  	float y = Rxyz[1][0]*vertex[i][0] + Rxyz[1][1]*vertex[i][1] + Rxyz[1][2]*vertex[i][2];
	  	float z = Rxyz[2][0]*vertex[i][0] + Rxyz[2][1]*vertex[i][1] + Rxyz[2][2]*vertex[i][2];
	  	vertex[i][0] = x;
	  	vertex[i][1] = y;
	  	vertex[i][2] = z;
	  }
		
}


// Proyectar los nodos a una superficie RCB
void mesh::proyectarRBC(float r)
{

		float c0 = 0.207;
        float c1 = 2.000;
        float c2 = -1.123;
        float R  = r;
        
        for(int i = 0;i<nNodos;i++)
        {
		float X = vertex[i][0];
		float Y = vertex[i][1];
		float a = ((X*X)+(Y*Y))/(R*R);
		float b = (c0 + c1*a + c2*(a*a));
		float Z = (0.5*R*(sqrt(fabs(1.0-a))))*b;
		vertex[i][0]=X;
		vertex[i][1]=Y;
		
		// Calcular magnitud del vector
		float x, y, z;
		x = vertex[i][0]*vertex[i][0];
		y = vertex[i][1]*vertex[i][1];
		z = vertex[i][2]*vertex[i][2];
		float mag = sqrt(x+y+z);
		if(vertex[i][2]/mag<0){
			vertex[i][2]=-Z;// + (X*X)/16;
		}else{
			vertex[i][2]=Z;// + (X*X)/16;
		}
	}
}


// Proyectar los nodos a una superficie esferica 
void mesh::proyectarEsfera(float r)
{
	for(int i=0; i<nNodos;i++)
	{
		float x, y, z;
		x = vertex[i][0]*vertex[i][0];
		y = vertex[i][1]*vertex[i][1];
		z = vertex[i][2]*vertex[i][2];
		float mag= sqrt(x+y+z);
		vertex[i][0] = vertex[i][0]*r/mag;
		vertex[i][1] = vertex[i][1]*r/mag;
		vertex[i][2] = vertex[i][2]*r/mag;
	}
}


// Funcion que agrega un nodo a la lista de nodos de la malla
// Retorna la posicion del nuevo nodo
int mesh::agregarNodo(float x, float y, float z)
{
	float **nuevaLista = new float*[nNodos+1];
	for(int i=0; i<nNodos+1; i++)
	{
		nuevaLista[i]=new float[3];
	}
	for(int i=0; i<nNodos; i++)
	{
		nuevaLista[i][0]=vertex[i][0];
		nuevaLista[i][1]=vertex[i][1];
		nuevaLista[i][2]=vertex[i][2];
	}
	nuevaLista[nNodos][0]=x;
	nuevaLista[nNodos][1]=y;
	nuevaLista[nNodos][2]=z;
	vertex = nuevaLista;
	nNodos=nNodos+1;
	return nNodos-1;
}

// Funcion que determina si un nodo existe en la lista
// retorna -1 si no existe
// retorna  N si existe
int mesh::existeNodo(float x, float y, float z)
{
	int existe = -1;
	for(int i=0;i<nNodos;i++)
	{
		bool a = (vertex[i][0]==x);
		bool b = (vertex[i][1]==y);
		bool c = (vertex[i][2]==z);
		if(a &&  b && c)
		{
			existe = i;
		}
	}
	return existe;
}


//Función para refinamiento de malla
void mesh::mesh_refine_tri4()
{
	// Dividir cada cara en 4 triangulos
	// Make new midpoints
	// a = (A+B)/2
	// b = (B+C)/2
	// c = (C+A)/2
	//        B
	//       /\
	//      /  \
	//    a/____\b       Construct new triangles
	//    /\    /\       [A,a,c]
	//   /  \  /  \      [a,B,b]
	//  /____\/____\     [c,b,C]
	// A      c     C    [a,b,c]
	// -------------------------------------------
	
	printf("Refinando...\n");
	// Crear la estructura para las nuevas celdas
	int **nuevaLista = new int*[4*nCeldas];
	for(int i=0; i<nCeldas*4;i++)
	{
		nuevaLista[i] = new int[3];
	}
	
	// Dividir todas las caras 
	for(int f=0;f<nCeldas;f++)
	{
		int NA, NB, NC;
		// Encontrar los vertices correspondientes a cada cara
		NA = faces[f][0];
		NB = faces[f][1];
		NC = faces[f][2];
		
		// Encontrar las coordenadas de cada vertice
		float A[3], B[3], C[3];
		A[0]=vertex[NA][0];
		A[1]=vertex[NA][1];
		A[2]=vertex[NA][2];
		
		B[0]=vertex[NB][0];
		B[1]=vertex[NB][1];
		B[2]=vertex[NB][2];
		
		C[0]=vertex[NC][0];
		C[1]=vertex[NC][1];
		C[2]=vertex[NC][2];
		
		// Hallar los puntos medios de cada coordenada
		float a[3], b[3], c[3];
		a[0] = (A[0]+B[0])/2.;
		a[1] = (A[1]+B[1])/2.;
		a[2] = (A[2]+B[2])/2.;
		
		b[0] = (B[0]+C[0])/2.;
		b[1] = (B[1]+C[1])/2.;
		b[2] = (B[2]+C[2])/2.;
		
		c[0] = (C[0]+A[0])/2.;
		c[1] = (C[1]+A[1])/2.;
		c[2] = (C[2]+A[2])/2.;
		
		// Hallar el indice de cada nodo de las nuevas celdas
		int Na, Nb, Nc;
		Na = posicionNodo(a[0], a[1], a[2]);
		Nb = posicionNodo(b[0], b[1], b[2]);
		Nc = posicionNodo(c[0], c[1], c[2]);
		
		// Crear las nuevas caras
		int C1[3]={NA, Na, Nc}, C2[3]={Na, NB, Nb};
		int C3[3]={Nc, Nb, NC}, C4[3]={Na, Nb, Nc};
		
		// Agregarlas a la nueva estructura de celdas
		nuevaLista[(f+1)*4-4][0] = C1[0];
		nuevaLista[(f+1)*4-4][1] = C1[1];
		nuevaLista[(f+1)*4-4][2] = C1[2];
		
		nuevaLista[(f+1)*4-3][0] = C2[0];
		nuevaLista[(f+1)*4-3][1] = C2[1];
		nuevaLista[(f+1)*4-3][2] = C2[2];
		
		nuevaLista[(f+1)*4-2][0] = C3[0];
		nuevaLista[(f+1)*4-2][1] = C3[1];
		nuevaLista[(f+1)*4-2][2] = C3[2];
		
		nuevaLista[(f+1)*4-1][0] = C4[0];
		nuevaLista[(f+1)*4-1][1] = C4[1];
		nuevaLista[(f+1)*4-1][2] = C4[2];
	}
	faces = nuevaLista;
	nCeldas = 4*nCeldas;
}

// Funcion auxiliar para obtener la posicion de los vertices
int mesh::posicionNodo(float x, float y, float z)
{
	//Decidir si existe o no existe en la lista
	int pos = existeNodo(x,y,z);
	if(pos>-1)
	{
		return pos;
	}
	return agregarNodo(x,y,z);
}

/**

*/
int mesh::guardarVTU(int t)
{
	FILE *archivo;/*El manejador de archivo*/
	char ruta[80];
	char numero[4];
	strcpy(ruta, "temp/esfera-");
	sprintf(numero,"%d",t);
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
	fprintf(archivo, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(archivo, "POINTS %i double\n",nNodos);
	// Escribir coordenadas de cada nodo
	for(int i = 0; i<nNodos ; i++){
		for(int j = 0; j<3 ; j++){
			fprintf(archivo, "%f ", vertex[i][j]);
		}
		fprintf(archivo, "\n");}
	// Escribir celdas 
	fprintf(archivo, "CELLS %d %d\n",nCeldas,nCeldas*4);
	for(int i = 0; i<nCeldas ; i++){
		fprintf(archivo, "3 ");
		for(int j = 0; j<3 ; j++){
			fprintf(archivo, "%i ", faces[i][j]);
		}
		fprintf(archivo, "\n");}
	fprintf(archivo,"CELL_TYPES %d\n",nCeldas);
	for(int i = 0; i<nCeldas ; i++){
		fprintf(archivo, "5 \n");}

	// Escribir velocidad 
	fprintf(archivo,"POINT_DATA %d\n", nNodos);
	fprintf(archivo,"VECTORS Velocidad double\n");
	for(int i=0;i<nNodos;i++)
	{
		fprintf(archivo, "%f %f %f\n", velocidad[i][0], velocidad[i][1], velocidad[i][2]);
	}
	
	// Escribir fuerza
	fprintf(archivo,"VECTORS Fuerza double\n", nNodos);
	for(int i=0;i<nNodos;i++)
	{
		fprintf(archivo, "%f %f %f\n", fuerza[i][0], fuerza[i][1], fuerza[i][2]);
	}
	
	// Escribir cambio de area
	fprintf(archivo,"CELL_DATA %d \n", nCeldas);
	fprintf(archivo,"SCALARS dArea double\n");
	fprintf(archivo,"LOOKUP_TABLE default\n");
	for(int i=0;i<nCeldas;i++)
	{
		fprintf(archivo, "%f \n", area[i]);
	}

	fclose(archivo);/*Cerramos el archivo*/
	return 0;
	}
}

void mesh::iniciarVelocidad()
{
	velocidad = new float*[nNodos];
	velocidad2 = new float*[nNodos];		
	for(int i = 0 ; i<nNodos ; i++)
	{
		velocidad[i] = new float[3];
		velocidad2[i] = new float[3];
	}
}

void mesh::iniciarFuerzas()
{
	fuerza = new float*[nNodos];	
	for(int i = 0 ; i<nNodos ; i++)
	{
		fuerza[i] = new float[3];
	}
}


void mesh::setVelocidad(int n, float ux, float uy, float uz)
{
	velocidad2[n][0] = velocidad[n][0];
	velocidad2[n][1] = velocidad[n][1]; 
	velocidad2[n][2] = velocidad[n][2];  
	
	velocidad[n][0] = ux;
	velocidad[n][1] = uy;
	velocidad[n][2] = uz;
}



void mesh::darPosNodo(int n, float pos[3])
{
	pos[0] = vertex[n][0];
	pos[1] = vertex[n][1];
	pos[2] = vertex[n][2];
}

void mesh::darFuerzaNodo(int n, float f[3])
{
	f[0] = fuerza[n][0];
	f[1] = fuerza[n][1];
	f[2] = fuerza[n][2];
}


void mesh::moverNodos(float dt, float dx)
{
	for(int u=0;u<nNodos;u++)
	{
		vertex[u][0] += (3/2)*(velocidad[u][0]*dt) - (1./2.0)*(velocidad2[u][0]);
		vertex[u][1] += (3/2)*(velocidad[u][1]*dt) - (1./2.0)*(velocidad2[u][1]);
		vertex[u][2] += (3/2)*(velocidad[u][2]*dt) - (1./2.0)*(velocidad2[u][2]);
	}
}

void mesh::calcularFuerzas(mesh referencia)
{
	int a,b,c;
	float va[3], vb[3], vc[3], vA[3], vB[3], vC[3];
	// Crear estructuras para pasar al algoritmo de fuerzas
	float ref[3][3], def[3][3], fuerzas[3][3];
	
	// Obtener indices de los nodos de cada elemento
	for(int e=0;e<nCeldas;e++)
	{
		a = faces[e][0];
		b = faces[e][1];
		c = faces[e][2];
		
		// Obtener las posiciones de cada vertice no deformado
		referencia.darPosNodo(a,va);
		referencia.darPosNodo(b,vb);	
		referencia.darPosNodo(c,vc);		
	
		// Obtener las posiciones de cada vertice deformado
		darPosNodo(a,vA);
		darPosNodo(b,vB);	
		darPosNodo(c,vC);

		ref[0][0] = va[0];
		ref[0][1] = va[1];
		ref[0][2] = va[2];
	
		ref[1][0] = vb[0];
		ref[1][1] = vb[1];
		ref[1][2] = vb[2];
	
		ref[2][0] = vc[0];
		ref[2][1] = vc[1];
		ref[2][2] = vc[2];
	
		def[0][0] = vA[0];
		def[0][1] = vA[1];
		def[0][2] = vA[2];
	
		def[1][0] = vB[0];
		def[1][1] = vB[1];
		def[1][2] = vB[2];
	
		def[2][0] = vC[0];
		def[2][1] = vC[1];
		def[2][2] = vC[2];

		rotacion(ref, def, fuerzas);
		
		//Agregar cada fuerza a cada nodo
		fuerza[a][0] += fuerzas[0][0];
		fuerza[a][1] += fuerzas[0][1];
		fuerza[a][2] += fuerzas[0][2];
		
		fuerza[b][0] += fuerzas[1][0];
		fuerza[b][1] += fuerzas[1][1];
		fuerza[b][2] += fuerzas[1][2];
		
		fuerza[c][0] += fuerzas[2][0];
		fuerza[c][1] += fuerzas[2][1];
		fuerza[c][2] += fuerzas[2][2];
	}
}

/**
+	Retorna el volumen del tetrahedro formado por los vertices de cada elemento
*   y el centro de la esfera.
*
*/

float mesh::darVolumenElemento(int i)
{
	float a[3], b[3], c[3], d[3], v1[3], v2[3], v3[3], temp[3], V;
	int A, B, C;
	
	A=faces[i][0];
	B=faces[i][1];
	C=faces[i][2];
	
	a[0]=vertex[A][0];
	a[1]=vertex[A][1];
	a[2]=vertex[A][2];
	
	b[0]=vertex[B][0];
	b[1]=vertex[B][1];
	b[2]=vertex[B][2];
	
	c[0]=vertex[C][0];
	c[1]=vertex[C][1];
	c[2]=vertex[C][2];
	
	d[0] = cX;
	d[1] = cY;
	d[2] = cZ;
	
	v1[0] = a[0] - d[0];
	v1[1] = a[1] - d[1];
	v1[2] = a[2] - d[2];
	
	v2[0] = b[0] - d[0];
	v2[1] = b[1] - d[1];
	v2[2] = b[2] - d[2];
	
	v3[0] = c[0] - d[0];
	v3[1] = c[1] - d[1];
	v3[2] = c[2] - d[2];
	
	cross(v2, v3, temp[0], temp[1], temp[2]);
	V = (temp[0]*v1[0] + temp[1]*v1[1] + temp[2]*v1[2])/6.0;
	return fabs(V);
}

void mesh::calcularVolumen()

{
	volumen = 0.0;
	for(int i = 0; i<nCeldas ; i++)
	{
		volumen += darVolumenElemento(i);
	}
	printf("Volumen: %f\n", volumen);
}

float mesh::darAreaElemento(int i)
{
	float a[3], b[3], c[3], v1[3], v2[3], temp[3];
	int A, B, C;
	A=faces[i][0];
	B=faces[i][1];
	C=faces[i][2];
	
	a[0]=vertex[A][0];
	a[1]=vertex[A][1];
	a[2]=vertex[A][2];
	
	b[0]=vertex[B][0];
	b[1]=vertex[B][1];
	b[2]=vertex[B][2];
	
	c[0]=vertex[C][0];
	c[1]=vertex[C][1];
	c[2]=vertex[C][2];
	
	v1[0] = b[0] - a[0];
	v1[1] = b[1] - a[1];
	v1[2] = b[2] - a[2];
	
	v2[0] = c[0] - a[0];
	v2[1] = c[1] - a[1];
	v2[2] = c[2] - a[2];
	
	cross(v1, v2, temp[0], temp[1], temp[2]);
	return norm(temp)/2.0;
}



/**
*	Función para realizar cross product entre dos vectores
*   @param float a[3], Vector a
*	@param float b[3], Vector b
*	@param float &resultado[3], Paso por parametros
*/
void cross(float a[3], float b[3], float &x, float &y, float &z)
{
	x = a[1]*b[2] - a[2]*b[1];
	y = a[2]*b[0] - a[0]*b[2];
	z = a[0]*b[1] - a[1]*b[0];
}


/**
*	Función para calcular la norma de un vector
*   @param float a[3], Vector a
*	@return float norma, Norma del vector
*/
float norm(float a[3])
{
	return sqrt( (a[0]*a[0]) + (a[1]*a[1]) + (a[2]*a[2]));
}


void mesh::calcularCambioArea(mesh ref)
{
	if(area==NULL)
	{
		// Contruir la estructura
		area = new float[nCeldas];
	}
	for(int i = 0 ; i < nCeldas; i++)
	{
		area[i] = darAreaElemento(i)-ref.darAreaElemento(i);
	}
}


void mesh::actualizarNodos(float **nodos)
{
	float x, y, z;
	for(int u = 0 ; u< nNodos ; u++)
	{
		x = nodos[u][0];
		y = nodos[u][1];
		z = nodos[u][2];
				
		vertex[u][0] = x;
		vertex[u][1] = y;
		vertex[u][2] = z;
	}
}
