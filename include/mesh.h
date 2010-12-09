#ifndef MESH_H
	#define MESH_H
	using namespace std;

	/**
	* @mesh.h
	* @author  Oscar Castillo <ol.castillo28@uniandes.edu.co>
	* @version 1.0
	*
	* @section LICENSE
	*
	* This program is free software; you can redistribute it and/or
	* modify it under the terms of the GNU General Public License as
	* published by the Free Software Foundation; either version 2 of
	* the License, or (at your option) any later version.
	* 
	* This program is distributed in the hope that it will be useful, but
	* WITHOUT ANY WARRANTY; without even the implied warranty of
	* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
	* General Public License for more details at
	* http://www.gnu.org/copyleft/gpl.html
	*
	* @section DESCRIPTION
	*
	* The mesh class represents a collection of point sand faces to reproduce
	* a 3D object representation.
	*/

	class mesh{

		private:
		
		float cX, cY, cZ;
		float **vertex;
		float **velocidad;
		float **velocidad2;
		float **fuerza;
		float *area;
		int **faces;
		int nNodos;	//Numero de nodos
		int nCeldas;	//Número de celdas
		float volumen; // Volumen encerrado por la capsula
	
		public:
		
		float** darNodos(){return vertex;}
		int** darCeldas(){return faces;}
		int darNumeroNodos(){return nNodos;}
		int darNumeroCeldas(){return nCeldas;}
		int posicionNodo(float x, float y, float z);
		int guardarVTU(int t);
		int agregarNodo(float x, float y, float z);
		int agregarCelda(int a, int b, int c);
		int existeNodo(float x, float y, float z);
		void mesh_refine_tri4();
		void proyectarEsfera(float r);
		void proyectarRBC(float r);
		void moverCentro(float x, float y, float z);
		void rotarEstructura(float alpha, float phi, float theta);
		void iniciarVelocidad();
		void iniciarFuerzas();
		void iniciarArea();
		void darPosNodo(int n, float pos[3]);
		void darFuerzaNodo(int n, float f[3]);
		void setVelocidad(int n, float ux, float uy, float uz);
		void moverNodos(float dt, float dx);
		void calcularFuerzas(mesh referencia);
		void actualizarNodos(float **);
		void calcularCambioArea(mesh ref);
		float darAreaElemento(int i);
		float darVolumenElemento(int i);
		void calcularVolumen();
		// Constructor
		mesh();
		// Destructor
		~mesh();
	};

#endif
