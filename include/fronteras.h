/**
* TODO: Implementar cada uno de los tipos de condiciones de frontera 
* tanto de velocidad como de presión
*/

typedef double Real;

#ifndef FRONTERAS_H
	#define FRONTERAS_H
	using namespace std;

	// Condiciones de frontera para veocidad
	void velNodoIzquierdo(Real g[19], Real f[19], Real U, Real V, Real W);
	void velNodoDerecho(Real g[19], Real f[19], Real U, Real V, Real W);
	void velNodoSuperior(Real g[19], Real f[19], Real U, Real V, Real W);
	void velNodoInferior(Real g[19], Real f[19], Real U, Real V, Real W);

	/**
	* TODO: Implementar condiciones de frontera ṕara presion
	*/

#endif
