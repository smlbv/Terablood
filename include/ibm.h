#ifndef IBM_H
	#define IBM_H

	#include "mesh.h"
	#include "fluid.h"

	using namespace std;

	void interpolation(fluid fluido, mesh membrana, int x, int y, int z);
	void spread(fluid fluido, mesh membrana, int x, int y, int z);
	float dirac_2(float *x);
	float dirac_3(float *x);
	float dirac_4(float *x);

#endif
