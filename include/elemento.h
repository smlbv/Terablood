#ifndef ELEMENTO_H
	#define ELEMENTO_H
	using namespace std;
	void fuerzas(float referencia[3][3], float deformado[3][3], float fuerzas[3][3]);
	void rotacion(float referencia[3][3], float deformado[3][3], float fuerzas[3][3]);
	void cross(float a[3], float b[3], float &x, float &y, float &z);
	float norm(float a[3]);
	void MdotV(float M[3][3], float v[3], float &x, float &y, float &z);
#endif
