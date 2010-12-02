#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

using namespace std;

/**
*	Función para calcular las fuerzas en un elemento triangular
*   @param float referencia[3][3], Vertices iniciales del elemento
*	@param float deformado[3][3], Vertices del elemento deformado en el plano
*	@return float** fuerzas, Apuntador hacia las fuerzas calculadas en cada nodo
*/

void fuerzas(float referencia[3][3], float deformado[3][3], float fuerzas[3][3])
{
		// Coordenadas de los tres vertices iniciales			
		float xi=referencia[0][0];
		float yi=referencia[0][1];
		float zi=referencia[0][2];

		float xj=referencia[1][0];
		float yj=referencia[1][1];
		float zj=referencia[1][2];

		float xk=referencia[2][0];
		float yk=referencia[2][1];
		float zk=referencia[2][2];

		// Coordenadas de los tres vertices del elemento deformado
		float Xi=deformado[0][0];
		float Yi=deformado[0][1];
		float Zi=deformado[0][2];
		
		float Xj=deformado[1][0];
		float Yj=deformado[1][1];
		float Zj=deformado[1][2];

		float Xk=deformado[2][0];
		float Yk=deformado[2][1];
		float Zk=deformado[2][2];

		// Desplazamientos de cada nodo
		float ui=Xi-xi;
		float vi=Yi-yi;
		float wi=Zi-zi;

		float uj=Xj-xj;
		float vj=Yj-yj;
		float wj=Zj-zj;

		float uk=Xk-xk;
		float vk=Yk-yk;
		float wk=Zk-zk;
		
		// Coeficientes de funciones de forma
		float ai = yj-yk;
		float bi = xk-xj;
		float ci = xj*yk - xk*yj;
		float Li = ai*xi + bi*yi + ci;
		
		float aj = yk - yi;
		float bj = xi - xk;
		float cj = xk*yi - xi*yk;
		float Lj = aj*xj + bj*yj + cj;
				
		float ak = yi - yj;
		float bk = xj - xi;
		float ck = xi*yj - xj*yi;
		float Lk = ak*xk + bk*yk + ck;

		// Derivadas parciales para calcular el vector [G] revision Sep/14/2010
		float dudx = ui*ai/Li + uj*aj/Lj + uk*ak/Lk;
		float dudy = ui*bi/Li + uj*bj/Lj + uk*bk/Lk;
		float dvdx = vi*ai/Li + vj*aj/Lj + vk*ak/Lk;
		float dvdy = vi*bi/Li + vj*bj/Lj + vk*bk/Lk;
		
		// Componentes del vector [G] revision Sep/14/2010
		float g11 = (1.0+dudx)*(1.0+dudx) + (dvdx)*(dvdx);
		float g12 = (1.0+dudx)*(dudy) + (1.0+dvdy)*(dvdx);
		float g21 = g12;
		float g22 = (1.0+dvdy)*(1.0+dvdy) + (dudy)*(dudy);
		float G[2][2] = {{g11,g12},{g21,g22}};

		// Calculo de lambda1 y lambda2 revision Sep/14/2010
		float l1 = sqrt((g11 + g22 + sqrt((g11-g22)*(g11-g22) + 4.*g12*g12))/2.);
		float l2 = sqrt((g11 + g22 - sqrt((g11-g22)*(g11-g22) + 4.*g12*g12))/2.);

		// Derivadas de la funcion Strain Energy respecto a lambda 1 y lambda 2 
		// Modelo de energia Skalak 1973 revision Sep/14/2010
		float I1 = (l1*l1) + (l2*l2) - 2.0;
		float I2 = (l1*l1)*(l2*l2)-1.0;
		float dI1dl1 = 2.0*l1; 
		float dI1dl2 = 2.0*l2;
		float dI2dl1 = 2.0*l1*(l2*l2);
		float dI2dl2 = 2.0*l2*(l1*l1);
		
		//    dwdl1 = (B/4)*(I1*dI1dl1 + dI1dl1 - dI2dl1) + (C/4)*(I2)*(dI2dl1)
		//    dwdl2 = (B/4)*(I1*dI1dl2 + dI1dl2 - dI2dl2) + (C/4)*(I2)*(dI2dl2)
		float ks = 0.01587;
		float dwdl1 = (ks/12.)*(2.*I1*dI1dl1 + 2.0*dI1dl1 -2.0*dI2dl1) + (ks/6.0)*I2*dI2dl1;
		float dwdl2 = (ks/12.)*(2.*I1*dI1dl2 + 2.0*dI1dl2 -2.0*dI2dl2) + (ks/6.0)*I2*dI2dl2;

		// Calculo de diferenciales sobre l1 y l2 respecto a desplazamientos de nodos
		// 1. Derivadas de [G] respecto a desplazamiento de nodos revision Sep/14/2010
		float dg11dui = 2.0*(1.0+ dudx)*(ai/Li);
		float dg11duj = 2.0*(1.0+ dudx)*(aj/Lj);
		float dg11duk = 2.0*(1.0+ dudx)*(ak/Lk);
		float dg11dvi = 2.0*dvdx*(ai/Li);
		float dg11dvj = 2.0*dvdx*(aj/Lj);
		float dg11dvk = 2.0*dvdx*(ak/Lk);
		
		float dg12dui = (1.0+dudx)*(bi/Li) + (ai/Li)*(dudy);
		float dg12duj = (1.0+dudx)*(bj/Lj) + (aj/Lj)*(dudy);
		float dg12duk = (1.0+dudx)*(bk/Lk) + (ak/Lk)*(dudy);
		float dg12dvi = (1.0+dvdy)*(ai/Li) + (bi/Li)*(dvdx);
		float dg12dvj = (1.0+dvdy)*(aj/Lj) + (bj/Lj)*(dvdx);
		float dg12dvk = (1.0+dvdy)*(ak/Lk) + (bk/Lk)*(dvdx);
		
		float dg22dui = 2.0*dudy*(bi/Li);
		float dg22duj = 2.0*dudy*(bj/Lj);
		float dg22duk = 2.0*dudy*(bk/Lk);
		
		float dg22dvi = 2.0*(1.0+dvdy)*(bi/Li);
		float dg22dvj = 2.0*(1.0+dvdy)*(bj/Lj);
		float dg22dvk = 2.0*(1.0+dvdy)*(bk/Lk);

		// 2. Calculo de las derivadas de lambda 1 y lambda 2 respecto desplazamientos
		// nodales revision Sep/14/2010
		// Formulacion Rolling John Hopkins
		float t0 = sqrt((g11-g22)*(g11-g22) + 4.*g12*g12);
		
		// Derivadas de lambda1  y lambda 2
		
		float dt0dui;
		if(t0>10e-3){
			dt0dui = ((((g11-g22)*(dg11dui-dg22dui))+(4.*g12*dg12dui))/(t0));}
		else{
			dt0dui = 0.0;}
		float dl1dui = ((sqrt(0.5)*0.5)/(sqrt(g11+g22+t0)))*(dg11dui + dg22dui + dt0dui);
		float dl2dui = ((sqrt(0.5)*0.5)/(sqrt(g11+g22-t0)))*(dg11dui + dg22dui - dt0dui);
		
		float dt0duj;
		if(t0>10e-3){
			dt0duj = ((((g11-g22)*(dg11duj-dg22duj))+(4.*g12*dg12duj))/(t0));}
		else{
			dt0duj = 0.0;}
		float dl1duj = ((sqrt(0.5)*0.5)/(sqrt(g11+g22+t0)))*(dg11duj + dg22duj + dt0duj);
		float dl2duj = ((sqrt(0.5)*0.5)/(sqrt(g11+g22-t0)))*(dg11duj + dg22duj - dt0duj);
		
		float dt0duk;
		if(t0>10e-3){
			dt0duk = ((((g11-g22)*(dg11duk-dg22duk))+(4.*g12*dg12duk))/(t0));}
		else{
			dt0duk = 0.0;}
		float dl1duk = ((sqrt(0.5)*0.5)/(sqrt(g11+g22+t0)))*(dg11duk + dg22duk + dt0duk);
		float dl2duk = ((sqrt(0.5)*0.5)/(sqrt(g11+g22-t0)))*(dg11duk + dg22duk - dt0duk);
		
		float dt0dvi;
		if(t0>10e-3){
			dt0dvi = ((((g11-g22)*(dg11dvi-dg22dvi))+(4.*g12*dg12dvi))/(t0));}
		else{
			dt0dvi = 0.0;}
		float dl1dvi = ((sqrt(0.5)*0.5)/(sqrt(g11+g22+t0)))*(dg11dvi + dg22dvi + dt0dvi);
		float dl2dvi = ((sqrt(0.5)*0.5)/(sqrt(g11+g22-t0)))*(dg11dvi + dg22dvi - dt0dvi);
		
		float dt0dvj;
		if(t0>10e-3){
			dt0dvj = ((((g11-g22)*(dg11dvj-dg22dvj))+(4.*g12*dg12dvj))/(t0));}
		else{
			dt0dvj = 0.0;}
		float dl1dvj = ((sqrt(0.5)*0.5)/(sqrt(g11+g22+t0)))*(dg11dvj + dg22dvj + dt0dvj);
		float dl2dvj = ((sqrt(0.5)*0.5)/(sqrt(g11+g22-t0)))*(dg11dvj + dg22dvj - dt0dvj);
		
		float dt0dvk;
		if(t0>10e-3){
			dt0dvk = ((((g11-g22)*(dg11dvk-dg22dvk))+(4.*g12*dg12dvk))/(t0));}
		else{
			dt0dvk = 0.0;}
		float dl1dvk = ((sqrt(0.5)*0.5)/(sqrt(g11+g22+t0)))*(dg11dvk + dg22dvk + dt0dvk);
		float dl2dvk = ((sqrt(0.5)*0.5)/(sqrt(g11+g22-t0)))*(dg11dvk + dg22dvk - dt0dvk);				
		
		// 3. Calculo de las derivadas de w respecto a los desplazamientos nodales
		// revision Sep/14/2010
		float dwdui = dwdl1*dl1dui + dwdl2*dl2dui;
		float dwdvi = dwdl1*dl1dvi + dwdl2*dl2dvi;
		float dwduj = dwdl1*dl1duj + dwdl2*dl2duj;
		float dwdvj = dwdl1*dl1dvj + dwdl2*dl2dvj;
		float dwduk = dwdl1*dl1duk + dwdl2*dl2duk;
		float dwdvk = dwdl1*dl1dvk + dwdl2*dl2dvk;

		// 4. Volumen del elemento revision Sep/14/2010
		float a0 = ((xj-xi)*(yk-yi) - (xk-xi)*(yj-yi))/2.0;
		float espesor = 3.5*0.0025;
		
		// 5. Calculo de las componentes de fuerza revision Sep/14/2010
		float fxi = dwdui*a0*espesor;
		float fyi = dwdvi*a0*espesor;
		float fzi = 0.0;
		float fxj = dwduj*a0*espesor;
		float fyj = dwdvj*a0*espesor;
		float fzj = 0.0;
		float fxk = dwduk*a0*espesor;
		float fyk = dwdvk*a0*espesor;
		float fzk = 0.0;
		
		fuerzas[0][0] = fxi;
		fuerzas[0][1] = fyi;
		fuerzas[0][2] = fzi;
		fuerzas[1][0] = fxj;
		fuerzas[1][1] = fyj;
		fuerzas[1][2] = fzj;
		fuerzas[2][0] = fxk;
		fuerzas[2][1] = fyk;
		fuerzas[2][2] = fzk;
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



/**
*	Función para calcular el producto Matriz dot Vector
*   @param M[3][3], Matriz de tamaño 3x3
*   @param V[3], Vector de tamaño 3
*   @param V[3], Vector de tamaño 3
*/
void MdotV(float M[3][3], float V[3], float &x, float &y, float &z)
{
	x = M[0][0]*V[0] + M[0][1]*V[1] + M[0][2]*V[2];
	y = M[1][0]*V[0] + M[1][1]*V[1] + M[1][2]*V[2];
	z = M[2][0]*V[0] + M[2][1]*V[1] + M[2][2]*V[2];
}

/**
*	Función para llevar los elementos de referencia y deformado a una misma base
*   @param float referencia[3][3], Vertices iniciales del elemento
*	@param float deformado[3][3], Vertices del elemento deformado en 3D
*	@return float** vertices, Apuntador hacia los vertices en un plano común
*/

void rotacion(float referencia[3][3], float deformado[3][3], float nfuerzas[3][3])
{

	// Coordenadas iniciales de los tres nodos respecto al sistema 
    // global de coordenadas i1, i2, i3
    float xi = referencia[0][0];
    float yi = referencia[0][1];
    float zi = referencia[0][2];
    
    float xj = referencia[1][0];
    float yj = referencia[1][1];
    float zj = referencia[1][2];
    
    float xk = referencia[2][0];
    float yk = referencia[2][1];
    float zk = referencia[2][2];
    
    // Coordenadas del elemento deformado de los tres nodos 
    // respecto al sistema global de coordenadas i1, i2, i3
    float Xi = deformado[0][0];
    float Yi = deformado[0][1];
    float Zi = deformado[0][2];
    
    float Xj = deformado[1][0];
    float Yj = deformado[1][1];
    float Zj = deformado[1][2];
    
    float Xk = deformado[2][0];
    float Yk = deformado[2][1];
    float Zk = deformado[2][2];
    
    // Unit vectors Local undeformed local coordinate axis 
    float m1 = sqrt((xj-xi)*(xj-xi) + (yj-yi)*(yj-yi) + (zj-zi)*(zj-zi));
	float ve1[3] = {(xj-xi)/m1, (yj-yi)/m1, (zj-zi)/m1};
    float m2 = sqrt((xk-xi)*(xk-xi) + (yk-yi)*(yk-yi) + (zk-zi)*(zk-zi));
	float ve4[3] = {(xk-xi)/m2, (yk-yi)/m2, (zk-zi)/m2};
    float ve3[3];
	cross(ve1, ve4, ve3[0], ve3[1], ve3[2]);
	float m3 = norm(ve3);

	ve3[0] = ve3[0]/m3;
	ve3[1] = ve3[1]/m3;
	ve3[2] = ve3[2]/m3;

    float ve2[3];
	cross(ve3,ve1,ve2[0], ve2[1], ve2[2]);
    
	// Formar la base e
	float e[3][3];
	e[0][0] = ve1[0];
	e[0][1] = ve1[1];
	e[0][2] = ve1[2];
	
	e[1][0] = ve2[0];
	e[1][1] = ve2[1];
	e[1][2] = ve2[2];

	e[2][0] = ve3[0];
	e[2][1] = ve3[1];
	e[2][2] = ve3[2];

	// Matriz de rotacion para la configuracion NO deformada 
	// [r] hacia las coordenadas globales
    float d1 = (xj-xi)/m1;
    float d2 = (yj-yi)/m1;
    float d3 = (zj-zi)/m1;
    float e1 = (xk-xi)/m2;
    float e2 = (yk-yi)/m2;
    float e3 = (zk-zi)/m2;
    float f1 = (d2*e3-d3*e2)/m3;
    float f2 = (d3*e1-d1*e3)/m3;
    float f3 = (d1*e2-d2*e1)/m3;
    float g1 = f2*d3 - f3*d2;
    float g2 = f3*d1 - f1*d3;
    float g3 = f1*d2 - f2*d1;
    float r[3][3] = {{d1, d2, d3},{g1, g2, g3},{f1, f2, f3}};

    // Unit vectors Local deformed coordinate axis
    float M1 = sqrt((Xj-Xi)*(Xj-Xi) + (Yj-Yi)*(Yj-Yi) + (Zj-Zi)*(Zj-Zi));
    float E1[3] = {(Xj-Xi)/M1, (Yj-Yi)/M1, (Zj-Zi)/M1};
    float M2 = sqrt((Xk-Xi)*(Xk-Xi) + (Yk-Yi)*(Yk-Yi) + (Zk-Zi)*(Zk-Zi));
    float E4[4] = {(Xk-Xi)/M2, (Yk-Yi)/M2, (Zk-Zi)/M2};
	float E3[3];    
	cross(E1,E4,E3[0], E3[1], E3[2]);
    float M3 = norm(E3);
    E3[0] = E3[0]/M3;
    E3[1] = E3[1]/M3;
    E3[2] = E3[2]/M3;
    float E2[3];
	cross(E3,E1,E2[0],E2[1],E2[2]);
    
	// Formar la base E
	float E[3][3];
	E[0][0] = E1[0];
	E[0][1] = E1[1];
	E[0][2] = E1[2];
	
	E[1][0] = E2[0];
	E[1][1] = E2[1];
	E[1][2] = E2[2];

	E[2][0] = E3[0];
	E[2][1] = E3[1];
	E[2][2] = E3[2];

	// Matriz de rotacion para la configuracion deformada 
	// [R] hacia la base global i
    d1 = (Xj-Xi)/M1;
    d2 = (Yj-Yi)/M1;
    d3 = (Zj-Zi)/M1;
    e1 = (Xk-Xi)/M2;
    e2 = (Yk-Yi)/M2;
    e3 = (Zk-Zi)/M2;
    f1 = (d2*e3-d3*e2)/M3;
    f2 = (d3*e1-d1*e3)/M3;
    f3 = (d1*e2-d2*e1)/M3;
    g1 = f2*d3 - f3*d2;
    g2 = f3*d1 - f1*d3;
    g3 = f1*d2 - f2*d1;
    float R[3][3] = {{d1, d2, d3},{g1, g2, g3},{f1, f2, f3}};

	// Vectores posicion de los nodos en coordenadas locales estado deformado
    // tiene origen en el nodo i deformado 
    float Xil = 0.0;
    float Yil = 0.0;
    float Zil = 0.0;
    float Xjl = Xj - Xi;
    float Yjl = Yj - Yi;
    float Zjl = Zj - Zi;
    float Xkl = Xk - Xi;
    float Ykl = Yk - Yi;
    float Zkl = Zk - Zi;


	// Transformar cada coordenada local del elemento deformado
    // a la base comun i,j,k 
	float Pi[3];
	float v1[3] = {Xil, Yil, Zil};
	MdotV(R, v1, Pi[0], Pi[1],Pi[2]);

	float Pj[3];    
	float v2[3]={Xjl, Yjl, Zjl};
	MdotV(R,v2,Pj[0],Pj[1],Pj[2]);

	float Pk[3];    
	float v3[3]={Xkl, Ykl, Zkl};
	MdotV(R,v3,Pk[0],Pk[1],Pk[2]);

	// Vectores posicion de los nodos en coordenadas locales estado NO deformado
    // tiene origen en el nodo i NO deformado 
    float xil = 0.0;
    float yil = 0.0;
    float zil = 0.0;
    float xjl = xj - xi;
    float yjl = yj - yi ;
    float zjl = zj - zi;
    float xkl = xk - xi;
    float ykl = yk - yi;
    float zkl = zk - zi;
    
    // Transformar cada coordenada local del elemento sin deformar a la base i,j,k 
	float V1[3] = {xil, yil, zil};
	float pi[3];
	MdotV(r,V1,pi[0],pi[1],pi[2]);

	float V2[3] = {xjl, yjl, zjl};
	float pj[3];
	MdotV(r,V2,pj[0],pj[1],pj[2]);

	float V3[3] = {xkl, ykl, zkl};
	float pk[3];
	MdotV(r,V3,pk[0],pk[1],pk[2]);

	// Vectores de desplazamientos en el plano
    float di[3];
	di[0] = Pi[0] - pi[0];
	di[1] = Pi[1] - pi[1];
	di[2] = Pi[2] - pi[2];

	float dj[3];
	dj[0] = Pj[0] - pj[0];
	dj[1] = Pj[1] - pj[1];
	dj[2] = Pj[2] - pj[2];
    
	float dk[3];
	dk[0] = Pk[0] - pk[0];
	dk[1] = Pk[1] - pk[1];
	dk[2] = Pk[2] - pk[2];
    
	// Llama rutina para calcular las fuerzas en el elemento deformado 
	float ref[3][3], def[3][3];
	ref[0][0] = pi[0];
	ref[0][1] = pi[1];
	ref[0][2] = pi[2];

	ref[1][0] = pj[0];
	ref[1][1] = pj[1];
	ref[1][2] = pj[2];

	ref[2][0] = pk[0];
	ref[2][1] = pk[1];
	ref[2][2] = pk[2];

	def[0][0] = Pi[0];
	def[0][1] = Pi[1];
	def[0][2] = Pi[2];

	def[1][0] = Pj[0];
	def[1][1] = Pj[1];
	def[1][2] = Pj[2];

	def[2][0] = Pk[0];
	def[2][1] = Pk[1];
	def[2][2] = Pk[2];
	fuerzas(ref, def, nfuerzas);

	// Transformar las fuerzas calculadas a la base original 
	float tR[3][3];
	for(int i = 0; i<3; i++ )
		for(int j = 0; j<3; j++ )
		{
			tR[i][j] = R[j][i];
		}
	float nfn1[3], fn1[3] = {nfuerzas[0][0], nfuerzas[0][1], nfuerzas[0][2]};
	float nfn2[3], fn2[3] = {nfuerzas[1][0], nfuerzas[1][1], nfuerzas[1][2]};
	float nfn3[3], fn3[3] = {nfuerzas[2][0], nfuerzas[2][1], nfuerzas[2][2]};
	MdotV(tR, fn1, nfn1[0], nfn1[1], nfn1[2]);
	MdotV(tR, fn2, nfn2[0], nfn2[1], nfn2[2]);
	MdotV(tR, fn3, nfn3[0], nfn3[1], nfn3[2]);

	nfuerzas[0][0] = nfn1[0];
	nfuerzas[0][1] = nfn1[1];
	nfuerzas[0][2] = nfn1[2];

	nfuerzas[1][0] = nfn2[0];
	nfuerzas[1][1] = nfn2[1];
	nfuerzas[1][2] = nfn2[2];

	nfuerzas[2][0] = nfn3[0];
	nfuerzas[2][1] = nfn3[1];
	nfuerzas[2][2] = nfn3[2];
}
