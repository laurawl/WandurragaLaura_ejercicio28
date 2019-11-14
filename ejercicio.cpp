#include<iostream>
#include <vector>
#include <math.h> 
#include<cmath>
#include <fstream>
#include <vector>
//Ideas de: https://github.com/Diegofisie/Metodos_tarea4/blob/master/ODE/ODE.cpp y https://josemontanac.github.io/Laboratorio-Metodos-Computacionales/2/Ordinary%20differential%20equations%20(ODE).slides.html#/ que son del laboratorio y la magistral
using namespace std;
float PI = 3.14159265359;
float g = 9.8;
float c = 0.7;  
float m  = 0.2; 
float ht = 0.001; 
float t = 3.0;
int N = (int) t/ht;
vector<float> fun[2];
vector<float> posx_0;
vector<float> vel_0;

void inc(vector<float> v0,vector<float> v1){
	float mag = sqrt(v1[0]*v1[0]+v1[1]*v1[1]);
	vector<float> mov;
	float x = -(c/m)*(mag)*v1[0];
	float y = -g-(c/m)*(mag)*v1[1];
	mov.push_back(x);
	mov.push_back(y);
	fun[0]=v1;
	fun[1]=mov;
}

void  movimiento(string file_exp){
	// asignando condiciones inicciales para la velocidad
	vel_0.push_back(10);
	vel_0.push_back(10);

	
	// Se define un vector para asignar un valor inicial del movimiento 
	vector<float> M[N][2];
	M[0][0]=posx_0;
	M[0][1]=vel_0;
    
    float x_total = 0;

    ofstream file;
    file.open(file_exp.c_str());
   
/*Ahora se debe implementar el metodo de Runge Kutta para resulver la eqn*/
	for (int i = 1; i < N; ++i){
		vector<float> k1[2];
        vector<float> k2[2];
        vector<float> k3[2];
        vector<float> k4[2];
		vector<float> valorM0_2;
		vector<float> valorM1_2;
		vector<float> valorM0_3;
		vector<float> valorM1_3;
        vector<float> valorM0_4;
		vector<float> valorM1_4;
        vector<float> pendiente_0;
		vector<float> pendiente_1;
		vector<float> M0;
		vector<float> M1;
       
        //Valor de K1 para la solucion 
		inc(M[i-1][0],M[i-1][1]);
		k1[0] = fun[0];
		k1[1] = fun[1];
		// Calculo de K2
		float valorM0_2_x = (M[i-1][0])[0]+(ht/2.0);
		float valorM0_2_y = (M[i-1][0])[1]+(ht/2.0);
		float valorM1_2_x = (M[i-1][1])[0]+((k1[1])[0])*(ht/2.0);
		float valorM1_2_y = (M[i-1][1])[1]+((k1[1])[1])*(ht/2.0);
		valorM0_2.push_back(valorM0_2_x);
		valorM0_2.push_back(valorM0_2_y);
		valorM1_2.push_back(valorM1_2_x);
		valorM1_2.push_back(valorM1_2_y);

		inc(valorM0_2,valorM1_2);
		k2[0] = fun[0];
		k2[1] = fun[1];
		//Se calcula el valor de k3
		
		float valorM0_3_x = (M[i-1][0])[0]+(ht/2.0);
		float valorM0_3_y = (M[i-1][0])[1]+(ht/2.0);
		float valorM1_3_x = (M[i-1][1])[0]+((k2[1])[0])*(ht/2.0);
		float valorM1_3_y = (M[i-1][1])[1]+((k2[1])[1])*(ht/2.0);
		valorM0_3.push_back(valorM0_3_x);
		valorM0_3.push_back(valorM0_3_x);
		valorM1_3.push_back(valorM1_3_x);
		valorM1_3.push_back(valorM1_3_y);

		inc(valorM0_2,valorM1_3);
		k3[0] = fun[0];
		k3[1] = fun[1];
		
		float valorM0_4_x = (M[i-1][0])[0]+(ht);
		float valorM0_4_y = (M[i-1][0])[1]+(ht);
		float valorM1_4_x = (M[i-1][1])[0]+((k3[1])[0])*(ht);
		float valorM1_4_y = (M[i-1][1])[1]+((k3[1])[1])*(ht);
		valorM0_4.push_back(valorM0_4_x);
		valorM0_4.push_back(valorM0_4_y);
		valorM1_4.push_back(valorM1_4_x);
		valorM1_4.push_back(valorM1_4_y);

		inc(valorM0_4,valorM1_4);
		k4[0] = fun[0];
		k4[1] = fun[1];
		float k1x = ((k1[0])[0]+2.0*(k2[0])[0]+2.0*(k3[0])[0]+(k4[0])[0])/6.0;
		float k1y = ((k1[0])[1]+2.0*(k2[0])[1]+2.0*(k3[0])[1]+(k4[0])[1])/6.0;
		float k2x = ((k1[1])[0]+2.0*(k2[1])[0]+2.0*(k3[1])[0]+(k4[1])[0])/6.0;
		float k2y = ((k1[1])[1]+2.0*(k2[1])[1]+2.0*(k3[1])[1]+(k4[1])[1])/6.0;

	 	pendiente_0.push_back(k1x);
	 	pendiente_0.push_back(k1y);
	 	pendiente_1.push_back(k2x);
	 	pendiente_1.push_back(k2y);
	 	M0.push_back((M[i-1][0])[0]+ht*pendiente_0[0]);
	 	M0.push_back((M[i-1][0])[1]+ht*pendiente_0[1]);
	 	M1.push_back((M[i-1][1])[0]+ht* pendiente_1[0]);
	 	M1.push_back((M[i-1][1])[1]+ht* pendiente_1[1]);
	 	M[i][0]=M0;
		M[i][1]=M1;
		file<<(M[i-1][0])[0]<<" "<<(M[i-1][0])[1]<<"\n";
    }
        file.close();
}
int main() {	
	posx_0.push_back(0.0);
	posx_0.push_back(0.0);
	movimiento("datos.dat");

}