#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

double fatorial(double n)
{
	double fat;
	for (fat=1;n>1;n--)
		fat = fat * n;

	return fat;
}	
	

//Altera-se Fxy para cada situação
double MetEuler(double h, double Yx0, double parada)
{
	double Y[2],X,Fxy;
	int n=0;
	
	printf("Yn+1 = Yn + h * f(Xn,Yn)\n\n");
	
	Y[0] = Yx0;
	X = 0;
	while (X < parada)
	{
		printf ("Y(%.1lf) = %lf\n",X,Y[n]);
		
		Fxy = -0.05 * Y[n];	//f(Xn,Yn)
		
		Y[n+1] = Y[n] + h * Fxy; 			//Yn+1 = Yn + h * Fxy
		Y[n] = Y[n+1];
		X += h;
	}
	return Y[0]; 
}


//Altera-se K1 e K2 para cada situação
double MetHeun(double h, double Yx0, double parada)
{
	double Y[2],X,K1,K2;
	int n = 0;
	
	printf("Yn+1 = Yn + h/2 * (K1 + K2)\n\n");
	
	Y[0] = Yx0;
	X = 0;
	while (X < parada)
	{
		printf ("Y(%.1lf) = %lf\n",X,Y[n]);
		
		K1 = -0.02*Y[n];				//f(Xn,Yn)    
		K2 = -0.02 * (Y[n]+(h * K1));  	//f(Xn+h , h*f(Xn,Yn))
		
		printf ("K1 = %lf\nK2 = %lf\n\n",K1,K2);
		Y[n+1] = Y[n] + h/2 * (K1 + K2); 
		Y[n] = Y[n+1];
		X += h;
	}
	return Y[0];
}

double MetHeun3(double h, double Yx0)
{
	double Y[2],X,K1,K2;
	int n = 0,anos=0;
	
	printf("Yn+1 = Yn + h/2 * (K1 + K2)\n\n");
	
	Y[0] = Yx0;
	X = 0;
	while (Y[n] > 0.92)
	{
		printf ("Y(%.1lf) = %lf\n",X,Y[n]);
		
		K1 = (-1.244*pow(10.0,-4.0))*Y[n];				//f(Xn,Yn)    
		K2 = (-1.244*pow(10.0,-4.0)) * (h * K1);  		//f(Xn+h , h*f(Xn,Yn))
		
		//printf ("K1 = %lf\nK2 = %lf\n\n",K1,K2);
		Y[n+1] = Y[n] + h/2 * (K1 + K2); 
		Y[n] = Y[n+1];
		X += h;
		anos++;
	}
	printf("Anos = %d\n",anos);
	return Y[n];
}

//Altera-se Fxy
double MetTaylor (double h, double Yx0, double parada, int ordem)
{
	double Y[2], X, Fxy[ordem];
	double k, G;
	int n=0, j;
	
	Y[0] = Yx0;
	X = 0;	
	while (X < parada)
	{
		
		printf ("Y(%.1lf) = %lf\n",X,Y[n]);
		//Modificar aqui
		Fxy[0] = Y[n];  
		Fxy[1] = 1;  
		for (k=1.0,j=0,G=0.0; j<ordem ; k+=1.0, j++)
			G += (pow(h,k)/fatorial(k)) * Fxy[j];

		Y[n+1] = Y[n] + G; 
		Y[n] = Y[n+1];
		X += h; 
	}
	return Y[0];
}

//Utilizar para ordem 3 ou ordem 4
//Alterar as expressões no vetor K[ordem] 
double MetRungeKutta (double h, double Yx0, double parada, int ordem)
{
	double Y[2], X, K[ordem];
	int n=0;
	
	Y[0] = Yx0;
	X = 0;
	
	//Ordem 3
	if (ordem == 3)
	{
		while (X < parada)
		{
			printf ("Y(%.1lf) = %lf\n",X,Y[n]);
			
			K[0] = h * (X - Y[n] + 2) ; //h*f(Xn,Yn)
			K[1] = h * (X+h/2.0 - (Y[n]+K[0]/2.0) + 2.0) ;	  		//h*f(Xn + h/2 , Yn + K1/2)
			K[2] = h * (X+3.0/4.0*h - (Y[n]+3.0/4.0*K[1]) + 2.0 );  //h*f(Xn + 3/4*h , Yn + 3/4*K2)
			
			printf("K1 = %lf\nK2 = %lf\nK3 = %lf\n",K[0],K[1],K[2]);
			
			Y[n+1] = Y[n] + (2.0/9.0)*K[0] + (1.0/3.0)*K[1] + (4.0/9.0)*K[2];
			
			
			
			X+=h;
			Y[n] = Y[n+1];
		}
	}
	//ordem 4
	else if (ordem == 4)
	{
		while (X < parada)
		{
			printf ("Y(%.1lf) = %lf\n",X,Y[n]);
			
			K[0] = h * 1; 	  	  //h*f(Xn,Yn)
			K[1] = h * 1;	  	  //h*f(Xn + h/2 , Yn + K1/2)
			K[2] = h * 1;		  //h*f(Xn + h/2 , Yn + K2/2)
			K[3] = h * 1;		  //h*f(Xn + h/2 , Yn + K3)
			
			Y[n+1] = Y[n] + 1.0/6.0*(K[0] + 2.0*K[1] + 2.0*K[2] + K[3]);
			
			X+=h;
			Y[n] = Y[n+1];
		}
	}

	return Y[0];	
}

int main()
{
	double resposta, h, Yx0, parada;
	int anos;
	int ordem;
	
	h = 0.2;
	Yx0 = 1;
	parada = 0.4;
	ordem = 2;
	resposta = MetTaylor(h,Yx0,parada,ordem);
	printf("Y(%.1lf) = %.4lf\n",parada,resposta);
	return 0;
}
