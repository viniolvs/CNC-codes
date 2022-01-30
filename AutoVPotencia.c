#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "AutoV.h"

int main()
{
	Matriz A; 
	Matriz *u,v,w[2]; // vetores
	double lambda[2]; lambda[0]=0; lambda[1]=0;
	double erro=100.0, e;
	double x;
	int k=0,n=0,passos=10;
	
	e=0.001;
	
	scanmatriz(&A,3,3);
	scanmatriz(&v,3,1);
	w[0] = v;
	encontraq(A);
	while (erro > e && k < passos)
	{
		u = multmat(A,w[0]);
		if (k>=2)
			lambda[n+1] = produtointerno(*u,v)/produtointerno(v,v);
		if (k>=3)
			erro = modulo((modulo(lambda[n]) - modulo(lambda[n+1])) / modulo(lambda[n+1]));
		x = 1.0/normaMax(*u);
		w[1] = produtoescalar(*u,x);
		w[0] = w[1];
		lambda[n] = lambda[n+1];
		k++;
	}
	
	printf ("Matriz A:\n");
	printmatriz(A);
	printf("Vetor v:\n");
	printmatriz(v);
	printf("Autovalor = %lf\n",lambda[n]);
	printf("Autovetor w:\n");
	printmatriz(w[0]);
	printf("Autovetor u:\n");
	printmatriz(*u);
	
	
	return 0;
}
