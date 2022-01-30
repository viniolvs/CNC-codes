#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "AutoV.h"

int main()
{
	Matriz A,B,C,D;
	Matriz *u,v,w[2]; // vetores
	int k,n=0,passos=6;
	int i,j,l,z;
	double lambda[2]; lambda[0]=0; lambda[1]=0;
	double lb;
	double erro=100.0, e, q;
	double x;
	
	z=3;
	e=0.001;
	
	scanmatriz(&A,z,z);
	scanmatriz(&v,z,1);
	
	double mu[A.lin];
	
	encontraq(A);
	q=-0.5;
	printmatriz(A);
	
	C = produtoescalar(identidade(A.lin),-1*q);
	printmatriz (C);
	
	D = somamat(A,C);
	printmatriz(D);
	
	//B  = inversa(D);
	scanmatriz(&B,A.lin,A.col);
	printmatriz(B);
	

	w[0] = v; k = 0; l=0;
	while (erro > e )//&& k < passos)
	{
		u = multmat(B,w[0]);
		if (k>=2)
		{
			lambda[n+1] = produtointerno(*u,v)/produtointerno(v,v);
				mu[l] = lambda[n+1];
				l++;
		}
		if (k>=3)
			erro = modulo((modulo(lambda[n]) - modulo(lambda[n+1])) / modulo(lambda[n+1]));
		x = 1.0/normaMax(*u);
		w[1] = produtoescalar(*u,x);
		w[0] = w[1];
		lambda[n] = lambda[n+1];

		k++;
	}
	
	for (i=0;i<A.lin;i++)
	{
		lb = (1.0/mu[i]) + q;
		printf("MU%d = %lf\n",i,mu[i]);
		printf("Autovalor %d = %lf\n",i,lb);
	}
	printf("Autovetor w:\n");
	printmatriz(w[0]);
	printf("Autovetor u:\n");
	printmatriz(*u);
	
	return 0;
}
