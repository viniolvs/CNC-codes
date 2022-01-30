#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "AutoV.h"

double modulo (double n)
{
	if (n<0)
		return -1.0*n;
	else 
		return n;
}

void scanmatriz (Matriz *x, int l, int c)
{
	int i,j;
	x->lin=l;
	x->col=c;
	for (i=0;i<l;i++)
		for (j=0;j<c;j++)
			scanf ("%lf",&x->m[i][j]);
}

void printmatriz (Matriz x)
{
	int i,j;
	for (i=0;i<x.lin;i++)
	{
		printf ("|");
		for (j=0;j<x.col;j++)
			printf ("%8.4lf ",x.m[i][j]);
		printf ("|\n");
	}
	printf ("\n");
}

Matriz identidade(int nm)
{
	Matriz I;
	int i,j;
	I.lin = I.col = nm;
	for (i=0 ; i<nm ; i++)
	{
		for (j=0 ; j<nm ; j++)
		{
			if (i==j)
				I.m[i][j] = 1.0;
			else 
				I.m[i][j] = 0.0;
		}	
	}
	return I;
}

Matriz inversa(Matriz A)
{
	Matriz I = identidade(A.lin);
	double p, pivo,m;

	int i,j,k,t,y;
	
	for (j = 0; j < A.col; j++)
    {
		pivo = A.m[j][j];
		//permuta a primeira linha com a segunda caso o pivo seja 0
		if (pivo == 0 && j<A.col-1)
		{
			double aux1[A.col],aux2[A.col];
			for (t=0, y=j; t<A.col ; t++)
			{
				aux1[t] = A.m[y][t];
				A.m[y][t] = A.m[y+1][t];
				A.m[y+1][t] = aux1[t];
				
				aux2[t] = I.m[y][t];
				I.m[y][t] = I.m[y+1][t];
				I.m[y+1][t] = aux2[t];
			}
			pivo = A.m[j][j];
		}
		p = 1.0;
		for (k = j; k < A.col; k++)
		{
			A.m[j][k] = (A.m[j][k]) / (pivo);
			I.m[j][k] = (I.m[j][k]) / (pivo);
		}
		for (i = 0; i < A.lin; i++)
		{
			if (i != j)
			{
				m = A.m[i][j];
				for (k = 0; k < A.col; k++)
				{
					A.m[i][k] = (A.m[i][k]) - (m * A.m[j][k]);
					I.m[i][k] = (I.m[i][k]) - (m * I.m[j][k]);
				}

			}
		}
	}
	return I;
	
}

Matriz *multmat(Matriz a, Matriz b)
{
	Matriz *c;
	c = (Matriz*)malloc(sizeof(Matriz));
	
	int i,j,k;
	if (a.col!=b.lin)
		printf ("Erro:Dimensoes da matriz nao coincidem.Impossivel calcular.\n");
	else 
	{
		c->lin=a.lin;
		c->col=b.col;
		for (i=0;i<a.lin;i++)
			for (j=0;j<b.col;j++)
				for (k=0;k<a.col;k++)
					c->m[i][j]+=a.m[i][k]*b.m[k][j];
	}
	return c;
}

Matriz somamat(Matriz a,Matriz b)
{
	int i,j;
	Matriz c;
	
	if (a.lin!=b.lin || a.col!=b.col)
	{
		printf ("Erro\n");
		exit (1);
	}
	else
	{
		c.lin=a.lin;
		c.col=a.col;
		for (i=0;i<a.lin;i++)
			for (j=0;j<a.col;j++)
				c.m[i][j] = a.m[i][j] + b.m[i][j];
	}
	return c;
}

double produtointerno(Matriz a, Matriz b)
{
	double prod=0;
	int i,j;
	
	for (i=0 ; i<a.lin ; i++)
		for (j=0 ; j<a.col ; j++)
			prod +=  a.m[i][j] * b.m[i][j];
	
	return prod;
}	

Matriz produtoescalar(Matriz u, double n)
{
	int i,j;
	Matriz c;
	c.lin = u.lin;
	c.col = u.col;
	for (i=0 ; i<u.lin ; i++)
		for (j=0 ; j<u.col ; j++)
			c.m[i][j] = n * u.m[i][j];
	return c;
}

double normaMax(Matriz u)
{
	double max = u.m[0][0];
	double mod;
	int i,j;
	for (i=0 ; i<u.lin ; i++)
	{
		for (j=0 ; j<u.col ; j++)
		{
			mod = modulo(u.m[i][j]);
			if ( mod > max)
				max = mod;
		}
	}
	return max;
}

double encontraq(Matriz A)
{
	double q,sum;
	double R[A.lin][2]; // vetor de intervalos
	int i,j,k,l;
	
	//zera o vetor R
	for (i=0 ; i<A.lin ; i++)
		for (j=0 ; j<2 ; j++)
			R[i][j] = 0.0;
	
	//encontra todos os intervalos de autovalores da Matriz e armazena em R
	for (k=0 ; k<2; k++)
	{
		for (i=0 ; i<A.lin ; i++)
		{
			R[i][k] = A.m[i][i];
			for (j=i+1,l=0,sum=0.0 ; l<A.col-1 ; l++)
			{
				
				if(j >= A.col)
					j = 0;
				sum+=modulo(A.m[i][j]);
				j++; 
			}
				
			if (k==0)
				R[i][k] -= sum;
			else 
				R[i][k] += sum; 
		}
	}
	//encontra q encontrando o maior intervalo e calculando sua média
	double maior=0.0;
	int aux=0;
	for (i=0 ; i<A.lin ; i++)
	{
		if ((R[i][0] - R[i][1]) > maior)
		{
			maior = (R[i][0] - R[i][1]);
			aux = i;
		}
	}
	q = (R[aux][0] + R[aux][1]) / 2.0;
	
	for (i=0 ; i<A.lin ; i++)
		printf("R%d: [%.1lf,%.1lf]\n",i+1,R[i][0],R[i][1]);
	
	return q;
}
