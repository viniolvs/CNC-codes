#define size 5

typedef struct
{
	double m[size][size];
	int lin,col;
}Matriz;

double modulo (double n);

void scanmatriz (Matriz *x, int l, int c);

void printmatriz (Matriz x);

Matriz identidade(int nm);

Matriz inversa(Matriz A);

Matriz *multmat(Matriz a, Matriz b);

Matriz somamat(Matriz a,Matriz b);

double produtointerno(Matriz a, Matriz b);

Matriz produtoescalar(Matriz u, double n);

double normaMax(Matriz u);

double encontraq(Matriz A);
