#include <stdio.h>
#include <stdlib.h>

// Li(x) da formula de lagrange
double base_lagrange (int i, double X, double *x, int n) {
    double Li = 1.0;
    for (int j = 0; j < n; j++) {
        if (j != i) {
            Li *= (X - x[j]) / (x[i] - x[j]);
        }
    }
    return Li;
}

double lagrange(double *x, double *y, int n, double X) {
    
    double P = 0.0;
    for (int i = 0; i < n; i++) {
        P += y[i] * base_lagrange(i, X, x, n);
    }
    return P;
}

int main(void) {

    int n;
    printf("Digite o numero de pontos: ");
    scanf("%d", &n);

    if (n < 2) {
        printf("necessário pelo menos 2 pontos.\n");
        return 1;
    }

    // alocacao dos vetores
    double *x = malloc(n * sizeof(double));
    double *y = malloc(n * sizeof(double));

    if (x == NULL || y == NULL) {
        printf("Erro de alocacao.\n");
        return 1;
    }

    // leitura dos pontos
    printf("Digite os valores de x:\n");
    for (int i = 0; i < n; i++) {
        scanf("%lf", &x[i]);
    }

    printf("Digite os valores de y:\n");
    for (int i = 0; i < n; i++) {
        scanf("%lf", &y[i]);
    }

    int grau;
    printf("Digite o grau desejado do polinomio interpolador: ");
    scanf("%d", &grau);

    if (grau < 1 || grau >= n) {
        printf("grau invalido, ele deve ser menor que o numero de pontos");
        free(x);
        free(y);
        return 1;
    }

    int m = grau + 1; // numero de pontos a serem usados para o polinomio

    // ponto onde queremos interpolar
    double X;
    printf("Digite o ponto X para interpolar: ");
    scanf("%lf", &X);

    // calcula P(X)
    double PX = lagrange(x, y, m, X);

    printf("O valor do polinomio interpolado em X = %.6f eh %.6f\n", X, PX);

    // libera memoria
    free(x);
    free(y);

    return 0;
}