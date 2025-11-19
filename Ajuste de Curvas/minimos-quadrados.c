#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void gauss(int n, double **A, double *b, double *c) {
    for (int k = 0; k < n - 1; k++) {
        if (fabs(A[k][k]) < 1e-12) {
            printf("Erro: pivo muito pequeno");
            exit(1);
        }
        for (int i = k + 1; i < n; i++) {
            // calculo do fator
            double fator = A[i][k] / A[k][k];
            for (int j = k; j < n; j++) {
                A[i][j] = A[i][j] - (fator * A[k][j]); 
            }
            b[i] = b[i] - (fator * b[k]);
        }
    }

    // retrossubstituicao de baixo para cima
    for (int i = n - 1; i >= 0; i--) {
        double soma = b[i];
        for (int j = i + 1; j < n; j++) {
            soma = soma - (A[i][j] * c[j]);
        }
        c[i] = soma / A[i][i];
    }
}

// entradas: matriz A, vetor de valores xi, numero de pontos, grau do polinomio  
void construir_matriz_A(double **A, double *x, int n, int m) {
    // dimensao da matriz = grau do polinomio + 1
    int M = m + 1;

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            // cada elemento da matriz A é um somatorio 
            double soma = 0.0;
            for (int k = 0; k < n; k++) {
                // A[i][j] = sum x^(2m - i - j)
                soma += pow(x[k], (M - 1 - i) + (M - 1 - j));
            }
            A[i][j] = soma;
        }
    }
}

// entradas: vetor b, veotor de valores xi, vetor de valores f(xi), numero de pontos, grau do polinomio
void construir_vetor_b(double *b, double *x, double *y, int n, int m) {
    //dimensao do vetor = grau do polinomio + 1
    int M = m + 1;

    for (int i = 0; i < M; i++) {
        // cada elemento do vetor b é definido por uma soma de f(xi) * x^n
        double soma = 0.0;
        for (int k = 0; k < n; k++) {
            soma += y[k] * pow(x[k], M - i - 1);
        }
        b[i] = soma;
    }

}

// entradas: vetor de coeficientes, vetor de valores xi, vetor de valores f(xi), numero de pontos, grau do polinomio
void ajuste_polinomial(double *c, double *x, double *y, int n, int m) {
    int M = m + 1;

    // alocacao de memoria
    double **A = malloc(M * sizeof(double *));
    for (int i = 0; i < M; i++) {
        A[i] = calloc(M, sizeof(double));
    }

    double *b = calloc(M, sizeof(double));

    // construcao da matriz e dos vetores
    construir_matriz_A(A, x, n, m);
    construir_vetor_b(b, x, y, n, m);

    printf("\nMatriz A (%dx%d):\n", M, M);
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++)
            printf("%12.6lf ", A[i][j]);
        printf("\n");
    }

    printf("\nVetor b (%d):\n", M);
    for (int i = 0; i < M; i++)
        printf("%12.6lf\n", b[i]);
    printf("\n");

    // resolve sistema por eliminacao de gauss
    gauss(M, A, b, c);

    // libera memoria
    for (int i = 0; i < M; i++) {
        free(A[i]);
    }
    free(A);
    free(b);
}

int main() {
    int n, m;
    printf("Digite o numero de pontos: ");
    scanf("%d", &n);

    double *x = malloc(n * sizeof(double));
    double *y = malloc(n * sizeof(double));

    printf("Digite os valores de x:\n");
    for (int i = 0; i < n; i++)
        scanf("%lf", &x[i]);

    printf("Digite os valores de y:\n");
    for (int i = 0; i < n; i++)
        scanf("%lf", &y[i]);

    printf("Digite o grau do polinomio: ");
    scanf("%d", &m);

    int M = m + 1;
    double *c = calloc(M, sizeof(double));

    // Executa ajuste
    ajuste_polinomial(c, x, y, n, m);

    // Resultado
    printf("\nCoeficientes do polinomio:\n");
    for (int i = 0; i < M; i++)
        printf("c[%d] = %lf\n", i, c[i]);

    free(x);
    free(y);
    free(c);

    return 0;
}