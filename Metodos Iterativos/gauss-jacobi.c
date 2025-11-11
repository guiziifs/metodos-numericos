#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define EPS 1e-12

// construcao da matriz C e do vetor g
void matC_vetg(double **A, double **C, double *b, double *g, int n) {
    for (int i = 0; i < n; i++) {
        if(fabs(A[i][i]) < EPS) {
            printf("Elemento diagonal 0 (divisao por zero)");
            exit(EXIT_FAILURE);
        }
        for (int j = 0; j < n; j++) {
            if (i == j) {
                C[i][j] = 0.0;
            } else {
                C[i][j] = - (A[i][j] / A[i][i]);
            }
        }
        g[i] = b[i] / A[i][i];
    }
}

double norma_linha(double **C, int n) {
    double max = 0.0;
    for (int i = 0; i < n; i++) {
        double soma = 0.0;
        for (int j = 0; j < n; j++) {
            soma += fabs(C[i][j]);
        }
        if (soma > max) {
            max = soma;
        }
    }
    return max;
}

double norma_coluna(double **C, int n) {
    double max = 0.0;
    for (int j = 0; j < n; j++) {
        double soma = 0.0;
        for (int i = 0; i < n; i++) {
            soma += fabs(C[i][j]);
        }
        if (soma > max) {
            max = soma;
        }
    }
    return max;
}

// verificacao da convergencia do metodo
int converge(double **C, int n) {
    double norm1 = norma_linha(C, n);
    double norm2 = norma_coluna(C, n);

    if (norm1 < 1.0 || norm2 < 1.0) {
        return 1;
    } else {
        return 0;
    }
}

// calcula o erro relativo entre iteracoes
double calc_erro(double *x, double *x_ant, int n) {
    double num = 0.0, den = 0.0;
    for (int i = 0; i < n; i++) {
        num += pow(x[i] - x_ant[i], 2);
        den += pow(x[i], 2);
    }
    return sqrt(num / den);
}

// entradas: A, b, x(n), precisao, n_max de iteracoes
void gauss_jacobi(double **A, double *b, double *x, double epslon, int iter_max, int n) {
    double **C = (double **)malloc(n * sizeof(double *));
    double *g = (double *)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        C[i] = (double *)malloc(n * sizeof(double));
    }

    // construir C e g
    matC_vetg(A, C, b, g, n);

    // verificar convergência
    if (!converge(C, n)) {
        printf("Encerrando execucao.\n");
        for (int i = 0; i < n; i++) free(C[i]);
        free(C);
        free(g);
        return;
    }

    // iteracao
    double *x_ant = (double *)malloc(n * sizeof(double));

    int k = 0;
    double erro = 1.0;
    while (erro > epslon && k < iter_max) {
        for (int i = 0; i < n; i++) {
            x_ant[i] = x[i];
        }

        for (int i = 0; i < n; i++) {
            double soma = 0.0;
            for (int j = 0; j < n; j++) {
                soma += C[i][j] * x_ant[j];
            }
            x[i] = soma + g[i];
        }

        erro = calc_erro(x, x_ant, n);
        k++;

    }

    printf("\nSolucao aproximada apos %d iteracoes:\n", k);
    for (int i = 0; i < n; i++)
        printf("x[%d] = %.4f\n", i, x[i]);

    // Libera memória
    for (int i = 0; i < n; i++) free(C[i]);
    free(C);
    free(g);
    free(x_ant);
}

int main() {
    int n, iter_max;
    double epslon;

    printf("Digite o numero de equacoes: ");
    scanf("%d", &n);

    // Alocação
    double **A = (double **)malloc(n * sizeof(double *));
    double *b = (double *)malloc(n * sizeof(double));
    double *x = (double *)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++)
        A[i] = (double *)malloc(n * sizeof(double));

    // Entrada da matriz A
    printf("Digite os coeficientes da matriz A linha por linha:\n");
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            scanf("%lf", &A[i][j]);

    // Entrada de b
    printf("Digite o vetor b:\n");
    for (int i = 0; i < n; i++)
        scanf("%lf", &b[i]);

    // Entrada do chute inicial
    printf("Digite o vetor inicial x0:\n");
    for (int i = 0; i < n; i++)
        scanf("%lf", &x[i]);

    printf("Digite a precisao desejada (epslon): ");
    scanf("%lf", &epslon);

    printf("Digite o numero maximo de iteracoes: ");
    scanf("%d", &iter_max);

    gauss_jacobi(A, b, x, epslon, iter_max, n);

    // Libera memória
    for (int i = 0; i < n; i++) free(A[i]);
    free(A);
    free(b);
    free(x);

    return 0;
}