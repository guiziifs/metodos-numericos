#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void gauss(double **A, double *b, double *x, int n) {
    // matriz aumentada [A | b]
    for (int k = 0; k < n - 1; k++) {

        int max = k;
 
        for(int i = k + 1; i < n; i++) {
            if (fabs(A[i][k]) > fabs(A[max][k])) {
                max = i;
            }
        }

        if (max != k) {
            double *temp = A[k];
            A[k] = A[max];
            A[max] = temp;

            double t = b[k];
            b[k] = b[max];
            b[max] = t;
        }

        // eliminacao
        for (int i = k + 1; i < n; i++) {
            if (fabs(A[k][k]) < 1e-12) {
                continue;
            }
            double fator = A[i][k] / A[k][k];
            for (int j = k; j < n; j++) {
                A[i][j] -= fator * A[k][j];
            }
            b[i] -= fator *b[k];
        }


    }

    // retrosubstituicao
    for (int i = n - 1; i >= 0; i--) {
        double soma = b[i];
        for (int j = i + 1; j < n; j++) {
            soma -= A[i][j] * x[j];
        }
        if (fabs(A[i][i]) < 1e-12) {
            x[i] = 0.0;
        } else {
            x[i] = soma / A[i][i];
        }
    }
}

int main() {
    int n;
    printf("Digite o número de equações: ");
    scanf("%d", &n);

    // alocacao da matriz A
    double **A = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) {
        A[i] = (double *)malloc(n * sizeof(double));
    }

    // alocacao dos vetores x e b
    double *x = (double *)malloc(n * sizeof(double));
    double *b = (double *)malloc(n * sizeof(double));

    // leitura da matriz
    printf("Digite os coeficientes da matriz A: \n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            scanf("%lf", &A[i][j]);
        }
    }

    // leitrua do vetor b
    printf("Digite os termos indep: \n");
    for(int i = 0; i < n; i++) {
        scanf("%lf", &b[i]);
    }

    gauss(A, b, x, n);

    printf("\nSolucao do sistema: \n");
    for(int i = 0; i < n; i++) {
        printf("x[%d] = %.4f\n", i, x[i]);
    }

    for (int i = 0; i < n; i++) {
        free(A[i]);
    }

    free(A);
    free(b);
    free(x);

    return 0;
}
