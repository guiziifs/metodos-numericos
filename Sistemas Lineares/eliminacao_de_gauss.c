#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void gauss(double **A, double *b, double *x, int n) {
    // matriz aumentada [A | b]
    for (int k = 0; k < n - 1; k++) {

        // pivotamento parcial 
        int max = k;
        
        // procura maior magnitude na coluna k 
        for(int i = k + 1; i < n; i++) {
            if (fabs(A[i][k]) > fabs(A[max][k])) {
                max = i;
            }
        }

        //  se uma linha abaixo de max tem maior pivo, troca a linha k com a linha max
        // se max == k significa que a linha atual ja e o maior valor
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
            // calculo do fator
            double fator = A[i][k] / A[k][k];
            for (int j = k; j < n; j++) {
                A[i][j] -= fator * A[k][j];
            }
            b[i] -= fator *b[k];
        }


    }

    // retrosubstituicao
    // comeca de baixo para cima
    for (int i = n - 1; i >= 0; i--) {
        // inicializa a soma com o termo indep
        double soma = b[i];
        // subtraimos com os valores ja conhecidos
        for (int j = i + 1; j < n; j++) {
            soma -= A[i][j] * x[j];
        }
        // se o pivo for zero, seta o resultado como 0 
        if (fabs(A[i][i]) < 1e-12) {
            x[i] = 0.0;
        } else {
            // divide a soma pelo pivo
            x[i] = soma / A[i][i];
        }
    }
}

int main() {
    int n;
    printf("Digite o numero de equacoes: ");
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
