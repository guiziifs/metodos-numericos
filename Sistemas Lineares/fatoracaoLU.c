#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define EPS 1e-12

void verifica_pivo_valido (double **A, int n) {
    for (int i = 0; i < n; i++) {
        if (fabs(A[i][i]) < EPS) {
            printf("Erro: pivo A[%d][%d] = %.6f e muito pequeno (possivel divisao por 0)\n", i, i, A[i][i]);
            exit(EXIT_FAILURE);
        }
    }
}

// fatora A em matriz Upper (triangular superior) e Lower (contem os multiplicadores da eliminacao de Gauss)
void fatoracaoLU(double **A, double **L, double **U, int n) {
    // Inicializa L como identidade e copia A para U
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            L[i][j] = (i == j) ? 1.0 : 0.0;
            U[i][j] = A[i][j];
        }
    }

    // Verifica pivôs iniciais (conforme pedido)
    verifica_pivo_valido(U, n);

    // Eliminação de Gauss com pivotamento parcial
    for (int k = 0; k < n - 1; k++) {
        // Escolhe o maior pivô na coluna
        int pivo = k;
        double max = fabs(U[k][k]);
        for (int i = k + 1; i < n; i++) {
            if (fabs(U[i][k]) > max) {
                max = fabs(U[i][k]);
                pivo = i;
            }
        }

        if (max < EPS) {
            printf("Erro: pivô nulo detectado na coluna %d (sistema singular ou mal condicionado).\n", k);
            exit(EXIT_FAILURE);
        }

        // Troca de linhas se necessário
        if (pivo != k) {
            double *temp = U[k];
            U[k] = U[pivo];
            U[pivo] = temp;

            // troca as partes já preenchidas de L (colunas 0..k-1)
            for (int j = 0; j < k; j++) {
                double tmp = L[k][j];
                L[k][j] = L[pivo][j];
                L[pivo][j] = tmp;
            }
        }

        // Eliminação
        for (int i = k + 1; i < n; i++) {
            double m = U[i][k] / U[k][k];
            L[i][k] = m;

            for (int j = k; j < n; j++) {
                U[i][j] -= m * U[k][j];
            }

            // força zero explícito na posição eliminada (evita resíduos numéricos)
            U[i][k] = 0.0;
        }
    }
}

// resolve a substituicao Ly = b
void substituicaoLy(double **L, double *b, double *y, int n) {
    for (int i = 0; i < n; i++) {
        double soma = 0.0;
        for (int j = 0; j < i; j++) {
            soma += L[i][j] * y[j];
        }
        y[i] = b[i] - soma;
    }
}

void substituicaoUx(double **U, double *x, double *y, int n) {
    for(int i = n - 1; i >= 0; i--) {
        double soma = 0.0;
        for (int j = i + 1; j < n; j++) {
            soma += U[i][j] * x[j];
        }

        if (fabs(U[i][i]) < EPS) {
            printf("Erro: pivô U[%d][%d] = %.6f é muito pequeno durante retro-substituição.\n", i, i, U[i][i]);
            exit(EXIT_FAILURE);
        }

        x[i] = (y[i] - soma) / U[i][i];
    }
}

int main() {
    int n;

    printf("Digite o numero de equacoes: ");
    if (scanf("%d", &n) != 1 || n <= 0) {
        printf("Entrada invalida.\n");
        return 1;
    }

    // ---------------
    // alocacao de mem
    // ---------------

    // vetores de ponteiros
    double **A = (double **)malloc(n * sizeof(double *));
    double **L = (double **)malloc(n * sizeof(double *));
    double **U = (double **)malloc(n * sizeof(double *));
    // aloca as linhas da matriz
    for (int i = 0; i < n; i++) {
        A[i] = (double *)malloc(n * sizeof(double));
        L[i] = (double *)malloc(n * sizeof(double)); 
        U[i] = (double *)malloc(n * sizeof(double));
    }

    // vetores
    double *b = (double *)malloc(n * sizeof(double));
    double *x = (double *)malloc(n * sizeof(double));
    double *y = (double *)malloc(n * sizeof(double));

    // ----------------
    // leituras
    // ----------------

    // matriz A
    printf("Digite os coeficientes da matriz A linha por linha: \n");
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            if (scanf("%lf", &A[i][j]) != 1) {
                printf("Erro na leitura de A[%d][%d].\n", i, j);
                return 1;
            }
        }
    }

    // vetor b
    printf("Digite os coeficientes do vetor b: \n");
    for(int i = 0; i < n; i++) {
        if (scanf("%lf", &b[i]) != 1) {
            printf("Erro na leitura de b[%d].\n", i);
            return 1;
        }
    }

    fatoracaoLU(A, L, U, n);

    substituicaoLy(L, b, y, n);
    substituicaoUx(U, x, y, n);

    printf("\nMatriz L (multiplicadores):\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) printf("%.4f ", L[i][j]);
        printf("\n");
    }

    printf("\nMatriz U (triangular superior):\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) printf("%.4f ", U[i][j]);
        printf("\n");
    }

    printf("\nVetor y (solucao de L*y = b):\n");
    for (int i = 0; i < n; i++) printf("y[%d] = %.4f\n", i, y[i]);

    printf("\nSolucao do sistema x (Ux = y):\n");
    for (int i = 0; i < n; i++) printf("x[%d] = %.4f\n", i, x[i]);

    // Libera memória
    for (int i = 0; i < n; i++) {
        free(A[i]); free(L[i]); free(U[i]);
    }
    free(A); free(L); free(U);
    free(b); free(y); free(x);

    return 0;
}
