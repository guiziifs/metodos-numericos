#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

// gera matriz de Vandermonde m×m usando apenas m pontos
void matrizVandermonde(double **A, double *x, int m) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            A[i][j] = pow(x[i], m - 1 - j);
        }
    }
}

// eliminação de gauss simples (sem pivotamento)
void gauss(double **A, double *b, double *x, int n) {

    for (int k = 0; k < n - 1; k++) {

        if (fabs(A[k][k]) < 1e-12) {
            printf("Erro: pivo zero (sem pivotamento).\n");
            return;
        }

        for (int i = k + 1; i < n; i++) {
            double fator = A[i][k] / A[k][k];

            for (int j = k; j < n; j++) {
                A[i][j] -= fator * A[k][j];
            }

            b[i] -= fator * b[k];
        }
    }

    // retro-substituição
    for (int i = n - 1; i >= 0; i--) {

        double soma = b[i];

        for (int j = i + 1; j < n; j++) {
            soma -= A[i][j] * x[j];
        }

        if (fabs(A[i][i]) < 1e-12) {
            printf("Erro: pivô zero na retro-substituição.\n");
            x[i] = 0.0;
        } else {
            x[i] = soma / A[i][i];
        }
    }
}

double avalia_polinomio(double *coef, int grau, double x0) {
    double resultado = 0.0;

    // for de 0 ate grau para pegar todos os coeficientes
    for (int i = 0; i <= grau; i++) {
        int exp = grau - i; //coef em ordem decrescente
        resultado += coef[i] * pow(x0, exp);
    }

    return resultado;
}


int main(void) {

    int n;
    printf("Digite o numero de pontos disponiveis: ");
    scanf("%d", &n);

    if (n < 2) {
        printf("E preciso pelo menos 2 pontos.\n");
        return 0;
    }

    // aloca vetores x e y
    double *x = malloc(n * sizeof(double));
    double *y = malloc(n * sizeof(double));

    printf("Digite os valores de x:\n");
    for (int i = 0; i < n; i++) {
        scanf("%lf", &x[i]);
    }

    printf("Digite os valores de f(x):\n");
    for (int i = 0; i < n; i++) {
        scanf("%lf", &y[i]);
    }

    // usuário escolhe grau
    int grau;
    printf("Digite o grau do polinomio interpolador (max = %d): ", n - 1);
    scanf("%d", &grau);

    if (grau < 1 || grau >= n) {
        printf("Grau invalido. Deve estar entre 1 e %d.\n", n - 1);
        return 0;
    }

    // número de equações = grau + 1
    int m = grau + 1;

    // aloca matriz Vandermonde m×m
    double **A = malloc(m * sizeof(double *));
    for (int i = 0; i < m; i++)
        A[i] = malloc(m * sizeof(double));

    // aloca vetor solução de tamanho m
    double *coef = calloc(m, sizeof(double));

    // vetor b usa apenas os primeiros m valores de y
    double *b = malloc(m * sizeof(double));
    for (int i = 0; i < m; i++)
        b[i] = y[i];

    // constrói Vandermonde 
    matrizVandermonde(A, x, m);

    // resolve sistema
    gauss(A, b, coef, m);

    printf("\nCoeficientes do polinomio (grau %d):\n", grau);
    for (int i = 0; i < m; i++) {
        printf("a%d = %lf\n", i, coef[i]);
    }

    // avaliar o polinomio em um ponto xp
    double xp;
    printf("\nDigite um ponto xp para avaliar o polinomio: ");
    scanf("%lf", &xp);

    double yp = avalia_polinomio(coef, grau, xp);
    printf("P(%lf) = %lf\n", xp, yp);

    // libera memoria
    for (int i = 0; i < m; i++) {
        free(A[i]);
    }

    free(A);
    free(x);
    free(y);
    free(b);
    free(coef);

    return 0;
}
