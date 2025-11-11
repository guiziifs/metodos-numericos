#include <stdio.h>
#include <math.h>

double func(double x) {
    return x*x*x - x - 2;
}

double secantes(double (*f)(double), double x0, double x1, double epslon, int *iteracao, int max_iteracoes) {
    double x2, f0, f1;

    while (*iteracao < max_iteracoes) {
        f0 = f(x0);
        f1 = f(x1);

        if (fabs(f1 - f0) < 1e-12) {
            printf("Divisão por zero\n");
            return NAN;
        }

        x2 = x1 - f1 * (x1 - x0) / (f1 - f0);

        (*iteracao)++;

        // parada
        if (fabs(x2 - x1) < epslon || fabs(f(x2)) < epslon) {
            return x2;
        }

        x0 = x1;
        x1 = x2;
    }

    printf("Numero max de iteracoes\n");
    return x1;
}

int main() {
    double x0, x1, epslon;
    int iteracoes = 0;

    printf("Digite o primeiro chute inicial (x0): ");
    scanf("%lf", &x0);

    printf("Digite o segundo chute inicial (x1): ");
    scanf("%lf", &x1);

    printf("Digite o criterio de parada (epslon): ");
    scanf("%lf", &epslon);

    double raiz = secantes(func, x0, x1, epslon, &iteracoes, 100);

    if (!isnan(raiz)) {
        printf("\nRaiz aproximada: %.4f\n", raiz);
        printf("Numero de iteraoes: %d\n", iteracoes);
    }

    return 0;
}