#include <stdio.h>
#include <math.h>

double func (double x) {
    return x*x - 2;
}

double dfunc (double x) {
    return 2*x;
}

double newton_raphson (double (*f)(double), double (*df)(double), double x0, double epslon, int *iteracao) {
    double x1;
    *iteracao = 0;

    while (1) {
        double fx = f(x0);
        double dfx = df(x0);

        if (fabs(dfx) < 1e-12) {
            printf("Derivada muito prox de zero");
            return NAN;
        }

        x1 = x0 - fx / dfx;
        (*iteracao)++;

        if (fabs(x1 - x0) < epslon || fabs(fx) < epslon) {
            break;
        }
        x0 = x1;
    }
    return x1;
}

int main() {
    double x0, epsilon;

    printf("Digite o chute inicial (x0): ");
    scanf("%lf", &x0);
    printf("Digite o criterio de parada (epsilon): ");
    scanf("%lf", &epsilon);

    int iter = 0;
    double raiz = newton_raphson(func, dfunc, x0, epsilon, &iter);

    if (!isnan(raiz)) {
        printf("\nRaiz aproximada: %.4f\n", raiz);
        printf("Numero de iteracoes: %d\n", iter);
    }

    return 0;
}