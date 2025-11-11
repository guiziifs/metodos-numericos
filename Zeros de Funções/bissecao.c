#include <stdio.h>
#include <math.h>

// aleterar a funcao correspondente
double func(double x) {
    return x*x*x - x - 2;
}

double bissec(double (*f)(double), double a, double b, double epslon, int *iteracao) {
    // certifica que o teorema de bolzano se cumpre
    if (f(a) * f(b) >= 0 ) {
        printf("Intervalo não satisfaz teorema de bolzano");
        return NAN;
    }

    double ponto_medio;
    *iteracao = 0;
    while ((b - a) / 2.0 > epslon) {
        ponto_medio = (a + b) / 2.0;
        double fm = f(ponto_medio);
        
        (*iteracao)++;

        // valor absoluto de f(ponto_medio) < epslon
        if (fabs(fm) < epslon) {
            return ponto_medio;
        }

        if (f(a) * fm < 0) {
            b = ponto_medio;
        } else {
            a = ponto_medio;
        }
    }
    return (a + b) / 2.0;
}

int main() {
    double a, b, epslon;

    printf("Digite o inicio do intervalo: ");
    scanf("%lf", &a);
    printf("Digite o fim do intervalo: ");
    scanf("%lf", &b);
    printf("Digite o cirtério de parada: ");
    scanf("%lf", &epslon);

    int i = 0;
    double raiz = bissec(func, a, b, epslon, &i);

    if (!isnan(raiz)) {
        printf("Raiz aproximada: %.4f\n", raiz);
        printf("Numero de iteracoes: %d", i);
    }

    return 0;

}