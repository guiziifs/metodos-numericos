#include <stdio.h>
#include <math.h>
#include <stdbool.h>

// aleterar a funcao correspondente
double func(double x) {
    return x*x - 2;
}

bool bolzano (double fa, double fb) {
    if (fa * fb >= 0) {
        return true; // nao satisfaz bolzano
    } else {
        return false; // satisfaz bolzano
    }
}

double falsa_posic(double (*f)(double), double a, double b, double epslon, int *iteracao) {
    // certifica que o teorema de bolzano se cumpre
    if (bolzano(f(a), f(b))) {
        printf("Intervalo não satisfaz teorema de bolzano");
        return NAN;
    }

    double fa = f(a);
    double fb = f(b);
    double xm, fxm;
    *iteracao = 0;

    do {
        xm = ((a*fb) - (b*fa)) / (fb - fa);
        fxm = f(xm);
        (*iteracao)++;

        if (fabs(fxm) < epslon) {
            return xm;
        }

        // atualiza intervalo
        if (fa * fxm < 0) {
            b = xm;
            fb = fxm;
        } else {
            a = xm;
            fa = fxm;
        }
    } while (fabs(b - a) > epslon);

    return xm;
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
    double raiz = falsa_posic(func, a, b, epslon, &i);

    if (!isnan(raiz)) {
        printf("Raiz aproximada: %.4f\n", raiz);
        printf("Numero de iteracoes: %d", i);
    }

    return 0;

}