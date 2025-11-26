#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// funcao f(x) para a integracao
double f(double x) {
    return pow(x, 2) * cos(x);
}

// entradas: funcao, intervalo a, intervalo b, numero de subintervalos
double trapezios_repetidos(double (*f)(double), double a, double b, int n) {

    // altura
    double h = (b - a) / n;

    // f(x0) + f(xn)
    double soma = f(a) + f(b); 

    // 2 * somatorio de f(xi)
    for (int i = 1; i < n; i++) {
        double xi = a + i * h;
        soma += 2.0 * f(xi);
    }

    return (h / 2.0) * soma;
}

int main(void) {

    double a, b;
    int n;

    printf("Digite o limite inferior a: ");
    scanf("%lf", &a);

    printf("Digite o limite superior b: ");
    scanf("%lf", &b);

    printf("Digite o numero de subintervalos n: ");
    scanf("%d", &n);

    if (n <= 0) {
        printf("Erro: n deve ser maior que zero.\n");
        return 1;
    }

    double resultado = trapezios_repetidos(f, a, b, n);

    printf("Aproximacao da integral: %.6lf\n", resultado);

    return 0;
}