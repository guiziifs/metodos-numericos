 #include <stdio.h>
 #include <stdlib.h>
 #include <math.h>

 double f(double x) {
    return pow(x, 2) * cos(x);
 }

 // entradas: funcao, intervalo a, intervalo b, numero de intervalos (par)
 double regra_simpson(double (*f)(double), double a, double b, int n) {

    // comprimento do intervalo
    double h = (b - a) / n;

    double soma = f(a) + f(b);

    for (int i = 1; i < n; i++) {
        double xi = a + (h*i);
        // se i e impar
        if(i % 2 == 1) {
            soma += 4 * f(xi);
        } else {
            // se for par
            soma += 2 * f(xi);
        }
    }

    return (h / 3.0) * soma;
 }

 int main(void) {

    double a, b;
    int n;

    printf("Digite o limite inferior a: ");
    scanf("%lf", &a);

    printf("Digite o limite superior b: ");
    scanf("%lf", &b);

    do {
        printf("Digite o numero de subintervalos (par): ");
        scanf("%d", &n);

    } while (n % 2 == 1);

    double I = regra_simpson(f, a, b, n);

    printf("\nAproximacao da integral = %.6lf\n", I);

    return 0;
}