#include <stdio.h>
#include <math.h>

void Taylor(int a, int b, double q, double h, long n, double epsilon, double an, double sum, double i, double* xvalues, double* yvalues);
double InterpolateLagrangePolynomial (double x, double* x_values, double* y_values, int size);
void TaylorOtrisovka (int a, int b, double q, double h, long n, double epsilon, double an, double sum, double i, double* xvalues, double* yvalues) ;
void Otrisovka (double l, int mashtab);

int main() {
    double i, q, an, sum;
    long n;
    int k=10;
    double epsilon = pow(10, -4);
    double a=0, b=2;
    double h = (double)(b - a) / k;
    double yvalues[6]={0.000000, 0.428435, 0.742168, 0.910282, 0.976399, 0.995371};
    double xvalues[6]={0.0, 0.4, 0.8, 1.2, 1.6, 2.0};
    printf("Task1.\n");
    Taylor(a, b, q, h, n, epsilon, an, sum, i, xvalues, yvalues);
    printf("\n\n");
    printf("Task2.\n");
    TaylorOtrisovka(a, b, q, h, n, epsilon, an, sum, i, xvalues, yvalues);
    return 0;
}
    

void Taylor (int a, int b, double q, double h, long n, double epsilon, double an, double sum, double i, double* xvalues, double* yvalues) {
    for (i = a; i <= b; i += h) {
        sum = 0, n = 0, an = i;
        while (fabs(an) > epsilon) {
            sum += an;
            q = -(i * i * (2 * n + 1)) / ((n + 1) * (2 * n + 3));
            an *= q;
            n++;
        }
        sum *= 2 / sqrt(M_PI);
        float l=(float) InterpolateLagrangePolynomial(i, xvalues, yvalues ,6);
        printf("x = %.1lf, y = %lf, Ln(x) = %lf, err = %lf\n", i, sum, l, fabs(l-sum));
    }  
    return;
}

double InterpolateLagrangePolynomial (double x, double* x_values, double* y_values, int size)
{
	double lagrange_pol = 0;
	double basics_pol;

	for (int i = 0; i < size; i++)
	{
		basics_pol = 1;
		for (int j = 0; j < size; j++)
		{
			if (j == i) continue;
			basics_pol *= (x - x_values[j])/(x_values[i] - x_values[j]);		
		}
		lagrange_pol += basics_pol*y_values[i];
	}
	return lagrange_pol;
}

void TaylorOtrisovka (int a, int b, double q, double h, long n, double epsilon, double an, double sum, double i, double* xvalues, double* yvalues) {
    printf("Lagrange graph.\n");
    for (i = a; i <= b; i += h) {
        sum = 0, n = 0, an = i;
        while (fabs(an) > epsilon) {
            sum += an;
            q = -(i * i * (2 * n + 1)) / ((n + 1) * (2 * n + 3));
            an *= q;
            n++;
        }
        sum *= 2 / sqrt(M_PI);
        float l=(float) InterpolateLagrangePolynomial(i, xvalues, yvalues ,6);
        double l1=l;
        double l2=fabs(l-sum);
        Otrisovka(l1, 50);
    }
    printf("\n");
    printf("Error plot.\n");
    for (i = a; i <= b; i += h) {
        sum = 0, n = 0, an = i;
        while (fabs(an) > epsilon) {
            sum += an;
            q = -(i * i * (2 * n + 1)) / ((n + 1) * (2 * n + 3));
            an *= q;
            n++;
        }
        sum *= 2 / sqrt(M_PI);
        float l=(float) InterpolateLagrangePolynomial(i, xvalues, yvalues ,6);
        double l1=l;
        double l2=fabs(l-sum);
        Otrisovka(l2, 60000);
    }
}

void Otrisovka (double l, int mashtab) {
    int pas;
        pas=(int) (mashtab*l);
        while (pas>0) {
          printf("*");
          pas=pas-1;
        }
        printf("*\n");  
}