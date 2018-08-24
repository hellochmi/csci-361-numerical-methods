#include <iostream>
#include <cmath>
using namespace std;

double trapezoid(double a, double b, int n, double (*f)(double x)){
	// perform validation checks n > 0, etc.
	double h = (b-a)/n;
	double sum = 0.5*(f(a) + f(b));
	for (int i = 1; i < n; ++i) {
		double x_i = a + i*h;
		sum += f(x_i);
		}
	return sum*h*M_PI;
}

double simpson(double a, double b, int n, double (*f)(double x)){
	// perform validation checks n > 0, etc.
	int i;
	double h = (b-a)/n;
	double sum1 = f(a) + f(b);
	double sum2 = 0.0;
	for (i = 1; i < n; i += 2) {
		double x_i = a +i*h;
		sum2 += f(x_i);
	}
	double sum3 = 0.0;
	for (i = 2; i < n; i += 2) {
		double x_i = a +i*h;
		sum3 += f(x_i);
	}
	double result = (sum1 +4.0*sum2 +2.0*sum3)*h/3.0;
	return result*M_PI;
}

double midpoint(double a, double b, int n, double (*f)(double x)){
	double h = (b-a)/n;
	double sum = 0.0;
	for (int i = 0; i < n; ++i) {
		double x_i = a + (i+0.5)*h;
		sum += f(x_i);
		}
	return sum*h*M_PI;
}

double triangle(double x){
	if (x <= 0 && x <= 0.5){
		return 4 * x;
	}
	if (x > 0.5 && x <= 1){
		return 4 * (1 - x);
	}
	return 1;
}

double romberg(double j, double k){
	// perform validation checks, jmax > 0 etc
	int jmax = 256; // input parameter
	double E[jmax]; // array of length at least jmax
	double **Rom; // memory must be allocated
	// extended trapezoid rule is called to compute E[j] for j=0, ..., jmax-1
	// set Rom[0][0]
	double *R0 = Rom[0];
	R0[0] = E[0];
	for (int j = 1; j < jmax; ++j) {
		double *tmp1 = Rom[j];
		double *tmp2 = Rom[j-1];
		tmp1[0] = E[j];
		for (int k = 1; k <= j; ++k) {
			double pow4k = (1 << (2*k)); // pow4k = 4^k
			tmp1[k] = (pow4k*tmp1[k-1] -tmp2[k-1]) /(pow4k-1.0);
		}
	}
}

double I(double x){
	return (1 + x*x)/(sqrt(1 - 0.5*(x*x)));
}

double I_2(double x){
	double a = 0.9;
	return (1 + x*x)/(sqrt(1 - ((a*a)*(x*x))));
}

double bessel(double x){
	return cos(200*sin(x) - x);
}

int main(){
	int n = 256;

	cout << "midpoint: " << midpoint(0.0, M_PI, n, &bessel)<<endl;;
	cout << "trapezoid: " << trapezoid(0.0, M_PI, n, &bessel)<<endl;;
	cout << "simpson: " << simpson(0.0, M_PI, n, &bessel)<<endl;
	cout << "romberg: " << romberg(1, 0)<<endl;

	//cout << simpson(0.0, 1.0, n, &I_2) << endl;
	return 0;
}