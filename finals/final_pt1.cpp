#include <iostream>
#include <cmath>
using namespace std;

// performs integration at the correct intervals for e1, e2
double trapezoid(double a, double b, int n, double (*e_1)(double x), double (*e_2)(double x)){
	double h = (b-a)/n;
	double sum = 0.5*(e_2(a) + e_1(b));
	for (int i = 1; i < n; ++i) {
		double x_i = a + i*h;
		if(x_i >= -0.5686 && x_i <= 1.0722)
		sum += e_2(x_i);
		if(x_i > 1.0722 && x_i <= 1.2324)
		sum += e_1(x_i);
		}
	return sum*h*2;
}

// general trapezoid integration function
double t(double a, double b, int n, double (*f)(double x)){
	double h = (b-a)/n;
	double sum = 0.5*(f(a) + f(b));
	for (int i = 1; i < n; ++i) {
		double x_i = a + i*h;
		sum += f(x_i);
		}
	return sum*h;
}

double e_1(double x){
	return sqrt(.5*(1-(pow(x-.2324,2))));
}

double e_2(double x){
	return sqrt(.25*(1-(pow(x-.4314,2))));
}

double closed(double x){
	return (asin(x) + x*sqrt(1-pow(x,2)));	
}
double integrate(double a, double b, double frac, double x_i){
	return frac*(closed(b-x_i)-closed(a-x_i));
}

int main(){
	cout << "A: " << trapezoid(-0.5686, 1.2324, 5000, &e_1, &e_2)<<endl;;
	cout << "e1 " << t(-.7676,1.2324,500,&e_1)*2<<endl;
	cout << "e2 " << t(.4314-1,.4314+1,500,&e_2)*2<<endl;

	double a = integrate(-.5686, 1.0722, .25, .4314);
	double b = integrate(1.0722, 1.2324, (sqrt(.5)/2), .2324);
	double c = 2*(a + b);
	cout << "a " << a << endl;
	cout << "b " << b << endl;
	cout << "c " << c << endl;

	return 0;
}