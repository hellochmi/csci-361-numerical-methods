#include <iostream>
#include <cmath>
using namespace std;

double I(double x){
	return (1/pow(1-(pow(x,-0.2324)),1/0.2324));
}

double I_prime(double x){
	return -pow(I(x),1.2324)/pow(x,1.2324);
}

double arc_length(double x){
	x = I_prime(x);
	return .5*(asinh(x) + x*sqrt(1+pow(x,2)));
}

double sum(double a, double b){
	return arc_length(b) - arc_length(a);
}

double test(double x){
	return sum(I(x), x);
}

double trapezoid(double a, double b, int n, double (*f)(double x)){
	double h = (b-a)/n;
	double sum = 0.5*(f(a) + f(b));
	for (int i = 1; i < n; ++i) {
		double x_i = a + i*h;
		sum += f(x_i);
		}
	return sum*h;
}

int main(){
	double alpha = 0.2324;
	double beta = 0.4314;

	double target_length = 1 + beta;

	double test_iterate = 19.7830;
	double test_length = test(test_iterate);
	int count = 0;

	while(test_length < target_length){
		test_iterate += .0001;
		test_length = test(test_iterate);
		count += 1;
	}

	cout << "target_length: " << target_length << endl;
	cout << "test_length: " << test_length << endl;
	cout << "a: " << I(test_iterate) << endl;
	cout << "b: " << test_iterate << endl;
	cout << "iterations: " << count << endl;

	cout << "area of curve: " << trapezoid(16.4470, 23.8804, 5000, &I) << endl;

	cout << "arc length proof: " << sum(16.447,23.8804) << endl;

	// find root such that (a,b) and (b,a) and arclength = 
	

	return 0;
}