#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

double a_func(double x, double h){
	return (h*h*x)-(h*h)-(2*x*x); 
}

double b_func(double x, double h){
	return (x*x)-(h/2)*x; 
}

double c_func(double x, double h){
	return (x*x)+(h/2)*x; 
}

int tridiagonal(const int n,
				const std::vector<double> & a,
				const std::vector<double> & b,
				const std::vector<double> & c,
				const std::vector<double> & rhs,
				std::vector<double> & x)
{
	const double tol = 1.0e-14;
	x.clear();
	if ((n < 1) || (a.size() < n) || (b.size() < n) || (c.size() < n)
				|| (rhs.size() < n)) return 1; // fail
	double alpha[n]; // temporary storage

	x.resize(n, 0.0);

	// initial equation i = 0
	int i = 0;
	double gamma = a[i];
	if (abs(gamma) <= tol){
		cout << "FAILED AT 1" << endl;
		return 1; // fail
	}
	x[i] = rhs[i]/gamma;
	alpha[i] = c[i]/gamma;

	// forward pass: elimination
	for (i = 1; i < n-1; ++i) {
		gamma = a[i] - b[i]*alpha[i-1];
		if (abs(gamma) <= tol){
			cout << "FAILED AT 2" << endl;
		 	return 1; // fail
		 }
		x[i] = (rhs[i] - b[i]*x[i-1])/gamma;
		alpha[i] = c[i]/gamma;
	}
	// solve final equation i = n-1
	i = n-1;
	gamma = a[i] - b[i]*alpha[i-1];
	if (std::abs(gamma) <= tol){
		cout << "FAILED AT 3" << endl;
		return 1;
	} // fail
	x[i] = (rhs[i] - b[i]*x[i-1])/gamma;

	// backward substitution
	for (i = n-2; i >= 0; --i) {
		x[i] -= alpha[i]*x[i+1];
	}
	return 0;
}

int bisection(int a, int b, std::vector<double> & x){
	double tol = 0.0001;
	int midpoint = (b - a)/2;

	double root = x[midpoint];
	if(abs(root) <= tol) return midpoint;

	if (x[a]*x[b] > 0) return bisection(midpoint, b, x);
	else return bisection(a, midpoint, x);
}

int main(){

	int n = 19000;

	// declare vectors
	std::vector<double> x(n);   // interval vector
	std::vector<double> rhs(n,0); // rhs vector
	std::vector<double> a(n);   // a vector
	std::vector<double> b(n);   // b vector
	std::vector<double> c(n);   // c vector

	// populate vectors
	for (int i = 0; i < n; i++){
		x[i] = i;
	}

	for(int i = 0; i < n; i++){
		a[i] = a_func(x[i], 1);
		b[i] = b_func(x[i], 1);
		c[i] = c_func(x[i], 1);
	}

	// change endpoints
	rhs[0] = 0.2324;
	rhs[n] = 0.4314;
	a[0]   = 1.0;
	a[n]   = -1.0;
	b[0]   = 0.0;
	b[n]   = 0.0;
	c[0]   = 0.0;
	c[n]   = 1.0;

	tridiagonal(n, a, b, c, rhs, x);

	// root finding tridiagonal algorithm

	// tolerance
	double m = 10000.0;

	// delcare vectors
	std::vector<double> x_root(m);     // new interval vector
	std::vector<double> rhs_root(m,0); // rhs vector
	std::vector<double> a_root(m);     // a vector
	std::vector<double> b_root(m);     // b vector
	std::vector<double> c_root(m);     // c vector

	// populate vectors
	for(int i = 0; i <= m; i++){
		x_root[i] = 6+(.0001*i);
	}

	for(int i = 0; i < m; i++){
		a_root[i] = a_func(x_root[i], 1);
		b_root[i] = b_func(x_root[i], 1);
		c_root[i] = c_func(x_root[i], 1);
	}

	// update endpoints

	rhs_root[0] = -0.170576;
	rhs_root[m] = 0.0666176;
	a_root[0]   = 1.0;
	a_root[m]   = 1.0;
	b_root[0]   = 0.0;
	b_root[m]   = 0.0;
	c_root[0]   = 0.0;
	c_root[m]   = 0.0;

	tridiagonal(m, a_root, b_root, c_root, rhs_root, x_root);

	// call bisection, add interval to x_low find x-value
	cout << 6+(.0001*bisection(0, m, x_root)) << endl;

	return 0;
}