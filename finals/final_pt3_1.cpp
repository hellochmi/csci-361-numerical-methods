//#define func_quadratic
//#define func_cubic
#include <iostream>
#include <cmath>
using namespace std;


#define func_cum_norm
#ifdef func_quadratic
double func(double x){
return x*x;
};
void func(double x, double &f, double &fprime){
f=x*x;
fprime=2.0*x;
};
#endif // func_quadratic
#ifdef func_cubic
double func(double x){
return x*x*x;
};
void func(double x, double &f, double &fprime){
f=x*x*x;
fprime=3.0*x*x;
};
#endif // func_cubic
#ifdef func_cum_norm
double cum_norm(double x)
{
const double root = sqrt(0.5);
return 0.5*(1.0 + erf(x*root));
}
double func(double x)
{
return cum_norm(x);
}
void func(double x, double &f, double &fprime)
{
const double pi = 4.0*atan2(1.0,1.0);
f = cum_norm(x);
fprime = exp(-0.5*x*x)/sqrt(2.0*pi);
}
#endif // func_cum_norm

/**

Newton Raphson Root Finding Algorithm

@author Michelle Chung

*/

int root_bisection(double target, double tol_f, double tol_x, int max_iter,
	double x_low, double x_high, double & x, int & num_iter){
	x = 0.0;
	num_iter = 0;
	double y_low = func(x_low);
	double diff_y_low = y_low - target;
	if (abs(diff_y_low) <= tol_f){
		x = x_low;
		return 0;
	};
	double y_high = func(x_high);
	double diff_y_high = y_high - target;
	if (abs(diff_y_high) <= tol_f){
		x = x_high;
		return 0;
	};
	if (diff_y_low * diff_y_high > 0.0){
		x = 0;
		return 1;
	}
	for (num_iter = 1; num_iter < max_iter; ++num_iter) {
		x = (x_low + x_high)/2.0;
		double y = func(x);
		double diff_y = y - target;
		if (abs(diff_y) <= tol_f){
			return 0;
		};
		if(diff_y * diff_y_low > 0.0){
			x_low = x;
		}
		else {
			x_high = x;
		};
		if (abs(x_high - x_low) <= tol_x){
			return 0;
		};
	};
	x = 0;
	num_iter = max_iter;
	return 1;
};


int root_NR(double target, double tol_f, double tol_x, int max_iter, double x0,
double & x, int & num_iter){
	const double tol_fprime = 1.0e-12;
	double f = 0;
	double fprime = 0;
	x = x0;
	for (num_iter = 1; num_iter < max_iter; ++num_iter) {
		func(x,f,fprime);
		double diff_f = f - target;
		if (abs(diff_f) <= tol_f){
			return 0;
		}
		if (abs(fprime) <= tol_fprime){
			x = 0;
			return 1;
		}
		double delta_x = diff_f/fprime;
		if (abs(delta_x) <= tol_x){
			return 0;
		}
		x -= delta_x;
	};
	x = 0;
	num_iter = max_iter;
	return 1;
};


int main(){
	double target = func(2.0);
	double x = 0;
	double x0 = 0;
	int num_iter = 0;

	root_bisection(func(2.0), 1.0e-6, 1.0e-6, 100, 0.0, 5.0, x, num_iter);
	cout<<x<<endl;
	cout<<num_iter<<endl;
	cout<<"STOP"<<endl;

	//x = 0;
	//num_iter = 0;

	root_NR(func(2.0), 1.0e-6, 1.0e-6, 100, x0, x, num_iter);
	cout<<x<<endl;
	cout<<num_iter<<endl;

	int ans3 = cum_norm(2);
	cout<<ans3<<endl;

	return 0;
};









