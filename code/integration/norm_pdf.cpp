#include <iostream>
using namespace std;

double cum_norm(double x)
{
	const double root = sqrt(0.5);
	return 0.5*(1.0 + erf(x*root));
}

int main(){
	return 0;
}