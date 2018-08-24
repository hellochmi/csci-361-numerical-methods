
//newton raphson alg

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