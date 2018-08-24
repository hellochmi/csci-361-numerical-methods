//bisection algorithm

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