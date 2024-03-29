functions{

  real dirichlet_multinomial2_lpmf(array[] int y, vector alpha) {
    	 real alpha_plus = sum(alpha);

    return lgamma(alpha_plus) + sum(lgamma(alpha + to_vector(y)))
                - lgamma(alpha_plus+sum(y)) - sum(lgamma(alpha));
  }

matrix vector_array_to_matrix(array[] vector x) {
		matrix[size(x), rows(x[1])] y;
		for (m in 1:size(x))
		  y[m] = x[m]';
		return y;
}

vector Q_sum_to_zero_QR(int N) {
    vector [2*N] Q_r;

    for(i in 1:N) {
      Q_r[i] = -sqrt((N-i)/(N-i+1.0));
      Q_r[i+N] = inv_sqrt((N-i) * (N-i+1));
    }
    return Q_r;
  }

  vector sum_to_zero_QR(vector x_raw, vector Q_r) {
    int N = num_elements(x_raw) + 1;
    vector [N] x;
    real x_aux = 0;

    for(i in 1:N-1){
      x[i] = x_aux + x_raw[i] * Q_r[i];
      x_aux = x_aux + x_raw[i] * Q_r[i+N];
    }
    x[N] = x_aux;
    return x;
  }

  real beta_regression_lpdf(array[] vector p, matrix X, matrix beta, vector phi){

		real lp = 0;
    matrix[num_elements(p[1]), num_elements(p[,1])] mu = (X * beta)';
		for(i in 1:cols(mu)) mu[,i] = softmax(mu[,i]);
		for(i in 1:cols(mu)) {

     	lp += beta_lpdf(p[i] | (mu[,i] .* phi), ((1.0 - mu[,i]) .* phi) );

		}
		return (lp);
	}

}
data{
	int N;
	int M;
	array[N] simplex[M] y;
	matrix[N, 2] X;

}
transformed data{
	int C = 2;
	vector[2*M] Q_r = Q_sum_to_zero_QR(M);
  real x_raw_sigma = inv_sqrt(1 - inv(M));
}
parameters{
	array[C] vector[M-1] beta_raw;
	vector[M] precision;

	// To exclude
  array[2] real prec_coeff;
  real<lower=0> prec_sd;
}
transformed parameters{
		array[C] vector[M] beta;
	  for(c in 1:C)	beta[c] =  sum_to_zero_QR(beta_raw[c], Q_r);


}
model{


	 y ~ beta_regression( X, vector_array_to_matrix(beta), exp(precision) );

	 precision ~ normal( beta[1] * prec_coeff[2] + prec_coeff[1], prec_sd);
   prec_sd ~ normal(0,2);

	 for(i in 1:C) beta_raw[i] ~ normal(0, x_raw_sigma * 5);


}
