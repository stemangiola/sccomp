data{
  int N; // Number of subjects
  int M; // Number of categories
  int C;
  int A;
  array[N] int exposure;
  matrix[N, C] X;
  matrix[A, A] XA;

  matrix[C,M] beta;

  real<lower=0> variability_multiplier;

}

parameters{

  array[2] real prec_coeff;
  real<lower=0> prec_sd;
}

generated quantities{
  array[N, M] int counts_uncorrected;
  matrix[N, M] counts;
  matrix[A,M] alpha;
  matrix[M,N] mu = (X * beta)';
  matrix[M,N] precision;

  matrix[A,M] beta_intercept_slope;
  // Vector of the generated exposure
  array[N] real generated_exposure;

  // matrix[A,M] alpha_intercept_slope;

  // All this because if A ==1 we have ocnversion problems
  // This works only with two discrete groups
  if(A == 1) beta_intercept_slope = to_matrix(beta[A,], A, M, 0);
  else beta_intercept_slope = (XA * beta[1:A,]);
  // if(A == 1)  alpha_intercept_slope = alpha;
  // else alpha_intercept_slope = (XA * alpha);

  // PRECISION REGRESSION
  for(a in 1:A) for(m in 1:M) alpha[a,m] = normal_rng( beta_intercept_slope[a,m] * prec_coeff[2] + prec_coeff[1], prec_sd);

  precision = (X[,1:A] * alpha)';

  // Precision adjustment
  precision = precision - log(variability_multiplier);

  for(i in 1:N) mu[,i] = softmax(mu[,i]);
  for(i in 1:cols(mu)) {
    counts_uncorrected[i,] = beta_binomial_rng(
      exposure[i],
      mu[,i] .* exp(precision[,i]),
      (1.0 - mu[,i]) .* exp(precision[,i])
      );
  }

  // Calculate the generated exposure
  for(n in 1:N) generated_exposure[n] = sum(counts_uncorrected[n]);
  for(n in 1:N) counts[n] = to_row_vector(counts_uncorrected[n]) / generated_exposure[n] * exposure[n];

}
