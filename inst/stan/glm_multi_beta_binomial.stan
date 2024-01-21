functions{
  vector Q_sum_to_zero_QR(int N) {
    vector [2*N] Q_r;

    for(i in 1:N) {
      Q_r[i] = -sqrt((N-i)/(N-i+1.0));
      Q_r[i+N] = inv_sqrt((N-i) * (N-i+1));
    }
    return Q_r;
  }

  row_vector sum_to_zero_QR(row_vector x_raw, vector Q_r) {
    int N = num_elements(x_raw) + 1;
    row_vector [N] x;
    real x_aux = 0;

    for(i in 1:N-1){
      x[i] = x_aux + x_raw[i] * Q_r[i];
      x_aux = x_aux + x_raw[i] * Q_r[i+N];
    }
    x[N] = x_aux;
    return x;
  }

  vector sum_to_zero_QR_vector(vector x_raw, vector Q_r) {
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

  array[] int rep_each(array[] int x, int K) {
    int N = size(x);
    array[N * K] int y;
    int pos = 1;
    for (n in 1:N) {
      for (k in 1:K) {
        y[pos] = x[n];
        pos += 1;
      }
    }
    return y;
  }

  row_vector average_by_col(matrix X) {
    int rows_X = rows(X);
    int cols_X = cols(X);
    row_vector[cols_X] means;
    
    for (j in 1:cols_X) {
      means[j] = mean(X[, j]);
    }
    
    return means;
  }

  real abundance_variability_regression(row_vector variability, row_vector abundance, array[] real prec_coeff, real prec_sd, int bimodal_mean_variability_association, real mix_p){

    real lp = 0;
    // If mean-variability association is bimodal such as for single-cell RNA use mixed model
    if(bimodal_mean_variability_association == 1){
      for(m in 1:cols(variability))
      lp += log_mix(mix_p,
      normal_lpdf(variability[m] | abundance[m] * prec_coeff[2] + prec_coeff[1], prec_sd ),
      normal_lpdf(variability[m] | abundance[m] * prec_coeff[2] + 1, prec_sd)  // -0.73074903 is what we observe in single-cell dataset Therefore it is safe to fix it for this mixture model as it just want to capture few possible outlier in the association
      );

      // If no bimodal
    } else {
      lp =  normal_lpdf(variability | abundance * prec_coeff[2] + prec_coeff[1], prec_sd );
    }

    return(lp);
  }

  real partial_sum_lpmf(int[] slice_y,
                        int start,
                        int end,
                        int[] exposure_array,
                        vector mu_array,
                        vector precision_array
                        ) {

return beta_binomial_lupmf(
    slice_y |
    exposure_array[start:end],
    (mu_array[start:end] .* precision_array[start:end]),
    (1.0 - mu_array[start:end]) .* precision_array[start:end]
  ) ;

}
}
data{
  int<lower=1> N;
  int<lower=1> M;
  int<lower=1> C;
  int<lower=1> A; // How many column in variability design\
  int<lower=1> A_intercept_columns; // How many intercept column in varibility design
  int<lower=1> B_intercept_columns; // How many intercept column in varibility design
  int<lower=1> Ar; // Rows of unique variability design
  array[N] int exposure;
  array[N,M] int y;
  matrix[N, C] X;
  matrix[Ar, A] XA; // The unique variability design
  matrix[N, A] Xa; // The variability design

  // Truncation
  int is_truncated;
  array[N,M] int truncation_up;
  array[N,M] int truncation_down;
  int<lower=1, upper=N*M> TNS; // truncation_not_size
  array[TNS] int<lower=1, upper=N*M> truncation_not_idx;
  int<lower=0, upper=1> is_vb;

  // Prior info
  array[2] real prior_prec_intercept;
  array[2] real prior_prec_slope;
  array[2] real prior_prec_sd;
  array[2] real prior_mean_intercept;
  array[2] real prior_mean_coefficients;

  // Exclude priors for testing purposes
  int<lower=0, upper=1> exclude_priors;
  int<lower=0, upper=1> bimodal_mean_variability_association;
  int<lower=0, upper=1> use_data;

  // Does the design icludes intercept
  int <lower=0, upper=1> intercept_in_design;

  // Random intercept
  int N_random_intercepts;
  int N_minus_sum;
  array[N_random_intercepts, 2] int paring_cov_random_intercept;
  int N_grouping;
  matrix[N, N_grouping] X_random_intercept;
  array[N_grouping, 2] int idx_group_random_intercepts;

  // LOO
  int<lower=0, upper=1> enable_loo;
}
transformed data{
  int<lower=1> grainsize = 1;    // Parallel chain
  vector[2*M] Q_r = Q_sum_to_zero_QR(M);
  real x_raw_sigma = inv_sqrt(1 - inv(M));
  matrix[N, C] Q_ast;
  matrix[C, C] R_ast;
  matrix[C, C] R_ast_inverse;
  array[N*M] int y_array;
  array[N*M] int truncation_down_array;
  array[N*M] int exposure_array;
  // EXCEPTION MADE FOR WINDOWS GENERATE QUANTITIES IF RANDOM EFFECT DO NOT EXIST
  int N_grouping_WINDOWS_BUG_FIX = max(N_grouping, 1);
  // thin and scale the QR decomposition
  Q_ast = qr_thin_Q(X) * sqrt(N - 1);
  R_ast_inverse = inverse(qr_thin_R(X) / sqrt(N - 1));
  // If I get crazy diagonal matrix omit it
  if(N_random_intercepts>0) { 
    if(max(R_ast_inverse)>1000 )
      print("sccomp says: The QR deconposition resulted in extreme values, probably for the correlation structure of your design matrix. Omitting QR decomposition.");
    Q_ast = X;
    R_ast_inverse = diag_matrix(rep_vector(1.0, C));
  }
  // Data vectorised
  y_array =  to_array_1d(y);
  truncation_down_array = to_array_1d(truncation_down);
  exposure_array = rep_each(exposure, M);
}
parameters{
  matrix[C, M-1] beta_raw_raw; // matrix with C rows and number of cells (-1) columns
  matrix[A, M] alpha; // Variability
  // To exclude
  array[2] real prec_coeff;
  real<lower=0> prec_sd;
  real<lower=0, upper=1> mix_p;
  // Random intercept // matrix with N_groupings rows and number of cells (-1) columns
  matrix[N_random_intercepts * (N_random_intercepts>0), M-1] random_intercept_raw;
  // sd of random intercept
  array[N_random_intercepts>0] real random_intercept_sigma_mu;
  array[N_random_intercepts>0] real random_intercept_sigma_sigma;
  row_vector[(M-1) * (N_random_intercepts>0)] random_intercept_sigma_raw;
  // If I have just one group
  array[N_random_intercepts>0] real zero_random_intercept;
}
transformed parameters{

  // Initialisation
  matrix[C,M] beta_raw;
  matrix[M, N] precision = (Xa * alpha)';
  matrix[C,M] beta;

  // Random effects
  matrix[N_minus_sum, M-1] random_intercept_minus_sum;
  row_vector[M-1] random_intercept_sigma;
  matrix[N_grouping, M-1] beta_random_intercept_raw;

  // locations distribution
  matrix[M, N] mu;

  // vectorisation
  vector[N*M] mu_array;
  vector[N*M] precision_array;

  for(c in 1:C)	beta_raw[c,] =  sum_to_zero_QR(beta_raw_raw[c,], Q_r);
  beta = R_ast_inverse * beta_raw; // coefficients on x

  // Calculate locations distribution
  mu = (Q_ast * beta_raw)';

  // random intercept
  if(N_random_intercepts>0 ){
    random_intercept_sigma = random_intercept_sigma_mu[1] + random_intercept_sigma_sigma[1] * random_intercept_sigma_raw;
    
    // Building the - sum, Loop across covariates
    for(a in 1:N_minus_sum){
      // Reset sum to zero
      row_vector[M-1] temp_random_intercept = rep_row_vector(0, M-1);
      // Loop across random intercept - 1
      for(n in 1:N_random_intercepts){
        if(paring_cov_random_intercept[n,1] == a)
        temp_random_intercept += random_intercept_raw[n];
      }
      // The sum to zero for each covariate
      random_intercept_minus_sum[a] = temp_random_intercept * -1;
    }
    // Build the beta_random_intercept_raw
    for(n in 1:N_grouping){
      
      // If primary parameter
      if(idx_group_random_intercepts[n,2]>0)
        beta_random_intercept_raw[idx_group_random_intercepts[n, 1]] =  random_intercept_raw[idx_group_random_intercepts[n, 2]]   .* exp(random_intercept_sigma / 3.0);
      
      // If sum to zero parameter
      else if(idx_group_random_intercepts[n,2]<0)
        beta_random_intercept_raw[idx_group_random_intercepts[n, 1]] = random_intercept_minus_sum[-idx_group_random_intercepts[n, 2]] .* exp(random_intercept_sigma / 3.0);
      
      
      // If a covariate has only one group
      else
        beta_random_intercept_raw[idx_group_random_intercepts[n, 1]] = rep_row_vector(zero_random_intercept[N_random_intercepts>0] * exp(random_intercept_sigma_mu[1] / 3.0), M-1) ;
    }

    // Update with summing mu_random_intercept
    mu = mu + append_row((X_random_intercept * beta_random_intercept_raw)', rep_row_vector(0, N));
  }

  // Calculate proportions
  for(n in 1:N)  mu[,n] = softmax(mu[,n]);

  // Convert the matrix m to a column vector in column-major order.
  mu_array = to_vector(mu);
  precision_array = to_vector(exp(precision));

}
model{

  // Fit main distribution
  if(use_data == 1){
    // target += beta_binomial_lpmf(
    //   y_array[truncation_not_idx] |
    //   exposure_array[truncation_not_idx],
    //   (mu_array[truncation_not_idx] .* precision_array[truncation_not_idx]),
    //   ((1.0 - mu_array[truncation_not_idx]) .* precision_array[truncation_not_idx])
    //   ) ;
      
    target +=  reduce_sum(
      partial_sum_lupmf,
      y_array[truncation_not_idx],
      grainsize,
      exposure_array[truncation_not_idx],
      mu_array[truncation_not_idx],
      precision_array[truncation_not_idx]
    );
  }

  // Priors
  if(exclude_priors == 0){
    
    // If interceopt in design or I have complex variability design
    // This would include the models 
    // composition ~ 1 + ...; composition ~ 0 + ...; 
    // variability ~ 1
    if(A == 1){
      target += abundance_variability_regression(
        alpha[1],
        beta[1], // average_by_col(beta[1:B_intercept_columns,]),
        prec_coeff,
        prec_sd,
        bimodal_mean_variability_association,
        mix_p
        );
    }
    else {
      // Loop across the intercept columns in case of a intercept-less design (covariate are intercepts)
      for(a in 1:A_intercept_columns)
        target += abundance_variability_regression(
          alpha[a],
          beta[a],
          prec_coeff,
          prec_sd,
          bimodal_mean_variability_association,
          mix_p
          );
          
        // Variability effect if the formula is more complex
        if(A>A_intercept_columns) for(a in (A_intercept_columns+1):A) alpha[a] ~ normal(beta[a] * prec_coeff[2], 2 );
    }

  }
  
  // If I don't have priors for overdispersion
  else{
     // Priors variability
     if(intercept_in_design || A > 1){
       for(a in 1:A_intercept_columns) alpha[a]  ~ normal( prec_coeff[1], prec_sd );
        if(A>A_intercept_columns) for(a in (A_intercept_columns+1):A) to_vector(alpha[a]) ~ normal ( 0, 2 );
     }
     // if ~ 0 + covariuate
     else {
       alpha[1]  ~ normal( prec_coeff[1], prec_sd );
     }
  }

  // // Priors abundance
  for(c in 1:B_intercept_columns) beta_raw_raw[c] ~ normal ( prior_mean_intercept[1], prior_mean_intercept[2] );
  if(C>B_intercept_columns) for(c in (B_intercept_columns+1):C) to_vector(beta_raw_raw[c]) ~ normal ( prior_mean_coefficients[1], prior_mean_coefficients[2]);

  // Hyper priors
  mix_p ~ beta(1,5);
  prec_coeff[1] ~ normal(prior_prec_intercept[1], prior_prec_intercept[2]);
  prec_coeff[2] ~ normal(prior_prec_slope[1],prior_prec_slope[2]);
  prec_sd ~ gamma(prior_prec_sd[1],prior_prec_sd[2]);

  // Random intercept
  if(N_random_intercepts>0){
    for(m in 1:(M-1))   random_intercept_raw[,m] ~ std_normal();
    random_intercept_sigma_raw ~ std_normal();
    random_intercept_sigma_mu ~ std_normal();
    random_intercept_sigma_sigma ~ std_normal();
     // If I have just one group
  zero_random_intercept ~ std_normal();
  }
}
generated quantities {
  matrix[A, M] alpha_normalised = alpha;
  // Rondom effect
  matrix[N_grouping_WINDOWS_BUG_FIX, M] beta_random_intercept;

  // LOO
  vector[TNS] log_lik = rep_vector(0, TNS);

  // These instructions regress out the effect of mean proportion to the overdispersion
  // This adjustment provide A overdispersion value that can be tested for a hypotheses for example differences between two conditions
  if(intercept_in_design){
    if(A > 1) for(a in 2:A) alpha_normalised[a] = alpha[a] - (beta[a] * prec_coeff[2] );
  }
  else{
    for(a in 1:A) alpha_normalised[a] = alpha[a] - (beta[a] * prec_coeff[2] );
  }

  // EXCEPTION MADE FOR WINDOWS GENERATE QUANTITIES IF RANDOM EFFECT DO NOT EXIST
  if(N_grouping==0) beta_random_intercept[1] = rep_row_vector(0.0, M);

  // Rondom effect
  else{
     beta_random_intercept[,1:(M-1)] = beta_random_intercept_raw;
  for(n in 1:N_grouping) beta_random_intercept[n, M] = -sum(beta_random_intercept_raw[n,]);
  }

  // LOO
  if(enable_loo==1)
    for (n in 1:TNS) {
      log_lik[n] = beta_binomial_lpmf(
        y_array[truncation_not_idx[n]] |
        exposure_array[truncation_not_idx[n]],
        (mu_array[truncation_not_idx[n]] .* precision_array[truncation_not_idx[n]]),
        ((1.0 - mu_array[truncation_not_idx[n]]) .* precision_array[truncation_not_idx[n]])
        ) ;
    }


}
