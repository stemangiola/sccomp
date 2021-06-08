#Intercept model
library(dplyr)
library(tibble)
data_generation_intercept_only = function(N, Tm, M, p_n) {
  # Sample N proportions for categories
  mu = runif(N, min = 0, max = 10)
  mu = mu / sum(mu)
  
  # Generate proportions
  I = matrix(c(rep(1, N)), nrow = N)
  beta = (p_n * (I - mu))
  alpha = (p_n * mu)
  proportion_begin = 0
  for (i in 1:length(I)) {
    proportion_i = rbeta(M, alpha[i], beta[i])
    proportion_begin = cbind(proportion_begin, proportion_i)
  }
  proportion = as.matrix(proportion_begin)[,1:N + 1]

  # Generate beta binomial cell count
  Cmn_md = matrix(c(rep(rep(0, N), M)), nrow = M)
  for (i in 1:M) {
    for (j in 1:N) {
      Cmn_md[i,j] =rbinom(1, Tm, proportion[i,j])
    }
  }
  count = as.data.frame(Cmn_md)
  colname=character(N)
  for (j in 1:N){
    col_id=j
    colname[j]=paste("categpry",col_id)
  }
  colnames(count)=colname
  rowname=character(M)
  for (i in 1:M){
    row_id=i
    rowname[i]=paste("subject",row_id)
  }
  rownames(count)=rowname
  return(count)
}
#example
data_generation_intercept_only(N = 10,
               Tm = 300,
               M = 5,
               p_n = exp(4))


# Covariates model


logsumexp <- function (x) {
  y = max(x)
  y + log(sum(exp(x - y)))
}

softmax <- function (x) {
  exp(x - logsumexp(x))
}
#input coefficients must be matrix
#input design matrix must be matrix that can left-multiple coefficient
data_generation = function(X, Tm, coefficients) {
  # Sample N proportions for categories
  N=ncol(coefficients)
  M=nrow(X)
  expected_values=X%*%coefficients
  mu_raw=(0)
  for (i in 1:M){
      unit=softmax(expected_values[i,1:N])
      mu_raw=rbind(mu_raw,unit)
    }
  mu=mu_raw[1:M+1,]
  # Generate precision, normally sampled from intercept
  intercept=coefficients[1,]
  precision=matrix(rep(NA,M*N),nrow = M)
  for (i in 1:N){
    precision[1:M,i]=exp(rnorm(M,mean=2.685 + intercept[i] * -0.744, sd=0.1))
  }
  # Generate proportions step 1, beta sampling 
  proportion = matrix(rep(NA,M*N),nrow=M, ncol=N)
  for (i in 1:M) {
    for (j in 1:N){
    proportion[i,j] = rbeta(1,mu[i,j]*precision[i,j],(1-mu[i,j])*precision[i,j])}}
  # Generate proportions step 2, make it to unit length 
  for (i in 1:M){ 
    proportion[i,]=proportion[i,]/rowSums(proportion)[i]}
  # Generate counts, binomial sampling 
  Cmn_md = matrix(c(rep(rep(NA, N), M)), nrow = M)
  for (i in 1:M) {
    for (j in 1:N){
    Cmn_md[i,j] =rbinom(1, Tm, proportion[i,j])}
  }
  count = as.data.frame(Cmn_md)
  # Name the data frame
  colname=character(N)
  for (j in 1:N){
    col_id=j
    colname[j]=paste("categpry",col_id)
  }
  colnames(count)=colname
  rowname=character(M)
  for (i in 1:M){
    row_id=i
    rowname[i]=paste("subject",row_id)
  }
  rownames(count)=rowname
  return(count)
}



# Calibration

# The same total counts per subject (M size vector, where M is the number of subjects)
bias=matrix(sample(-20:20,30),nrow=30)
Tm=matrix(rep(1000,30),nrow = 30)+bias

# Number of subjects 30, number of categories 20
# Setup coefficient to have same intercept (for simplicity), and zero slope
beta_0=matrix(rep(runif(1,min =  -3.193120, max=2.030292),20),ncol=20)

for (k in 1:100){
  # Create an empty data frame 
  df = data.frame(matrix(,30,20))
    # Name the data frame
  colname=character(20)
  for (j in 1:20){
    col_id=j
    colname[j]=paste("categpry",col_id)
  }
  colnames(df)=colname
  rowname=character(30)
  for (i in 1:30){
    row_id=i
    rowname[i]=paste("subject",row_id)
  }
  rownames(df)=rowname
  
  # Design matrix would have an intercept column and a factor of interest between -1 and 1
  X1=matrix(rep(1,30),nrow = 30)
  X2=matrix(runif(30,min = -1, max=1),nrow=30,ncol=1)
  # Design matrices are in size of 30*2
  X=cbind(X1,X2)
  
  # create the second row of coefficient
  beta_1=matrix(rep(runif(1,min =  -3.193120, max=2.030292),20),ncol=20)
  coefficients=rbind(beta_0,beta_1)
  
  df=data_generation(X,
                  Tm ,
                  coefficients)
  
  # label the output in each iteration
  assign(paste0("DF", k),df)
  }
