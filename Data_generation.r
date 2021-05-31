library(dplyr)
library(tibble)
logsumexp <- function (x) {
  y = max(x)
  y + log(sum(exp(x - y)))
}

softmax <- function (x) {
  exp(x - logsumexp(x))
}
#input coefficients must be matrix
#input design matrix must be matrix that can left-multiple coefficient
data_generation = function(X, Tm, coefficients, p_n) {
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
  # Generate proportions step 1 
  proportion = matrix(rep(NA,M*N),nrow=M, ncol=N)
  Cmn_md = matrix(c(rep(rep(0, N), M)), nrow = M)
  for (i in 1:M) {
    for (j in 1:N){
  # Generate proportions step 2
    proportion[i,j] = rbeta(1,mu[i,j]*p_n,(1-mu[i,j])*p_n)
  # Generate counts
    Cmn_md[i,j] =rbinom(1, Tm, proportion[i,j])}
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
# Example 1
trial=matrix(runif(40,min = 0,max=10),nrow=5,ncol=2)
data_generation(X=trial,
               Tm = 50,
               coefficients=matrix(rep(0.3,20),nrow=2,ncol=10),
               p_n = exp(3))
# Example 2
X= matrix(c(1, 1, 1, 1, 1, 0, 1, 0, 1, 0), ncol = 2)
coefficients =matrix(c(1, 2, 3, 4, 5, 1, 2 ,3, 4, 5, 1, 1, 1, -1, -1, -1, 1, 1, 2, -2), nrow=2, byrow = T)
data_generation(X,
                Tm = 200,
                coefficients,
                p_n = exp(3))