data 
{
  vector[99] log_y;
  vector[99] dll;
  int material[99];
  vector[99] t;
  real lambda_mu;
  real lambda_sigma;
  real b_mu;
  real b_sigma;
  real sigma_mu;
  real sigma_sigma;
}

parameters {vector[3] lambda;vector[3] b;vector[3] sigma;}

model 
{
  vector[99] log_y1; 
  for (i in 1:99)
  {
    log_y1[i]=b[material[i]]-lambda[material[i]]*t[i];
    if (log_y[i]<=dll[i]) {target+=normal_lcdf(dll[i]|log_y1[i],sigma[material[i]]);} 
    else 
    {
        target+=log_diff_exp(normal_lcdf(log_y[i]|log_y1[i],sigma[material[i]]),
                     normal_lcdf(log_y[i]-0.25|log_y1[i],sigma[material[i]]));
    }
  }
 
  b~normal(b_mu,b_sigma);
  lambda~normal(lambda_mu,lambda_sigma);
  sigma[1]~normal(sigma_mu,sigma_sigma) T[0,];
  sigma[2]~normal(sigma_mu,sigma_sigma) T[0,];
  sigma[3]~normal(sigma_mu,sigma_sigma) T[0,];
  
}

