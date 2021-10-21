//
  // This Stan program defines the Bayesian preference model, with a vector of 
// preferences 'y' modeled as Bernoulli distributed with parameter 'prob'.
//
  // Learn more about model development with Stan at:
  //
  //    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//
  
  
  
  // The input data  
data {
  int<lower=1> n_pairs; // number of pairs compared  
  int y[n_pairs]; // user preferences 
  int index_i[n_pairs]; // i index - 1st solution for comparison 
  int index_j[n_pairs]; // j index - 2nd solution for comparison 
  int<lower=2> n_points; // number of feasible points 
  real obj_1[n_points]; // normalized DO-level at Fortuna objective
  real obj_2[n_points]; // normalized DO-level at the state line objective
 
  vector<lower=0>[2] eta; // Dirichlet hyperparameter 
  real<lower = 1> a; // hyperparameter for the utility function  
  real<lower = 1> p; // L-p norm parameter 
  real z_star[2]; // ideal point for the four objectives 
}


// The parameters accepted by the model
parameters {
  simplex[2] weights; // weights parameter vector 
}



transformed parameters {
  real<lower=0, upper=1.0> prob[n_pairs];    // derived prob. 
  real<lower=0, upper=1.0> u1;               // local variable 
  real<lower=0, upper=1.0> u2;               // local variable 
  
  for (n in 1:n_pairs) {
    u1 = exp(-a * (
      weights[1] * pow(fabs(z_star[1] - obj_1[index_i[n]]), p) +
        weights[2] * pow(fabs(z_star[2] - obj_2[index_i[n]]), p)  
      
    ));
    u2 = exp(-a * (
      weights[1] * pow(fabs(z_star[1] - obj_1[index_j[n]]), p) +
        weights[2] * pow(fabs(z_star[2] - obj_2[index_j[n]]), p) 
       
    ));
    prob[n] = u1 / (u1 + u2);
  }
  
}

// The model to be estimated. We model the output 'y' to be Bernoulli 
// distributed with parameter 'prob'
model {
  weights ~ dirichlet(eta);
  for (n in 1:n_pairs)
    y[n] ~ bernoulli(prob[n]);
}
