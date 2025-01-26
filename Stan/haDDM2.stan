// generated with brms 2.14.4
functions {
 
   real wiener_diffusion_lpdf(real y, int dec, real alpha, 
                              real tau, real beta, real delta) { 
     if (dec == 1) {
       return wiener_lpdf(y | alpha, tau, beta, delta);
     } else {
       return wiener_lpdf(y | alpha, tau, 1 - beta, - delta);
     }
   }
}

data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=0,upper=1> dec[N];  // decisions
  
  int<lower=1> K_bias;  // number of population-level effects
  matrix[N, K_bias] X_bias;  // population-level design matrix
  
  vector[N] value_right;  // value
  vector[N] value_left;  // value
  real<lower=0,upper=1> gaze_right[N];  // gaze
  real<lower=0,upper=1> gaze_left[N];  // gaze
  
  int<lower=1> N_sbj;  // number of grouping levels
  int<lower=1> M;  // number of coefficients per level
  int<lower=1> J[N];  // grouping indicator per observation
  vector[N] Z; // group-level predictor values
  
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  real min_Y = min(Y);
}
parameters {
  //real Intercept;  // temporary intercept for centered predictors
  real<lower=0> Intercept_delta;
  real Intercept_theta;
  real Intercept_bs;  // temporary intercept for centered predictors
  real Intercept_ndt;  // temporary intercept for centered predictors
  vector[K_bias] Intercept_bias;  // temporary intercept for centered predictors
  
  //delta
  vector<lower=0>[M] sd_0;  // group-level standard deviations
  vector[N_sbj] z_0[M];  // standardized group-level effects
  //theta
  vector<lower=0>[M] sd_1;  // group-level standard deviations
  vector[N_sbj] z_1[M];  // standardized group-level effects
  //boundary-separation
  vector<lower=0>[M] sd_2;  // group-level standard deviations
  vector[N_sbj] z_2[M];  // standardized group-level effects
  //non-decision time
  vector<lower=0>[M] sd_3;  // group-level standard deviations
  vector[N_sbj] z_3[M];  // standardized group-level effects
  //starting point
  vector<lower=0>[M] sd_4;  // group-level standard deviations
  vector[N_sbj] z_4[M];  // standardized group-level effects
}
transformed parameters {
  vector[N_sbj] r_1_0;  // actual group-level effects
  vector[N_sbj] r_1_1;  // actual group-level effects
  vector[N_sbj] r_2_bs_1;  // actual group-level effects
  vector[N_sbj] r_3_ndt_1;  // actual group-level effects
  vector[N_sbj] r_4_bias_1;  // actual group-level effects
  
  r_1_0 = (sd_0[1] * (z_0[1]));
  r_1_1 = (sd_1[1] * (z_1[1]));
  r_2_bs_1 = (sd_2[1] * (z_2[1]));
  r_3_ndt_1 = (sd_3[1] * (z_3[1]));
  r_4_bias_1 = (sd_4[1] * (z_4[1]));
}
model {
  // likelihood including all constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] delta = Intercept_delta + rep_vector(0.0, N);
    vector[N] theta = Intercept_theta + rep_vector(0.0, N);
    vector[N] mu;
    vector[N] bs = Intercept_bs + rep_vector(0.0, N);
    vector[N] ndt = Intercept_ndt + rep_vector(0.0, N);
    vector[N] bias = X_bias * Intercept_bias;
    
    for (n in 1:N) {
      // add more terms to the linear predictor
      delta[n] += r_1_0[J[n]] * Z[n];
    }
    
    for (n in 1:N) {
      // add more terms to the linear predictor
      theta[n] += r_1_1[J[n]] * Z[n];
    }
    
    for (n in 1:N) {
      // apply the inverse link function
      theta[n] = inv_logit(theta[n]);
    }
    
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] = gaze_right[n]*(delta[n]*(value_right[n] - value_left[n]*theta[n])) + gaze_left[n]*(delta[n]*(theta[n]*value_right[n] - value_left[n]));
    }
    
    for (n in 1:N) {
      bs[n] += r_2_bs_1[J[n]] * Z[n];
    }
    for (n in 1:N) {
      ndt[n] += r_3_ndt_1[J[n]] * Z[n];
    }
    for (n in 1:N) {
      bias[n] += r_4_bias_1[J[n]] * Z[n]; // add more terms to the linear predictor
    }
    
    for (n in 1:N) {
      // apply the inverse link function
      bs[n] = exp(bs[n]);
    }
    for (n in 1:N) {
      // apply the inverse link function
      ndt[n] = exp(ndt[n]);
    }
    for (n in 1:N) {
      bias[n] = inv_logit(bias[n]);
    }
    
    for (n in 1:N) {
      target += wiener_diffusion_lpdf(Y[n] | dec[n], bs[n], ndt[n], bias[n], mu[n]);
    }
  }
  // priors including all constants
  target += std_normal_lpdf(Intercept_delta);
  target += std_normal_lpdf(Intercept_theta);
  target += normal_lpdf(Intercept_bs | -0.6, 1.3);
  target += normal_lpdf(Intercept_ndt | -1, 0.5);
  //delta
  target += student_t_lpdf(sd_0 | 3, 0, 2.5) - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(z_0[1]);
  //theta
  target += student_t_lpdf(sd_1 | 3, 0, 2.5) - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(z_1[1]);
  //bs
  target += student_t_lpdf(sd_2 | 3, 0, 2.5) - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(z_2[1]);
  //ndt
  target += student_t_lpdf(sd_3 | 3, 0, 2.5) - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(z_3[1]);
    //bias
  target += std_normal_lpdf(Intercept_bias);
  target += student_t_lpdf(sd_4 | 3, 0, 2.5) - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(z_4[1]);
}
generated quantities {
  vector[N] log_lik;
  // initialize linear predictor term
    vector[N] delta = Intercept_delta + rep_vector(0.0, N);
    vector[N] theta = Intercept_theta + rep_vector(0.0, N);
    vector[N] mu;
    vector[N] bs = Intercept_bs + rep_vector(0.0, N);
    vector[N] ndt = Intercept_ndt + rep_vector(0.0, N);
    vector[N] bias = X_bias * Intercept_bias;
    
    for (n in 1:N) {
      // add more terms to the linear predictor
      delta[n] += r_1_0[J[n]] * Z[n];
    }
    
    for (n in 1:N) {
      // add more terms to the linear predictor
      theta[n] += r_1_1[J[n]] * Z[n];
    }
    
    for (n in 1:N) {
      // apply the inverse link function
      theta[n] = inv_logit(theta[n]);
    }
    
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] = gaze_right[n]*(delta[n]*(value_right[n] - value_left[n]*theta[n])) + gaze_left[n]*(delta[n]*(theta[n]*value_right[n] - value_left[n]));
    }
    
    for (n in 1:N) {
      bs[n] += r_2_bs_1[J[n]] * Z[n];
    }
    for (n in 1:N) {
      ndt[n] += r_3_ndt_1[J[n]] * Z[n];
    }
    for (n in 1:N) {
      bias[n] += r_4_bias_1[J[n]] * Z[n]; // add more terms to the linear predictor
    }
    
    for (n in 1:N) {
      // apply the inverse link function
      bs[n] = exp(bs[n]);
    }
    for (n in 1:N) {
      // apply the inverse link function
      ndt[n] = exp(ndt[n]);
    }
    for (n in 1:N) {
      bias[n] = inv_logit(bias[n]);
    }
   
   for (n in 1:N) {
     log_lik[n] = wiener_diffusion_lpdf(Y[n] | dec[n], bs[n], ndt[n], bias[n], mu[n]);
   }
}
