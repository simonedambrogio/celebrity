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
  int<lower=1,upper=9> cond[N];  // decisions
  
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
  real Intercept_theta_C;
  real Intercept_theta_N;
  real Intercept_theta_S;
  
  real Intercept_bs;  // temporary intercept for centered predictors
  real Intercept_ndt;  // temporary intercept for centered predictors
  vector[K_bias] Intercept_bias;  // temporary intercept for centered predictors
  
  //delta
  vector<lower=0>[M] sd_0;  // group-level standard deviations
  vector[N_sbj] z_0[M];  // standardized group-level effects
  //theta
  vector<lower=0>[M] sd_1_1;  // group-level standard deviations
  vector[N_sbj] z_1_1[M];  // standardized group-level effects
  vector<lower=0>[M] sd_1_2;  // group-level standard deviations
  vector[N_sbj] z_1_2[M];  // standardized group-level effects
  vector<lower=0>[M] sd_1_3;  // group-level standard deviations
  vector[N_sbj] z_1_3[M];  // standardized group-level effects
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
  // delta
  vector[N_sbj] r_1_0;  // actual group-level effects
  // theta
  vector[N_sbj] r_1_1_1;  // actual group-level effects
  vector[N_sbj] r_1_1_2;  // actual group-level effects
  vector[N_sbj] r_1_1_3;  // actual group-level effects
  // bs
  vector[N_sbj] r_2_bs_1;  // actual group-level effects
  // ndt
  vector[N_sbj] r_3_ndt_1;  // actual group-level effects
  // starting-point
  vector[N_sbj] r_4_bias_1;  // actual group-level effects
  
  //delta
  r_1_0 = (sd_0[1] * (z_0[1]));
  //theta
  r_1_1_1 = (sd_1_1[1] * (z_1_1[1]));
  r_1_1_2 = (sd_1_2[1] * (z_1_2[1]));
  r_1_1_3 = (sd_1_3[1] * (z_1_3[1]));
  //bs
  r_2_bs_1 = (sd_2[1] * (z_2[1]));
  //ndt
  r_3_ndt_1 = (sd_3[1] * (z_3[1]));
  //starting-point
  r_4_bias_1 = (sd_4[1] * (z_4[1]));
}
model {
  // likelihood including all constants
  if (!prior_only) {

    vector[N] delta = Intercept_delta + rep_vector(0.0, N);
    vector[N] theta_C = Intercept_theta_C + rep_vector(0.0, N);
    vector[N] theta_N = Intercept_theta_N + rep_vector(0.0, N);
    vector[N] theta_S = Intercept_theta_S + rep_vector(0.0, N);
    vector[N] mu;
    vector[N] theta_L;
    vector[N] theta_R;
    vector[N] bs = Intercept_bs + rep_vector(0.0, N);
    vector[N] ndt = Intercept_ndt + rep_vector(0.0, N);
    vector[N] bias = X_bias * Intercept_bias;
    
    //delta
    for (n in 1:N) {
      delta[n] += r_1_0[J[n]] * Z[n]; // add more terms to the linear predictor
    }
    // theta
    for (n in 1:N) {
      theta_C[n] += r_1_1_1[J[n]] * Z[n]; // add more terms to the linear predictor
    }
    for (n in 1:N) {
      theta_N[n] += r_1_1_2[J[n]] * Z[n]; // add more terms to the linear predictor
    }
    for (n in 1:N) {
      theta_S[n] += r_1_1_3[J[n]] * Z[n]; // add more terms to the linear predictor
    }
    for (n in 1:N) {
      // apply the inverse link function
      theta_C[n] = inv_logit(theta_C[n]);
    }
    for (n in 1:N) {
      // apply the inverse link function
      theta_N[n] = inv_logit(theta_N[n]);
    }
    for (n in 1:N) {
      // apply the inverse link function
      theta_S[n] = inv_logit(theta_S[n]);
    }
    
    
    for (n in 1:N) {
      if( cond[n] == 1){ // Celebrity - Celebrity
        theta_L = theta_C; theta_R = theta_C;
      } else if( cond[n] == 2 ) { // Non Celebrity - Non Celebrity
        theta_L = theta_N; theta_R = theta_N;
      } else if( cond[n] == 3 ){ // Snack - Snack
        theta_L = theta_S; theta_R = theta_S;
      } else if( cond[n] == 4 ){ // Celebrity - Non Celebrity
        theta_L = theta_C; theta_R = theta_N;
      } else if( cond[n] == 5 ){ // Non Celebrity - Celebrity
        theta_L = theta_N; theta_R = theta_C;
      } else if( cond[n] == 6 ){ // Celebrity - Snack
        theta_L = theta_C; theta_R = theta_S;
      } else if( cond[n] == 7 ){ // Snack - Celebrity
        theta_L = theta_S; theta_R = theta_C;
      } else if( cond[n] == 8 ){ // Non Celebrity - Snack
        theta_L = theta_N; theta_R = theta_S;
      } else if( cond[n] == 9 ){ // Snack - Non Celebrity
        theta_L = theta_S; theta_R = theta_N;
      }
      mu[n] = gaze_right[n]*(delta[n]*(value_right[n] - value_left[n]*theta_L[n])) + gaze_left[n]*(delta[n]*(theta_R[n]*value_right[n] - value_left[n]));
    }
    
    //bs
    for (n in 1:N) {
      bs[n] += r_2_bs_1[J[n]] * Z[n]; // add more terms to the linear predictor
    }
    for (n in 1:N) {
      ndt[n] += r_3_ndt_1[J[n]] * Z[n];
    }
    for (n in 1:N) {
      bias[n] += r_4_bias_1[J[n]] * Z[n]; // add more terms to the linear predictor
    }
    
    // apply the inverse link function
    for (n in 1:N) {
      bs[n] = exp(bs[n]);
    }
    for (n in 1:N) {
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
  //delta
  target += std_normal_lpdf(Intercept_delta);
  target += student_t_lpdf(sd_0 | 3, 0, 2.5) - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(z_0[1]);
  //theta
  target += std_normal_lpdf(Intercept_theta_C);
  target += std_normal_lpdf(Intercept_theta_N);
  target += std_normal_lpdf(Intercept_theta_S);
  target += student_t_lpdf(sd_1_1 | 3, 0, 2.5) - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(z_1_1[1]);
  target += student_t_lpdf(sd_1_2 | 3, 0, 2.5) - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(z_1_2[1]);
  target += student_t_lpdf(sd_1_3 | 3, 0, 2.5) - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(z_1_3[1]);
  //bs
  target += normal_lpdf(Intercept_bs | -0.6, 1.3);
  target += student_t_lpdf(sd_2 | 3, 0, 2.5) - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(z_2[1]);
  //ndt
  target += normal_lpdf(Intercept_ndt | -1, 0.5);
  target += student_t_lpdf(sd_3 | 3, 0, 2.5) - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(z_3[1]);
  //bias
  target += std_normal_lpdf(Intercept_bias);
  target += student_t_lpdf(sd_4 | 3, 0, 2.5) - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(z_4[1]);
}
generated quantities {
  vector[N] log_lik;
  vector[N] delta = Intercept_delta + rep_vector(0.0, N);
    vector[N] theta_C = Intercept_theta_C + rep_vector(0.0, N);
    vector[N] theta_N = Intercept_theta_N + rep_vector(0.0, N);
    vector[N] theta_S = Intercept_theta_S + rep_vector(0.0, N);
    vector[N] mu;
    vector[N] theta_L;
    vector[N] theta_R;
    vector[N] bs = Intercept_bs + rep_vector(0.0, N);
    vector[N] ndt = Intercept_ndt + rep_vector(0.0, N);
    vector[N] bias = X_bias * Intercept_bias;
    
    //delta
    for (n in 1:N) {
      delta[n] += r_1_0[J[n]] * Z[n]; // add more terms to the linear predictor
    }
    // theta
    for (n in 1:N) {
      theta_C[n] += r_1_1_1[J[n]] * Z[n]; // add more terms to the linear predictor
    }
    for (n in 1:N) {
      theta_N[n] += r_1_1_2[J[n]] * Z[n]; // add more terms to the linear predictor
    }
    for (n in 1:N) {
      theta_S[n] += r_1_1_3[J[n]] * Z[n]; // add more terms to the linear predictor
    }
    for (n in 1:N) {
      // apply the inverse link function
      theta_C[n] = inv_logit(theta_C[n]);
    }
    for (n in 1:N) {
      // apply the inverse link function
      theta_N[n] = inv_logit(theta_N[n]);
    }
    for (n in 1:N) {
      // apply the inverse link function
      theta_S[n] = inv_logit(theta_S[n]);
    }
    
    
    for (n in 1:N) {
      if( cond[n] == 1){ // Celebrity - Celebrity
        theta_L = theta_C; theta_R = theta_C;
      } else if( cond[n] == 2 ) { // Non Celebrity - Non Celebrity
        theta_L = theta_N; theta_R = theta_N;
      } else if( cond[n] == 3 ){ // Snack - Snack
        theta_L = theta_S; theta_R = theta_S;
      } else if( cond[n] == 4 ){ // Celebrity - Non Celebrity
        theta_L = theta_C; theta_R = theta_N;
      } else if( cond[n] == 5 ){ // Non Celebrity - Celebrity
        theta_L = theta_N; theta_R = theta_C;
      } else if( cond[n] == 6 ){ // Celebrity - Snack
        theta_L = theta_C; theta_R = theta_S;
      } else if( cond[n] == 7 ){ // Snack - Celebrity
        theta_L = theta_S; theta_R = theta_C;
      } else if( cond[n] == 8 ){ // Non Celebrity - Snack
        theta_L = theta_N; theta_R = theta_S;
      } else if( cond[n] == 9 ){ // Snack - Non Celebrity
        theta_L = theta_S; theta_R = theta_N;
      }
      mu[n] = gaze_right[n]*(delta[n]*(value_right[n] - value_left[n]*theta_L[n])) + gaze_left[n]*(delta[n]*(theta_R[n]*value_right[n] - value_left[n]));
    }
    
    //bs
    for (n in 1:N) {
      bs[n] += r_2_bs_1[J[n]] * Z[n]; // add more terms to the linear predictor
    }
    for (n in 1:N) {
      ndt[n] += r_3_ndt_1[J[n]] * Z[n];
    }
    for (n in 1:N) {
      bias[n] += r_4_bias_1[J[n]] * Z[n]; // add more terms to the linear predictor
    }
    
    // apply the inverse link function
    for (n in 1:N) {
      bs[n] = exp(bs[n]);
    }
    for (n in 1:N) {
      ndt[n] = exp(ndt[n]);
    }
    for (n in 1:N) {
      bias[n] = inv_logit(bias[n]);
    }
    
   for (n in 1:N) {
     log_lik[n] = wiener_diffusion_lpdf(Y[n] | dec[n], bs[n], ndt[n], bias[n], mu[n]);
   }
}
