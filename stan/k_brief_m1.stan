data{
    int<lower=1> N_subj;
    int<lower=1> N_Y_obs_total;
    int k[N_subj];
    int k_vec;
    vector[N_Y_obs_total] cog_obs;
    int<lower=1, upper=N_subj> cog_obs_subj[N_Y_obs_total];  
}
parameters{
    vector[N_subj] true_cog;
    real<lower=0> sigma_cog;
    vector[k_vec] a_cog;
    real<lower=0> sigma;
}
model{
    vector[N_subj] mu;
    sigma ~ exponential( 1 );

    sigma_cog ~ exponential(1);
    a_cog ~ normal(0, 1);
    true_cog ~ normal(0, 5);
    
    
    
    for (i in 1:N_Y_obs_total) {
      cog_obs[i] ~ normal(true_cog[cog_obs_subj[i]], sigma_cog);
    }
    
    for (j in 1:N_subj) {
      mu[j] = a_cog[k[j]];
    }
    
    true_cog ~ student_t(4, mu, sigma);
}
generated quantities{
    vector[N_subj] log_lik;
    vector[N_subj] mu;
    for ( j in 1:N_subj ) {
        mu[j] = a_cog[k[j]];
    }
    for ( j in 1:N_subj ) log_lik[j] = student_t_lpdf( true_cog[j] | 4, mu[j] , sigma );
}





