data{
    int<lower=1> N_subj;
    int<lower=1> N_Y_obs_total;
    vector[N_Y_obs_total] cog_obs;
    vector[N_subj] bf;
    vector[N_subj] shannon;
    int<lower=1, upper=N_subj> cog_obs_subj[N_Y_obs_total];  
    
    int edu[N_subj];
    vector[5] alpha; // parameters dirichlet prior 
}
parameters{
    vector[N_subj] true_cog;
    real<lower=0> sigma_cog;
    real a_cog;
    real<lower=0> sigma;
    real bB;
    real bS;
    real bE;
    simplex[5] delta; // will be appended to 0 below
}
model{
    vector[6] delta_j;
    delta ~ dirichlet( alpha );
    delta_j = append_row(0, delta);
    
    vector[N_subj] mu;
    sigma ~ exponential( 1 );
    bB ~ normal(0, 0.5);
    bS ~ normal(0, 0.5);
    bE ~ normal(0, 0.5);
    sigma_cog ~ exponential(1);
    a_cog ~ normal(0, 1);
    true_cog ~ normal(0, 5);

    
    
    
    for (i in 1:N_Y_obs_total) {
      cog_obs[i] ~ normal(true_cog[cog_obs_subj[i]], sigma_cog);
    }
    
    for (j in 1:N_subj) {
      mu[j] = a_cog + bE * sum(delta_j[1:edu[j]]) + bB * bf[j] + bS * shannon[j];
    }
    
    true_cog ~ student_t(4, mu, sigma);
}
generated quantities{
    vector[N_subj] log_lik;
    vector[N_subj] mu;
    vector[6] delta_j;
    delta_j = append_row(0, delta);
    for ( j in 1:N_subj ) {
        mu[j] = a_cog + bE * sum(delta_j[1:edu[j]]) + bB * bf[j] + bS * shannon[j];
    }
    for ( j in 1:N_subj ) log_lik[j] = student_t_lpdf( true_cog[j] | 4, mu[j] , sigma );
}





