data{
    int<lower=1> N_subj;
    vector[N_subj] cog;
    vector[N_subj] shannon;
    vector[N_subj] age;    
    int edu[N_subj];
    vector[5] alpha; // parameters dirichlet prior 
}
parameters{
    real a_cog;
    real<lower=0> sigma;
    real bS;
    real bA;
    real bE;
    simplex[5] delta; // will be appended to 0 below
}
model{
    vector[6] delta_j;
    delta ~ dirichlet(alpha);
    delta_j = append_row(0, delta);
    
    vector[N_subj] mu;
    sigma ~ exponential(1);
    a_cog ~ normal(0, 1);
    bE ~ normal(0, 0.5);
    bS ~ normal(0, 0.5);
    bA ~ normal(0, 0.5);
    for (i in 1:N_subj) {
        mu[i] = a_cog + bS * shannon[i] + bE * sum(delta_j[1:edu[i]]) + bA * age[i];
    }
    cog ~ student_t(4, mu, sigma);
}
generated quantities{
    vector[N_subj] log_lik;
    vector[N_subj] mu;
    vector[6] delta_j;
    delta_j = append_row(0, delta);
    for ( i in 1:N_subj ) {
        mu[i] = a_cog + bS * shannon[i] + bE * sum(delta_j[1:edu[i]]) + bA * age[i];
    }
    for ( i in 1:N_subj ) log_lik[i] = student_t_lpdf( cog[i] | 4, mu[i] , sigma );
}