data{
    int<lower=1> N_subj;
    vector[N_subj] cog;
    vector[N_subj] shannon;
    vector[N_subj] age;
}
parameters{
    real a_cog;
    real<lower=0> sigma;
    real bS;
    real bA;
}
model{
    vector[N_subj] mu;
    sigma ~ exponential( 1 );
    a_cog ~ normal( 0 , 1 );
    bA ~ normal(0, 0.5);
    bS ~ normal(0, 0.5);
    for ( i in 1:N_subj ) {
        mu[i] = a_cog + shannon[i] * bS + bA * age[i];
    }
    cog ~ student_t( 4, mu , sigma );
}
generated quantities{
    vector[N_subj] log_lik;
    vector[N_subj] mu;
    for ( i in 1:N_subj ) {
        mu[i] = a_cog + shannon[i] * bS + bA * age[i];
    }
    for ( i in 1:N_subj ) log_lik[i] = student_t_lpdf( cog[i] | 4, mu[i] , sigma );
}

