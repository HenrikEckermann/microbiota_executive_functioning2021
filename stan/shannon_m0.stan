data{
    int<lower=1> N_subj;
    vector[N_subj] cog;
    vector[N_subj] age;
}
parameters{
    real a_cog;
    real bA;
    real<lower=0> sigma;
}
model{
    vector[N_subj] mu;
    sigma ~ exponential( 1 );
    bA ~ normal(0, 0.5);
    a_cog ~ normal( 0 , 1 );
    for (j in 1:N_subj) {
      mu[j] = a_cog + bA * age[j];
    }
    cog ~ student_t( 4, mu , sigma );
}
generated quantities{
    vector[N_subj] log_lik;
    vector[N_subj] mu;
    for ( i in 1:N_subj ) {
        mu[i] = a_cog + bA * age[i];
    }
    for ( i in 1:N_subj ) log_lik[i] = student_t_lpdf( cog[i] | 4, mu[i] , sigma );
}

