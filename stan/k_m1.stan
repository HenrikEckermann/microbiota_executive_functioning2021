data{
    int<lower=1> N;
    vector[N] cog;
    int k[N];
    vector[N] age;
    int k_vec;
}
parameters{
    vector[k_vec] a;
    real<lower=0> sigma;
    real bA;
}
model{
    vector[N] mu;
    bA ~ normal(0, 0.5);
    sigma ~ exponential( 1 );
    a ~ normal( 0 , 1 );
    for ( i in 1:N ) {
        mu[i] = a[k[i]] + bA * age[i];
    }
    cog ~ student_t( 4, mu , sigma );
}
generated quantities{
    vector[N] log_lik;
    vector[N] mu;
    for ( i in 1:N ) {
        mu[i] = a[k[i]] + bA * age[i];
    }
    for ( i in 1:N ) log_lik[i] = student_t_lpdf( cog[i] | 4, mu[i] , sigma );
}

