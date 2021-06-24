data{
    int<lower=1> N;
    vector[N] cog;
    vector[N] age;
}
parameters{
    real a;
    real<lower=0> sigma;
    real bA;
}
model{
    vector[N] mu;
    sigma ~ exponential( 1 );
    a ~ normal( 0 , 1 );
    bA ~ normal(0, 0.5);
    
    for ( i in 1:N ) {
        mu[i] = a + bA * age[i];
    }
    cog ~ student_t( 4, mu , sigma );

}
generated quantities{
    vector[N] log_lik;
    vector[N] mu;
    for ( i in 1:N ) {
        mu[i] = a + bA * age[i];
    }
    for ( i in 1:N ) log_lik[i] = student_t_lpdf( cog[i] | 4, mu[i] , sigma );
}

