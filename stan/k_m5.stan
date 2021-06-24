data{
    int<lower=1> N;
    vector[N] cog;
    vector[N] bf;
    int sex[N];
    int k[N];
    vector[N] age;
    int k_vec;
    
    int edu[N];
    vector[5] alpha; // parameters dirichlet prior 
}
parameters{
    vector[k_vec] a;
    vector[2] s;
    real<lower=0> sigma;
    real bA;
    real bE;
    real bB;
    simplex[5] delta; // will be appended to 0 below
}
model{
    vector[6] delta_j;
    delta ~ dirichlet(alpha);
    delta_j = append_row(0, delta);
    
    vector[N] mu;
    sigma ~ exponential(1);
    a ~ normal(0, 1);
    s ~ normal(0, 1);
    bE ~ normal(0, 0.5);
    bA ~ normal(0, 0.5);
    bB ~ normal(0, 0.5);
    
    for (i in 1:N) {
        mu[i] = a[k[i]] + s[sex[i]] + bE * sum(delta_j[1:edu[i]]) + bA * age[i] + bB * bf[i];
    }
    cog ~ student_t(4, mu, sigma);
}
generated quantities{
    vector[N] log_lik;
    vector[N] mu;
    vector[6] delta_j;
    delta_j = append_row(0, delta);
    for ( i in 1:N ) {
        mu[i] = a[k[i]] + s[sex[i]] + bE * sum(delta_j[1:edu[i]]) + bA * age[i] + bB * bf[i];
    }
    for ( i in 1:N ) log_lik[i] = student_t_lpdf( cog[i] | 4, mu[i] , sigma );
}