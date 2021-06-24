#####################################################################
##########               README                            ##########
#####################################################################

# In this script I perform 2 steps:
# 1 prior predictive checks:
# The goal is not to find perfect priors for any specific model (as I fit 
# many different models per time point and per covariate structure).
# Rather the goal is to find priors that allow for any possible assocation
# while only being slightly restricting parameter space. Given the amount of
# data we have, the chosen priors will play barely a role as long as I do not 
# prevent the model from searching in realistic parameter space.
# 2 simulate data and see whether our stan model works as expected.



library(tidyverse)
# the rethinking package is convenient for prior predictive simulation and keeps
# the model close to how it is written in stan
library(rethinking)




# 1 Prior predictive checks

# simulate data to fit a model
N <- 150 
shannon <- rnorm(n = N, mean = 0, sd = 1)
ef <- rnorm(n = N, mean = 0.25 * shannon, sd = 1)
d <- list(
  shannon = shannon,
  ef = ef
)
# fit a simple model that includes only slope for shannon
# try different priors as needed
m0 <- ulam(
  alist(
    ef ~ dstudent(4, mu, sigma),
    mu <- a + bS * shannon,
    a ~ dnorm(0, 1),
    bS ~ dnorm(0, 0.5), 
    sigma ~ dexp(1)
  ),
  data = d
)
# plot regression lines only based on the priors:
prior <- extract.prior(m0)
mu0 <- link(m0, post = prior, data = list(shannon = c(-2, 2))) %>% 
  as_tibble() %>%
  pivot_longer(everything(), names_to = "var", values_to = "mu") %>%
  mutate(
    var = ifelse(var == "V1", -2, 2),
    g = rep(1:1000, each = 2))
# plot only the first 200 lines of the posterior to avoid cluttering:
ggplot(head(mu0, 200), aes(var, mu, group = g)) +
  geom_point() +
  geom_path()

                                                                  
# this prior setup fits with our goal of only very mildly restricting parameter 
# space. Many regression lines are still unrealistic given what we expect in 
# terms of effect size. However, the data will easily overwhelm any prior 
# that is not highly restrictive. Therefore, this way, we can apply a similar 
# set of priors to all models. E.g. this will apply to all beta coefficients. 
# Next I will include the ordinal variable in addition and see if that model 
# can recover the parameters used for simulation with these priors and our 
# sample size:

# 2 simulate data

N <- 150
n_sim <- 1:100
posterior_summary <- list()
for (i in seq_along(n_sim)) {
  
  sex <- sample(c(1, 2), size = N, replace = TRUE)
  edu <- round(runif(n = N, 1, 8), 0) 
  age_ef <- rnorm(n = N)
  bf <- rnorm(n = N)
  shannon <- rnorm(n = N, mean = 0.25 * edu - 0.5 * bf, sd = 1)
  # we simulate that sex does not matter, edu, age and bf are all positively 
  # related 
  ef <- rnorm(
    n = N, 
    mean = 0.125 * edu + age_ef * 0.25 + bf * 0.25 + shannon * 0.25, sd = 1
  )
  d <- list(
    edu = edu,
    ef = ef,
    bf = bf,
    age = age_ef,
    sex = sex,
    shannon = shannon,
    alpha = rep(2, 7)
  )

  m1 <- ulam(
    alist(
      ef ~ dstudent(4, mu, sigma),
      mu <- a[sex] + bE * sum(delta_j[1:edu]) + bA * age + bS * shannon + bF *bf,
      a[sex] ~ dnorm(0, 1),
      c(bE, bA, bS, bF) ~ dnorm(0, 0.5), 
      vector[8]: delta_j <<- append_row(0, delta),
      simplex[7]: delta ~ dirichlet(alpha),
      sigma ~ dexp(1)
    ),
    iter = 3000,
    data = d
  )
  post <- extract.samples(m1)
  sex_diff <- as_tibble(post$a, rownames = NA) %>% 
                mutate(diff = V1 - V2) %>%
                summarise(
                  mean = mean(diff), 
                  sd = sd(diff), 
                  lower = quantile(diff, 0.025),
                  upper = quantile(diff, 0.975)
                )  %>%
              mutate(parameter = "sexdiff") %>%
              select(parameter, everything())
  posterior_summary[[i]] <- precis(m1) %>% as.matrix() %>% 
    as.data.frame() %>% 
    rownames_to_column("parameter") %>% 
    select(parameter, mean, sd, lower = "5.5%", upper = "94.5%") %>%
    bind_rows(sex_diff)
}

posterior_summary
# the model reliably correctly estimates sign of the effect. Also effect size 
# estimates are very good.

# same for the model that implements the multilevel model:
library(cmdstanr)

file <- here::here("stan/shannon_brief_m4.stan")
m2 <- cmdstan_model(file)

posterior_summary2 <- list()
n_sim <- 1:100
for (i in seq_along(n_sim)) {
  # since the brief models are normed we leave out sex and age:
  edu <- round(runif(n = N, 1, 6), 0) 
  bf <- rnorm(n = N)
  shannon <- rnorm(n = N, mean = 0.25 * edu + 0.5 * bf, sd = 1)
  mu_ef <- rnorm(
    n = N, 
    mean = 0.2 * edu + bf * 0.25 + shannon * 0.25, sd = 1
  )
  N_subj <- 1:N
  ef <- map2_dfr(mu_ef, N_subj, function(mu, id) {
    tibble(
      subject_id = id, 
      mu_ef = mu,
      # in line with our data we assume here low repeated parent report variation
      ef1 = rnorm(n = 1, mean = mu_ef, sd = 0.4),  
      ef2 = rnorm(n = 1, mean = mu_ef, sd = 0.4) 
      )
    }) 
  cog_obs <- select(ef, subject_id, ef1, ef2) %>%
    pivot_longer(
      all_of(c("ef1", "ef2")), 
      names_to = "time", 
      values_to = "brief")

  datlist <- list(
    shannon = shannon,
    edu = edu,
    bf = bf,
    N_subj = N,
    cog_obs =  cog_obs$brief,
    cog_obs_subj = cog_obs$subject_id,
    N_Y_obs_total = dim(cog_obs)[1],
    alpha = rep(2, 5)
  )

  fit <- m2$sample(
    data = datlist,
    seed = 123,
    chains = 4,
    parallel_chains = 2,
    refresh = 500
  )
  post_sum <- fit$summary() %>% 
                select(variable, mean, q5, q95) %>%
                as_tibble(rownames = NA)
  # how well are the betas estimated in this model?
  betas <- post_sum %>% filter(str_detect(variable, "^b\\w$"))
  # how well does the model identify the true ef based on the two measurements?
  ef_restore <- post_sum %>% filter(str_detect(variable, "^true_cog\\[\\d+\\]$")) %>%
    mutate(subject_id = as.integer(str_extract(variable, "\\d+"))) %>%
    full_join(ef, by = "subject_id") %>% select(subject_id, mu_ef, mean, q5, q95)

  
  posterior_summary2[[i]] <- list(betas = betas, ef = ef_restore)
}

posterior_summary2
# depending on how high we set the SD (ef), the model identifies true mean well.
# Even if we simulate the more difficult case of higher SD (which would NOT be
# in line with our data where correlation is > .8), the model identifies effect
# signs correctly