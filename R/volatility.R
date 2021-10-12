library(tidyverse)
library(phyloseq)
library(vegan)
library(glue)
library(microbiome)
library(cmdstanr)
library(posterior)
library(bayesplot)
standardize <- rethinking::standardize
source(here::here("R/mb_helper.R"))

load(here::here("rdata/bibo_16s.Rds"))
load(here::here("rdata/16s_meta.Rds"))
load(here::here("data/breastfeeding/anybf.Rds"))
df_bf <- dplyr::rename(df_bf, subject_id = id)
meta <- meta %>% full_join(df_bf, by = "subject_id")

# calculate Aitchison distance for the samples in our study 
ids <- sd_to_df(bibo) %>% filter(!is.na(time)) %>% .$sample_id
bibo_clr <- prune_samples(ids, bibo)
bibo_clr <- transform(bibo_clr, transform = "clr")
asv <- otu_to_df(bibo_clr) %>% column_to_rownames("sample_id")
ait <- vegdist(asv, method = "euclidean")
ait <- as.matrix(ait)
# fitler out mock samples 
sids <- sd_to_df(bibo) %>% 
  filter(!is.na(time)) %>% 
  .$subject_id %>% 
  unique()

# calculate volatility sequentially
vol <- map_dfr(sids, function(sid) {
  
  # each subject has maximal 5 samples 
  s1 <- glue("a{sid}")
  s2 <- glue("b{sid}")
  s3 <- glue("c{sid}")
  s4 <- glue("d{sid}")
  s5 <- glue("e{sid}")
  
  # we can only extract volatility out of the ait matrix for non missing pairs 
  vol1 = ifelse(s1 %in% ids & s2 %in% ids, ait[s1, s2], NA)
  vol2 = ifelse(s2 %in% ids & s3 %in% ids, ait[s2, s3], NA)
  vol3 = ifelse(s3 %in% ids & s4 %in% ids, ait[s3, s4], NA)
  vol4 = ifelse(s4 %in% ids & s5 %in% ids, ait[s4, s5], NA)

  
  tibble(
    subject_id = rep(sid, 4),
    time = c("1-2", "2-3","3-4", "4-5"),
    vol = c(vol1, vol2, vol3, vol4)
  )
  
  }) %>% arrange(subject_id, time)



# # doublecheck if this worked
# as.data.frame(ait) %>% rownames_to_column("sid") %>%
#   filter(str_detect(sid, "202")) %>%
#   select(sid, contains("202"))
# vol


vol_by_comp <- group_by(vol, time) %>% nest()

# visualize distributions of volatility per time point pair 
map(vol_by_comp[[2]], function(df) {
  ggplot(df, aes(vol)) +
    geom_density()
})




#####################################################################
##########             Digit  Span Fitting                 ##########
#####################################################################



if (!file.exists(here::here("rdata/volatility_models.Rds"))) {
  outcomes <- c("total_lns", "total_bw", "total_f")
  times <- c("1-2", "2-3", "3-4", "4-5")
  volatility_models <- map(outcomes, function(y) {
    map2(vol_by_comp[[2]], times, function(df, time_var) {
        d <- select(meta, edu, age10, sex, bf, subject_id, all_of(y), time) %>%
          filter(time == 28) %>% # time does not matter here but we need
                                 # one line y
          mutate(
            subject_id = as.character(subject_id),
            edu = as.integer(as.factor(edu))
          ) %>%
          right_join(df, by = "subject_id") %>%
          na.omit() %>%
          mutate(
            subject_id = as.integer(as.factor(subject_id)),
            vol = standardize(vol)
          ) %>%
          select(subject_id, all_of(y), edu, vol, age10, sex, bf) 

      
  
      
        dat_list <- list(
          cog = standardize(d[[y]]),
          vol = d[["vol"]], 
          edu = d[["edu"]],
          alpha = rep(2, 5),
          N_subj = length(d[[y]]),
          age = standardize(d$age10),
          sex = as.integer(as.factor(d$sex)),
          bf = standardize(d$bf)
        )

      
    
      file <- here::here("stan/volatility.stan")
      mod <- cmdstan_model(file)
      fit <- mod$sample(
        data = dat_list,
        seed = 123,
        chains = 4,
        parallel_chains = 2,
        refresh = 500
      )

      fit$summary() %>% filter(rhat > 1.05)
      post <- as_draws_df(fit$draws())
      post_diff <- post %>% mutate(diff = `a_cog[1]` - `a_cog[2]`) %>%
        summarise(
          mean_diff = mean(diff), 
          sd_diff = sd(diff), 
          lower = quantile(diff, 0.025), 
          upper = quantile(diff, 0.975))


      
        list(
          "time" = time_var,
          "y" = y,
          "post" = as.data.frame(post),
          "post_diff" = post_diff
        )
    })
  })
  save(volatility_models, file = here::here("rdata/volatility_models.Rds"))
 } else {
  load(file = here::here("rdata/volatility_models.Rds"))
}



#####################################################################
##########               BRIEF Fitting                     ##########
#####################################################################



if (!file.exists(here::here("rdata/volatility_models_brief.Rds"))) {
    vol_models_brief <- map2(vol_by_comp[[2]], times, function(df, time_var) {
      
      y <- c("brief_total8_t", "brief_total10_t")
      d <- meta %>% filter(time == 28) %>%
        select(
          subject_id,
          all_of(y),
          edu,
          age10,
          bf
        ) %>% 
        mutate(
          subject_id = as.character(subject_id),
          edu = as.integer(as.factor(edu))
        ) %>%
        left_join(df, by = "subject_id") %>%
        na.omit() %>%  
        mutate(subject_id = as.integer(as.factor(subject_id))) %>%
        arrange(subject_id) 



      cog_lm <- d %>% select(subject_id, all_of(y)) %>%
        pivot_longer(
          all_of(y), 
          names_to = "time", 
          values_to = "brief") %>%
          mutate(brief = standardize(brief))


      
      dat_list <- list(
        N_subj = length(d$subject_id),
        N_Y_obs_total = dim(cog_lm)[1],
        cog_obs = cog_lm$brief,
        cog_obs_subj = cog_lm$subject_id,
        edu = d[["edu"]],
        alpha = rep(2, 5),
        vol = standardize(d[["vol"]]),
        bf = standardize(d$bf)
      )


      file <- here::here("stan/volatility_brief.stan")
      mod <- cmdstan_model(file)

      fit <- mod$sample(
        data = dat_list,
        seed = 123,
        chains = 4,
        parallel_chains = 2,
        refresh = 500
      )

      fit$summary() %>% filter(rhat > 1.05)
      post_sum <-fit$summary() %>% select(variable, mean, q5, q95)
      post <- as_draws_df(fit$draws())

      list(
        "time" = time_var,
        "post_sums" = post_sum,
        "post" = as.data.frame(post)
      )
    })

  save(vol_models_brief, file = here::here("rdata/volatility_models_brief.Rds"))
 } else {
  load(file = here::here("rdata/volatility_models_brief.Rds"))
}






#####################################################################
##########                Plotting                         ##########
#####################################################################


# extract betas from DS models 
b_post <- map_dfr(volatility_models, function(listobj) {
  map_dfr(listobj, function(listobj2) {
    listobj2$post %>% as_tibble() %>%
      select(volatility = bV) %>%
      mutate(time = listobj2$time, y = listobj2$y)
  })
})

b_post_grouped <- b_post %>% group_by(time, y) %>%
  nest()

# reshape 
b_post1 <- bind_cols(b_post_grouped$data[c(1:4)])
b_post2 <- bind_cols(b_post_grouped$data[c(5:8)])
b_post3 <- bind_cols(b_post_grouped$data[c(9:12)])

b_post <- map(list(b_post1, b_post2, b_post3), function(post) {
  colnames(post) <- glue("Time {unique(b_post_grouped$time)}")
  post
})


# extract slopes for BRIEF models 
vol_post <- map_dfr(vol_models_brief, function(listobj) {
    listobj$post %>% as_tibble() %>%
      select(volatility = bV) %>%
      mutate(time = listobj$time, y = listobj$y)}) %>%
    group_by(time) %>%
    nest()
b_post4 <- bind_cols(vol_post$data)
colnames(b_post4) <- colnames(b_post[[1]])
b_post[[4]] <- b_post4






# plots betas 
b_plots_vol <- map2(b_post, list("DS Total LNS", "DS Total Backwards", "DS Total Forwards", "BRIEF"),  function(post, y) {
  color_scheme_set("gray")
  post %>%
    mcmc_areas(
      pars = vars(contains("Time ")),
      prob = 0.95) +
      ggtitle(y) +
      geom_vline(aes(xintercept = 0), linetype = "dashed") +
      theme_bw(base_size = 20)
})


b_plots_vol

save(b_plots_vol, file = here::here("rdata/b_plots_vol.Rds"))

select(meta, contains("total")) %>%
  GGally::ggpairs()
  