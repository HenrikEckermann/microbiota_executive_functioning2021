library(tidyverse)
library(patchwork)
library(phyloseq)
library(microbiome)
library(tidybayes)
library(cluster)
library(fpc)
library(DirichletMultinomial)
library(rethinking)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(loo)
map <- purrr::map

# helper functions for pseq objects 
source(here::here("R/mb_helper.R"))
# helper function for RF analyses 
source(here::here("R/ml_helper.R"))

# load data 
load(file = here::here("rdata/16s_meta.rds"))
load(here::here("rdata/bibo_16s.Rds"))
# agglomerate at genus level 
pseq <- tax_glom(bibo, taxrank = "Genus")
# for the RF analyses we need relative abundance 
pseq_rel <- microbiome::transform(pseq, transform = "compositional")
# for the PAM method we need to use clr abundance 
pseq_clr <- microbiome::transform(pseq, transform = "clr")















#####################################################################
##########             1 Random Forest analyses            ##########
#####################################################################

# this objects contains per timepoint which sample_ids belong to that timepoint
identifier_by_time <- sample_data(bibo) %>% 
  sd_to_df() %>%
  select(sample_id, subject_id, time) %>%
  filter(!is.na(time)) %>% 
  group_by(time) %>% 
  nest()

# use the mean of brief 
meta <- mutate(meta, brief = (brief_total8_t + brief_total10_t) / 2)
# the analyses will be performed per outcome Y and per timepoint
# therefore, I create an object that eases mapping the analyses workflow later
Y <- c("total_lns", "total_f", "total_bw", "brief")

data_per_y <- map(Y, function(y) {
  # extract complete data according to y
  otus <- map2(identifier_by_time$time, identifier_by_time$data, function(time_var, df_temp) {
    
    otu <- pseq_rel %>% 
      otu_to_df() %>% 
      filter(sample_id %in% df_temp$sample_id) %>%  
      full_join(df_temp, otu, by = "sample_id")
      
    df_y <- meta %>% select(all_of(y), subject_id, time) %>%
      filter(time == time_var) %>% 
      select(-time) %>%
      na.omit() %>%
      mutate(subject_id = as.character(subject_id))
    
    # listwise deletion!
    combined <- full_join(
      otu, df_y, 
      by = "subject_id") %>% 
        na.omit() %>%
        arrange(subject_id) %>%
        mutate(time = time_var)
    combined
  })
})

# now we can run the analyses using the helper functions from the loaded script.
# This analyses ran on a compute cluster of the Donders Institute for Brain, 
# Cognition and Behavior. Running the algorithm took several days. Make sure to
# adjust parameters depending on your goal: 

model_sums <- map2(Y, data_per_y, function(y, data_per_time) {
  map(data_per_time, function(data) {
    if (!file.exists(here::here(glue("rdata/rf_{data$time[1]}_{y}_genus.Rds")))) {
      # select only bacterial genera
      X <- data %>% 
        select(matches("\\d+")) %>% 
        colnames()
      # 10x 4 fold CV without any tuning
      model_and_data0 <- future_map(1:10, function(rep) {
        rf_cv(
          data = data,
          features = X,
          y = y,
          k = 4,
          ntree = 5000
        )
      }, .options = furrr_options(seed = TRUE)) %>% flatten()
        
      # hyperparameter tuning
      pars <- tune_rf(
        data = data,
        features = X,
        y = y,
        regression = TRUE,
        ntree = 5000,
        tune.parameters = c("mtry", "sample.fraction")
      )
        
        
      # 10x 4 fold CV with tuning 
      model_and_data1 <- future_map(1:10, function(rep) {
        rf_cv(
          data = data,
          features = X,
          y = y,
          k = 4,
          ntree = 5000,
          regression = TRUE,
          mtry = pars$recommended.pars[1, "mtry"],
          sample.fraction = pars$recommended.pars[1, "sample.fraction"],
        )}, .options = furrr_options(seed = TRUE)) %>% flatten()

        # we obtain out of bag error and pearson correlation of pred and actual 
        # values of y
        oob0 = get_oob(model_and_data0)
        oob1 = get_oob(model_and_data1)
        pearson0 = get_pearson(model_and_data0, y)
        pearson1 = get_pearson(model_and_data1, y)
        
        save(oob0, oob1, pearson0, pearson1, pars, file = here::here(glue("rdata/rf_{data$time[1]}_{y}_genus.Rds")))
        
      } else {
        load(here::here(glue("rdata/rf_{data$time[1]}_{y}_genus.Rds")))
      }
    
    return(list(
      "time" = data$time[1],
      "Y" = y,
      "pars" = pars,
      "oob0" = oob0,
      "oob1" = oob1,
      "pearson0" = pearson0,
      "pearson1" = pearson1
    ))
  })
})

# summarise results to later use the median correlation for p value calculation
# note that pearson0 will be ignored and was solely to compare how tuning 
# performs without tuning.

results_df <- map_dfr(model_sums, function(models_per_y) {
  map_dfr(models_per_y, function(mlist) {
    bind_rows(
      mlist$oob0 %>% 
        mutate(statistic = "oob0"), 
      mlist$oob1 %>% mutate(statistic = "oob1"), 
      mlist$pearson0 %>% 
        mutate(statistic = "pearson0"), 
      mlist$pearson1 %>% 
        mutate(statistic = "pearson1")) %>% 
        mutate(y = mlist$Y, time = mlist$time)
  })
})

# obtain p value for nperm number of permutations
nperm <- 1000
nulldist <- map_dfr(1:nperm, function(nulliter) {
  if (!file.exists(here::here(glue("rdata/perm{nulliter}.Rds")))) {
    rf_null <- map2_dfr(Y, data_per_y, function(y, data_per_time) {
      map(data_per_time, function(data) {
        
        # replace y with the permuted y, the rest stays the same.
        y_perm <- sample(data[[y]], size = length(data[[y]]), replace = FALSE)
        data[[y]] <- y_perm
        
        X <- data %>% 
          select(matches("\\d+")) %>% 
          colnames()
        
        pars <- tune_rf(
          data = data,
          features = X,
          y = y,
          regression = TRUE,
          ntree = 5000,
          tune.parameters = c("mtry", "sample.fraction")
        )
          
          
        model_and_data <- future_map(1:10, function(rep) {
          rf_cv(
            data = data,
            features = X,
            y = y,
            k = 4,
            ntree = 5000,
            regression = TRUE,
            mtry = pars$recommended.pars[1, "mtry"],
            sample.fraction = pars$recommended.pars[1, "sample.fraction"],
          )}, .options = furrr_options(seed = TRUE)) %>% flatten()
        
        oob = get_oob(model_and_data) %>% mutate(statistic = "oob")
        pearson = get_pearson(model_and_data, y) %>% mutate(statistic = "pearson")
        
        bind_rows(oob, pearson) %>%
          mutate(y = y, time = data$time[1])
      })
    })
    
    save(rf_null, file = here::here(glue("rdata/perm{nulliter}.Rds")))
  } else {
      load(here::here(glue("rdata/perm{nulliter}.Rds")))
    }
    
    rf_null
})
# calculate p value per time and outcome 
null_per_test <- nulldist %>% group_by(time, y, statistic) %>% nest()
alt_per_test <- results_df %>% 
  filter(statistic %in% c("oob1", "pearson1")) %>% # only use the tuned models
  group_by(time, y, statistic) %>% 
  nest()
alt_per_test$data <- map2(null_per_test$data, alt_per_test$data, function(null, alt) {
  
  alt$pvalue <- mean(alt$median < null$median)
  alt
})
# correct for multiple testing 
summary_df <- alt_per_test %>% 
  unnest(cols = c(data)) %>%
  ungroup()
summary_df$q  <- qvalue::qvalue(summary_df$pvalue)$qvalue
summary_df %>% 
  filter(statistic == "pearson1") %>%
  select(time, y, median, pvalue, q, -statistic) %>%
  arrange(time, y)

# summary_df is table 3 in the paper
save(summary_df, file = here::here("rdata/rf_summary_df.Rds"))


# we also tried above algorithm with 2 different feature selection methods:
# 1st method: obtain the core microbiota before performing above procedure:

pseq_infant <- subset_samples(pseq_rel, time %in% c(28, 75, 105))
pseq_adult <- subset_samples(pseq_rel, !time %in% c(28, 75, 105))
core_infancy <- core_members(pseq_infant, detection = 1/100, prevalence = 10/100)
core_adult <- core_members(pseq_adult, detection = 1/100, prevalence = 10/100)

model_sums_core <- map2(Y, data_per_y, function(y, data_per_time) {
  map(data_per_time, function(data) {
    if (!file.exists(here::here(glue("rdata/rf_{data$time[1]}_{y}_genus_core.Rds")))) {
      # select only core genera
      if (data$time[1] %in% c(28, 75, 105)) {
        X <- core_infancy
      } else if (data$time[1] %in% c(2193, 3655)) {
        X <- core_adult
      }
        
      # hyperparameter tuning
      pars <- tune_rf(
        data = data,
        features = X,
        y = y,
        regression = TRUE,
        ntree = 5000,
        tune.parameters = c("mtry", "sample.fraction")
      )
        
        
      # 10x 4 fold CV with tuning 
      model_and_data <- future_map(1:4, function(rep) {
        rf_cv(
          data = data,
          features = X,
          y = y,
          k = 4,
          ntree = 5000,
          regression = TRUE,
          mtry = pars$recommended.pars[1, "mtry"],
          sample.fraction = pars$recommended.pars[1, "sample.fraction"],
        )}, .options = furrr_options(seed = TRUE)) %>% flatten()

        # we obtain out of bag error and pearson correlation of pred and actual 
        # values of y
        oob = get_oob(model_and_data)
        pearson = get_pearson(model_and_data, y)
        
        save(oob, pearson, pars, file = here::here(glue("rdata/rf_{data$time[1]}_{y}_genus_core.Rds")))
        
      } else {
        load(here::here(glue("rdata/rf_{data$time[1]}_{y}_genus_core.Rds")))
      }
    
    return(list(
      "time" = data$time[1],
      "Y" = y,
      "pars" = pars,
      "oob" = oob,
      "pearson" = pearson
    ))
  })
})

# summarise results to later use the median correlation for p value calculation
# note that pearson0 will be ignored and was solely to compare how tuning 
# performs without tuning.

results_df_core <- map_dfr(model_sums_core, function(models_per_y) {
  map_dfr(models_per_y, function(mlist) {
    bind_rows(
      mlist$oob %>% mutate(statistic = "oob"), 
      mlist$pearson %>% 
        mutate(statistic = "pearson")) %>% 
        mutate(y = mlist$Y, time = mlist$time)
  })
})


# obtain p value for nperm number of permutations
nperm <- 1000
nulldist_core <- map_dfr(1:nperm, function(nulliter) {
  if (!file.exists(here::here(glue("rdata/perm{nulliter}_core.Rds")))) {
    rf_null <- map2_dfr(Y, data_per_y, function(y, data_per_time) {
      map(data_per_time, function(data) {
        
        # replace y with the permuted y, the rest stays the same.
        y_perm <- sample(data[[y]], size = length(data[[y]]), replace = FALSE)
        data[[y]] <- y_perm
        
        # select only core genera
        if (data$time[1] %in% c(28, 75, 105)) {
          X <- core_infancy
        } else if (data$time[1] %in% c(2193, 3655)) {
          X <- core_adult
        }
        
        pars <- tune_rf(
          data = data,
          features = X,
          y = y,
          regression = TRUE,
          ntree = 5000,
          tune.parameters = c("mtry", "sample.fraction")
        )
          
          
        model_and_data <- future_map(1:4, function(rep) {
          rf_cv(
            data = data,
            features = X,
            y = y,
            k = 4,
            ntree = 5000,
            regression = TRUE,
            mtry = pars$recommended.pars[1, "mtry"],
            sample.fraction = pars$recommended.pars[1, "sample.fraction"],
          )}, .options = furrr_options(seed = TRUE)) %>% flatten()
        
        oob = get_oob(model_and_data) %>% mutate(statistic = "oob")
        pearson = get_pearson(model_and_data, y) %>% mutate(statistic = "pearson")
        
        bind_rows(oob, pearson) %>%
          mutate(y = y, time = data$time[1])
      })
    })
    
    save(rf_null, file = here::here(glue("rdata/perm{nulliter}_core.Rds")))
  } else {
      load(here::here(glue("rdata/perm{nulliter}_core.Rds")))
    }
    
    rf_null
})


# calculate p value per time and outcome 
null_per_test_core <- nulldist_core %>% group_by(time, y, statistic) %>% nest()
alt_per_test_core <- results_df_core %>% 
  filter(statistic %in% c("oob", "pearson")) %>% # only use the tuned models
  group_by(time, y, statistic) %>% 
  nest()
alt_per_test_core$data <- map2(null_per_test_core$data, alt_per_test_core$data, function(null, alt) {
  
  alt$pvalue <- mean(alt$median < null$median)
  alt
})

# correct for multiple testing 
summary_df_core <- alt_per_test_core %>% 
  unnest(cols = c(data)) %>%
  ungroup()
summary_df_core
summary_df_core$q  <- qvalue::qvalue(summary_df_core$pvalue)$qvalue
summary_df_core %>% 
  filter(statistic == "pearson") %>%
  select(time, y, median, pvalue, q, -statistic) %>%
  arrange(time, y)

save(summary_df_core, file = here::here("rdata/rf_summary_df_core.Rds"))





# now using the feature selection algorithm 

model_sums_fl <- map2(Y, data_per_y, function(y, data_per_time) {
  map(data_per_time, function(data) {
    if (!file.exists(here::here(glue("rdata/rf_{data$time[1]}_{y}_genus_fl.Rds")))) {
        
      X <- data %>% 
        select(matches("\\d+")) %>% 
        colnames()
          
      model_and_data_fl <- rf_cv(
        data,
        features = X,
        y = y,
        k = 4,
        ntree = 5000,
        regression = TRUE
      )
        
      selected_feat <- select_features(
        model_and_data_fl,
        id_name = "genera",
        n_features = 40 # we tried 30, 40 and 50
      )
      
      # hyperparameter tuning
      pars <- tune_rf(
        data = data,
        features = selected_feat[["genera"]],
        y = y,
        regression = TRUE,
        ntree = 5000,
        tune.parameters = c("mtry", "sample.fraction")
      )
        
      # 10x 4 fold CV with tuning 
      model_and_data <- future_map(1:4, function(rep) {
        rf_cv(
          data = data,
          features = selected_feat[["genera"]],
          y = y,
          k = 4,
          ntree = 5000,
          regression = TRUE,
          mtry = pars$recommended.pars[1, "mtry"],
          sample.fraction = pars$recommended.pars[1, "sample.fraction"],
        )}, .options = furrr_options(seed = TRUE)) %>% flatten()
      

        # we obtain out of bag error and pearson correlation of pred and actual 
        # values of y
        oob = get_oob(model_and_data)
        pearson = get_pearson(model_and_data, y)
        
        save(oob, pearson, pars, file = here::here(glue("rdata/rf_{data$time[1]}_{y}_genus_fl.Rds")))
        
      } else {
        load(here::here(glue("rdata/rf_{data$time[1]}_{y}_genus_fl.Rds")))
      }
    
    return(list(
      "time" = data$time[1],
      "Y" = y,
      "pars" = pars,
      "oob" = oob,
      "pearson" = pearson
    ))
  })
})

# summarise results to later use the median correlation for p value calculation
# note that pearson0 will be ignored and was solely to compare how tuning 
# performs without tuning.

results_df_fl <- map_dfr(model_sums_fl, function(models_per_y) {
  map_dfr(models_per_y, function(mlist) {
    bind_rows(
      mlist$oob %>% mutate(statistic = "oob"), 
      mlist$pearson %>% 
        mutate(statistic = "pearson")) %>% 
        mutate(y = mlist$Y, time = mlist$time)
  })
})

# obtain p value for nperm number of permutations
nperm <- 1000
nulldist_fl <- map_dfr(1:nperm, function(nulliter) {
  if (!file.exists(here::here(glue("rdata/perm{nulliter}_fl.Rds")))) {
    rf_null <- map2_dfr(Y, data_per_y, function(y, data_per_time) {
      map(data_per_time, function(data) {
        
        # replace y with the permuted y, the rest stays the same.
        y_perm <- sample(data[[y]], size = length(data[[y]]), replace = FALSE)
        data[[y]] <- y_perm
        
        X <- data %>% 
          select(matches("\\d+")) %>% 
          colnames()
          
        model_and_data_fl <- rf_cv(
          data,
          features = X,
          y = y,
          k = 4,
          ntree = 5000,
          regression = TRUE
        )
          
        selected_feat <- select_features(
          model_and_data_fl,
          id_name = "genera",
          n_features = 40 # we tried 30, 40 and 50
        )
        
        pars <- tune_rf(
          data = data,
          features = selected_feat[["genera"]],
          y = y,
          regression = TRUE,
          ntree = 5000,
          tune.parameters = c("mtry", "sample.fraction")
        )
          
          
        model_and_data <- future_map(1:4, function(rep) {
          rf_cv(
            data = data,
            features = selected_feat[["genera"]],
            y = y,
            k = 4,
            ntree = 5000,
            regression = TRUE,
            mtry = pars$recommended.pars[1, "mtry"],
            sample.fraction = pars$recommended.pars[1, "sample.fraction"],
          )}, .options = furrr_options(seed = TRUE)) %>% flatten()
        
        oob = get_oob(model_and_data) %>% mutate(statistic = "oob")
        pearson = get_pearson(model_and_data, y) %>% mutate(statistic = "pearson")
        
        bind_rows(oob, pearson) %>%
          mutate(y = y, time = data$time[1])
      })
    })
    
    save(rf_null, file = here::here(glue("rdata/perm{nulliter}_fl.Rds")))
  } else {
      load(here::here(glue("rdata/perm{nulliter}_fl.Rds")))
    }
    
    rf_null
})
# calculate p value per time and outcome 
null_per_test_fl <- nulldist_fl %>% group_by(time, y, statistic) %>% nest()
alt_per_test_fl <- results_df_fl %>% 
  filter(statistic %in% c("oob", "pearson")) %>% # only use the tuned models
  group_by(time, y, statistic) %>% 
  nest()
alt_per_test_fl$data <- map2(null_per_test_fl$data, alt_per_test_fl$data, function(null, alt) {
  
  alt$pvalue <- mean(alt$median < null$median)
  alt
})
# correct for multiple testing 
summary_df_fl <- alt_per_test_fl %>% 
  unnest(cols = c(data)) %>%
  ungroup()
summary_df_fl$q  <- qvalue::qvalue(summary_df_fl$pvalue)$qvalue
summary_df_fl %>% 
  filter(statistic == "pearson") %>%
  select(time, y, median, pvalue, q, -statistic) %>%
  arrange(time, y)

save(summary_df_fl, file = here::here("rdata/rf_summary_df_fl.Rds"))






















#####################################################################
##########   2 linear models ef ~ shannon + covariates     ##########
#####################################################################


if (!file.exists(here::here("rdata/16s_shannon_models.Rds"))) {
  outcomes <- c("total_lns", "total_bw", "total_f")
  times <- c(28, 75, 105, 2193, 3655)
  shannon_models <- map(outcomes, function(y) {
    map(times, function(time_var) {
      # first create a dataset corresponding to the time point
      # after joining the dataframes, I prepare them for the input to stan
        ad_temp <- select(meta, subject_id, diversity_shannon, time) %>%
          filter(time == time_var)
        df_lm <- meta %>% filter(time == time_var) %>%
          select(
            subject_id,
            all_of(y),
            edu,
            age10,
            sex,
            bf
          ) %>% 
          mutate(
            subject_id = as.character(subject_id),
            edu = as.integer(as.factor(edu))
          ) %>%
          left_join(ad_temp, by = "subject_id") %>%
          select(-time) %>%
          na.omit() %>% # complete case analysis
          mutate(
            subject_id = as.integer(as.factor(subject_id)),
            shannon = standardize(diversity_shannon)
          ) %>%
          select(subject_id, all_of(y), edu, shannon, age10, sex, bf) %>% 
          arrange(subject_id) 

        # we explore different covariate structures, each requiring different 
        # datasets
        dat_list0 <- list(
          cog = standardize(df_lm[[y]]),
          N_subj = length(df_lm[[y]]),
          age = standardize(df_lm$age10)
        )

        dat_list1 <- list(
          cog = standardize(df_lm[[y]]),
          shannon = df_lm[["shannon"]],
          N_subj = length(df_lm[[y]]),
          age = standardize(df_lm$age10)
        )

        dat_list2 <- list(
          cog = standardize(df_lm[[y]]),
          shannon = df_lm[["shannon"]],
          edu = df_lm[["edu"]],
          alpha = rep(2, 5), # see stan model 
          N_subj = length(df_lm[[y]]),
          age = standardize(df_lm$age10)
        )
        
        dat_list3 <- list(
          cog = standardize(df_lm[[y]]),
          edu = df_lm[["edu"]],
          alpha = rep(2, 5),
          N_subj = length(df_lm[[y]]),
          age = standardize(df_lm$age10)
        )
        
        dat_list4 <- list(
          cog = standardize(df_lm[[y]]),
          shannon = df_lm[["shannon"]],
          edu = df_lm[["edu"]],
          alpha = rep(2, 5),
          N_subj = length(df_lm[[y]]),
          age = standardize(df_lm$age10),
          sex = as.integer(as.factor(df_lm$sex))
        )
        
        dat_list5 <- list(
          cog = standardize(df_lm[[y]]),
          shannon = df_lm[["shannon"]],
          edu = df_lm[["edu"]],
          alpha = rep(2, 5),
          N_subj = length(df_lm[[y]]),
          age = standardize(df_lm$age10),
          sex = as.integer(as.factor(df_lm$sex)),
          bf = standardize(df_lm$bf)
        )

        # for each covariate structure, we use a different stan model that we 
        # we will fit to then extract the posterior distribution for summary
        file0 <- here::here("stan/shannon_m0.stan")
        mod0 <- cmdstan_model(file0)
        fit0 <- mod0$sample(
          data = dat_list0,
          seed = 123,
          chains = 4,
          parallel_chains = 2,
          refresh = 500
        )
        fit0$summary() %>% filter(rhat > 1.05)
        post_sum0 <-fit0$summary() %>% select(variable, mean, q5, q95)
        loo0 <- loo(fit0$draws("log_lik"))



        file1 <- here::here("stan/shannon_m1.stan")
        mod1 <- cmdstan_model(file1)
        fit1 <- mod1$sample(
          data = dat_list1,
          seed = 123,
          chains = 4,
          parallel_chains = 2,
          refresh = 500
        )

        fit1$summary() %>% filter(rhat > 1.05)
        post_sum1 <-fit1$summary() %>% select(variable, mean, q5, q95)
        loo1 <- loo(fit1$draws("log_lik"))



        file2 <- here::here("stan/shannon_m2.stan")
        mod2 <- cmdstan_model(file2)
        fit2 <- mod2$sample(
          data = dat_list2,
          seed = 123,
          chains = 4,
          parallel_chains = 2,
          refresh = 500
        )

        fit2$summary() %>% filter(rhat > 1.05)
        post_sum2 <-fit2$summary() %>% select(variable, mean, q5, q95)
        loo2 <- loo(fit2$draws("log_lik"))
        
        
        file3 <- here::here("stan/shannon_m3.stan")
        mod3 <- cmdstan_model(file3)
        fit3 <- mod3$sample(
          data = dat_list3,
          seed = 123,
          chains = 4,
          parallel_chains = 2,
          refresh = 500
        )

        fit3$summary() %>% filter(rhat > 1.05)
        post_sum3 <-fit3$summary() %>% select(variable, mean, q5, q95)
        loo3 <- loo(fit3$draws("log_lik"))
        
        file4 <- here::here("stan/shannon_m4.stan")
        mod4 <- cmdstan_model(file4)
        fit4 <- mod4$sample(
          data = dat_list4,
          seed = 123,
          chains = 4,
          parallel_chains = 2,
          refresh = 500
        )

        fit4$summary() %>% filter(rhat > 1.05)
        post_sum4 <-fit4$summary() %>% select(variable, mean, q5, q95)
        
      post4 <- as_draws_df(fit4$draws())
      post4_diff <- post4 %>% mutate(diff = `a_cog[1]` - `a_cog[2]`) %>%
        summarise(
          mean_diff = mean(diff), 
          sd_diff = sd(diff), 
          lower = quantile(diff, 0.025), 
          upper = quantile(diff, 0.975))
        
        loo4 <- loo(fit4$draws("log_lik"))
        
        
        file5 <- here::here("stan/shannon_m5.stan")
        mod5 <- cmdstan_model(file5)
        fit5 <- mod5$sample(
          data = dat_list5,
          seed = 123,
          chains = 4,
          parallel_chains = 2,
          refresh = 500
        )

        fit5$summary() %>% filter(rhat > 1.05)
        post_sum5 <-fit5$summary() %>% select(variable, mean, q5, q95)
        
      post5 <- as_draws_df(fit5$draws())
      post5_diff <- post5 %>% mutate(diff = `a_cog[1]` - `a_cog[2]`) %>%
        summarise(
          mean_diff = mean(diff), 
          sd_diff = sd(diff), 
          lower = quantile(diff, 0.025), 
          upper = quantile(diff, 0.975))
        
        loo5 <- loo(fit5$draws("log_lik"))

        loo_comp <- loo_compare(loo0, loo1, loo2, loo3, loo4, loo5)

        list(
          "time" = time_var,
          "y" = y,
          "post_sums" = list(post_sum0, post_sum1, post_sum2, post_sum3, post_sum4, post_sum5),
          "loo" = list(loo0, loo1, loo2, loo3, loo4, loo5),
          "post_diff" = list(post4_diff, post5_diff),
          "model_comp" = loo_comp,
          "posterior" = post5
        )
    })
  })
  save(shannon_models, file = here::here("rdata/16s_shannon_models.Rds"))
 } else {
  load(file = here::here("rdata/16s_shannon_models.Rds"))
}

# This file can be used for model comparison, to inspect summaries of the 
# different covariate structures and to work with the posterior distribution 
# of the model shown in the paper. This will be true also for the other linear 
# models in this script, e.g. here where we fit a similar model for the BRIEF:


# For BRIEF I use a multilevle structure to estimate the "true" BRIEF given our 
# 2 measurements:
if (!file.exists(here::here("rdata/16s_shannon_models_brief.Rds"))) {
    shannon_models_brief <- map(times, function(time_var) {
      
      y <- c("brief_total8_t", "brief_total10_t")
      df_lm <- meta %>% filter(time == time_var) %>%
        select(
          subject_id,
          all_of(y),
          edu,
          diversity_shannon,
          age10,
          bf
        ) %>% 
        mutate(
          subject_id = as.character(subject_id),
          edu = as.integer(as.factor(edu))
        ) %>%
        na.omit() %>%  
        mutate(subject_id = as.integer(as.factor(subject_id))) %>%
        arrange(subject_id)



      cog_lm <- df_lm %>% select(subject_id, all_of(y)) %>%
        pivot_longer(
          all_of(y), 
          names_to = "time", 
          values_to = "brief") %>%
          mutate(brief = standardize(brief))



      dat_list0 <- list(
        N_subj = length(df_lm$subject_id),
        N_Y_obs_total = dim(cog_lm)[1],
        cog_obs = cog_lm$brief,
        cog_obs_subj = cog_lm$subject_id
      )


      dat_list1 <- list(
        shannon = standardize(df_lm[["diversity_shannon"]]),
        N_subj = length(df_lm$subject_id),
        N_Y_obs_total = dim(cog_lm)[1],
        cog_obs = cog_lm$brief,
        cog_obs_subj = cog_lm$subject_id
      )

      dat_list2 <- list(
        N_subj = length(df_lm$subject_id),
        N_Y_obs_total = dim(cog_lm)[1],
        cog_obs = cog_lm$brief,
        cog_obs_subj = cog_lm$subject_id,
        shannon = standardize(df_lm[["diversity_shannon"]]),
        edu = df_lm[["edu"]],
        alpha = rep(2, 5)
      )
      
      dat_list3 <- list(
        N_subj = length(df_lm$subject_id),
        N_Y_obs_total = dim(cog_lm)[1],
        cog_obs = cog_lm$brief,
        cog_obs_subj = cog_lm$subject_id,
        edu = df_lm[["edu"]],
        alpha = rep(2, 5)
      )
      
      dat_list4 <- list(
        N_subj = length(df_lm$subject_id),
        N_Y_obs_total = dim(cog_lm)[1],
        cog_obs = cog_lm$brief,
        cog_obs_subj = cog_lm$subject_id,
        edu = df_lm[["edu"]],
        alpha = rep(2, 5),
        shannon = standardize(df_lm[["diversity_shannon"]]),
        bf = standardize(df_lm$bf)
      )


      file0 <- here::here("stan/shannon_brief_m0.stan")
      mod0 <- cmdstan_model(file0)
      fit0 <- mod0$sample(
        data = dat_list0,
        seed = 123,
        chains = 4,
        parallel_chains = 2,
        refresh = 500
      )
      fit0$summary() %>% filter(rhat > 1.05)
      post_sum0 <-fit0$summary() %>% select(variable, mean, q5, q95)
      loo0 <- loo(fit0$draws("log_lik"))



      file1 <- here::here("stan/shannon_brief_m1.stan")
      mod1 <- cmdstan_model(file1)
      fit1 <- mod1$sample(
        data = dat_list1,
        seed = 123,
        chains = 4,
        parallel_chains = 2,
        refresh = 500
      )

      fit1$summary() %>% filter(rhat > 1.05)
      post_sum1 <-fit1$summary() %>% select(variable, mean, q5, q95)
      loo1 <- loo(fit1$draws("log_lik"))


      file2 <- here::here("stan/shannon_brief_m2.stan")
      mod2 <- cmdstan_model(file2)

      fit2 <- mod2$sample(
        data = dat_list2,
        seed = 123,
        chains = 4,
        parallel_chains = 2,
        refresh = 500
      )

      fit2$summary() %>% filter(rhat > 1.05)
      post_sum2 <-fit2$summary() %>% select(variable, mean, q5, q95)
      loo2 <- loo(fit2$draws("log_lik"))



      file3 <- here::here("stan/shannon_brief_m3.stan")
      mod3 <- cmdstan_model(file3)

      fit3 <- mod3$sample(
        data = dat_list3,
        seed = 123,
        chains = 4,
        parallel_chains = 2,
        refresh = 500
      )

      fit3$summary() %>% filter(rhat > 1.05)
      post_sum3 <-fit3$summary() %>% select(variable, mean, q5, q95)
      loo3 <- loo(fit3$draws("log_lik"))




      file4 <- here::here("stan/shannon_brief_m4.stan")
      mod4 <- cmdstan_model(file4)

      fit4 <- mod4$sample(
        data = dat_list4,
        seed = 123,
        chains = 4,
        parallel_chains = 2,
        refresh = 500
      )

      fit4$summary() %>% filter(rhat > 1.05)
      post_sum4 <-fit4$summary() %>% select(variable, mean, q5, q95)
      loo4 <- loo(fit4$draws("log_lik"))
      
      post4 <- as_draws_df(fit4$draws())


      loo_comp <- loo_compare(loo0, loo1, loo2, loo3, loo4)

      list(
        "time" = time_var,
        "post_sums" = list(post_sum0, post_sum1, post_sum2, post_sum3, post_sum4),
        "loo" = list(loo0, loo1, loo2, loo3, loo4),
        "model_comp" = loo_comp,
        "posterior" = post4
      )
    })

  save(shannon_models_brief, file = here::here("rdata/16s_shannon_models_brief.Rds"))
 } else {
  load(file = here::here("rdata/16s_shannon_models_brief.Rds"))
}





# create plots for covariates and betas separately for the paper 

# extract covariates from DS models 
cov_post <- map_dfr(shannon_models, function(listobj) {
  map_dfr(listobj, function(listobj2) {

    listobj2$posterior %>% as_tibble() %>%
      mutate(male = `a_cog[2]` - `a_cog[1]`) %>%
      select(male, breastfeeding = bB, education = bE, age = bA) %>%
      mutate(time = listobj2$time, y = listobj2$y)
  })
})
cov_post_grouped <- cov_post %>% group_by(y) %>%
  nest() 

# extract betas from DS models 
b_post <- map_dfr(shannon_models, function(listobj) {
  map_dfr(listobj, function(listobj2) {
    listobj2$posterior %>% as_tibble() %>%
      select(shannon = bS) %>%
      mutate(time = listobj2$time, y = listobj2$y)
  })
})
b_post_grouped <- b_post %>% group_by(time, y) %>%
  nest()

# reshape 
b_post1 <- bind_cols(b_post_grouped$data[c(1:5)])
b_post2 <- bind_cols(b_post_grouped$data[c(6:10)])
b_post3 <- bind_cols(b_post_grouped$data[c(11:15)])
b_post <- map(list(b_post1, b_post2, b_post3), function(post) {
  colnames(post) <- glue("Day {unique(b_post_grouped$time)}")
  post
})

# extract covariates from BRIEF models 
cov_post2 <- map_dfr(shannon_models_brief, function(listobj) {
    listobj$posterior %>% as_tibble() %>%
      select(breastfeeding = bB, education = bE, shannon = bS) %>%
      mutate(time = listobj$time, y = listobj$y)
})
# reshape 
cov_post2_grouped <- cov_post2 %>% group_by(time) %>%
  nest()
cov_post2 <- bind_cols(cov_post2_grouped$data)
# now merge betas of DS and BRIEF together 
brief_shannon <- cov_post2 %>% select(contains("shannon"))
colnames(brief_shannon) <- colnames(b_post[[1]])
b_post[[4]] <- brief_shannon

# rename the 6 and 10 year columns
b_post <- map(b_post, function(post) {
  select(post, 
    "T1" = "Day 28", 
    "T2" = "Day 75", 
    "T3" = "Day 105",
    "T4" = "Day 2193",
    "T5" = "Day 3655"
  )
})


# plots betas
b_plots_infancy <- map2(b_post, list("DS Total LNS", "DS Total Backwards", "DS Total Forwards", "BRIEF"),  function(post, y) {
  post %>%
    mcmc_areas(
      pars = vars(contains("Day ")),
      prob = 0.95) +
      ggtitle(y) +
      geom_vline(aes(xintercept = 0), linetype = "dashed") +
      theme_bw(base_size = 20)
})

b_plots_childhood <- map2(b_post, list("DS Total LNS", "DS Total Backwards", "DS Total Forwards", "BRIEF"),  function(post, y) {
  post %>%
    mcmc_areas(
      pars = vars(contains("Year ")),
      prob = 0.95) +
      ggtitle(y) +
      geom_vline(aes(xintercept = 0), linetype = "dashed") +
      theme_bw(base_size = 20)
})


# average posterior of the covariates of BRIEF and DS together 
cov_brief <- cov_post2 %>% select(-contains("shannon")) %>%
  pivot_longer(contains("breastfeeding"), names_to = "time", values_to = "breastfeeding") %>%
  pivot_longer(contains("education"), names_to = "time2", values_to = "education") %>% select(breastfeeding, education)
# we use 0 as a placeholder to create the plot because the BRIEF was normed
cov_brief$male <- 0
cov_brief$age <- 0
cov_brief <- select(cov_brief, male, breastfeeding, education, age)
p_data_list <- as.list(cov_post_grouped$data)
p_data_list[[4]] <- cov_brief

# plots covariates 
cov_plots <- map2(
  as.list(c(as.character(cov_post_grouped$y), "BRIEF")), 
  p_data_list, function(y, data) {
  outcome <- ifelse(y == "total_lns", "DS Letter Number Sequencing", 
        ifelse(y == "total_bw", "DS Backwards", ifelse(y == "total_f", "DS Forwards", "BRIEF")))
  # may 2021: Carolina asked me to rename for publication:
  if (y == "BRIEF") {
    selector <- c("breastfeeding", "maternal education") 
    } else {
      selector <- c("child sex: male", "breastfeeding", "child age", "maternal education")
    } 
  # selector <- c("male", "breastfeeding", "education", "age")

  
  data %>% 
    select(
      "child sex: male" = "male", 
      "breastfeeding",
      "child age" = "age",
      "maternal education" = "education"
    ) %>%
    mcmc_areas(
      pars = vars(all_of(selector)),
      prob = 0.95
    ) +
    ggtitle(outcome) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    theme_bw(base_size = 20) 
})

# above plots are shown in the paper
save(b_plots_infancy, b_plots_childhood, cov_plots, file = here::here("rdata/shannon_plots.Rds"))






















#####################################################################
##########             3 PAM clustering                    ##########
#####################################################################

# I need to check for each timepoint separately if there is clustering
# so I first apply PAM for each otu table and then store plots to
# evaluate if there is clustering present. The function takes the meta
# data for the specific timepoint and a complete otu list as input
determine_k <- function(pseq, time_d, k = 10, method = pamkCBI, sil = "pam") {
    # use only samples for correponding timepoint 
    selector <- sample_data(pseq) %>% 
      as_tibble(rownames = NA) %>% 
      rownames_to_column("sample_id") %>%
      filter(time %in% time_d) %>%
      .$sample_id
    otu_clr <- prune_samples(selector, pseq) %>% otu_table() %>% t()
    
    # silhoutte index + calinski h index
    if (sil == "pam") {
      sil_i <- map(1:k, function(k) {
          model_pam <- pam(x = otu_clr, k = k)
          c(
            k = k, 
            sil = model_pam$silinfo$avg.width, 
            calinh = calinhara(otu_clr, model_pam$clustering))
          })      
    } else if (sil == "km") {
      sil_i <- map(1:k, function(k) {
        model_km <- kmeans(otu_clr, centers = k)
        si <- silhouette(model_km$cluster, dist(otu_clr))
        c(k = k,
          sil = ifelse(is.null(dim(si)), NA, mean(si[, "sil_width"])),
          calinh = calinhara(otu_clr, model_km$cluster))
      })
    }
    si_df <- sil_i %>% map_df(bind_rows)
    
    
    # here I calculate predictive strength for each time point
    ps_k <- prediction.strength(
        xdata = otu_clr, 
        Gmin = 2, 
        Gmax = k, 
        M = 50,
        clustermethod = method,
        classification = "centroid",
        cutoff = 0.8,
        distances = FALSE,
        count = FALSE)
    
    
    p <- tibble(k = 1:k, ps = ps_k$mean.pred) %>% 
        ggplot(aes(k, ps)) +
        scale_x_continuous(breaks = 1:k) +
        geom_hline(aes(yintercept = 0.9), linetype = "dashed", color = "#8da0cb", size = 1.5) +
        geom_hline(aes(yintercept = 0.5), linetype = "dashed", color = "#fc8d62", size = 1.5) +
        geom_point(color = "#8da0cb", size = 3) +
        geom_line(color = "#8da0cb", size = 1.5) +
        geom_point(data = si_df, aes(k, sil), color = "#fc8d62", size = 3) +
        geom_line(data = si_df, aes(k, sil), color = "#fc8d62", size = 1.5) +
        ggtitle(time_d) +
        ylab("") + xlab("") +
        theme_bw()
    
    p_cal <- ggplot(si_df, aes(k, calinh)) +
        scale_x_continuous(breaks = 1:k) +
        geom_point(color = "#66c2a5", size = 3) +
        geom_line(color = "#66c2a5", size = 1.5) +
        ggtitle(time_d) +
        ylab("") + xlab("") +
        theme_bw()    
        
    list(p, p_cal)
}




times <- c(28, 75, 105, 2193, 3655)
if(!file.exists(here::here("rdata/pam.Rds"))) {
  out_pam <- map(
    times, function(time) {
      determine_k(pseq = pseq_clr, time_d = time, k = 8, method = pamkCBI)
    }
  )
  save(out_pam, file = here::here("rdata/pam.Rds"))
 } else {
  load(here::here("rdata/pam.Rds"))
}
out_pam




if(!file.exists(here::here("rdata/kmeans.Rds"))) {
  out_km <- map(
    times, function(time) {
      determine_k(pseq = pseq_clr, time_d = time, k = 10, method = kmeansCBI, sil = "km")
    }
  )
  save(out_km, file = here::here("rdata/kmeans.Rds"))
 } else {
  load(here::here("rdata/kmeans.Rds"))
}


# explore whether we get clusters when we use all infant or childhood samples

times <- list(
  infancy = c(28, 75, 105),
  childhood = c(2193, 3655)
)

if(!file.exists(here::here("rdata/pam_all.Rds"))) {
  out_pam_all <- map(
    times, function(time) {
      determine_k(pseq = pseq_clr, time_d = time, k = 8, method = pamkCBI)
    }
  )
  save(out_pam_all, file = here::here("rdata/pam_all.Rds"))
 } else {
  load(here::here("rdata/pam_all.Rds"))
}

out_pam_all














#####################################################################
##########              4 DMM clustering                   ##########
#####################################################################


times <- c(28, 75, 105, 2193, 3655)
if (!file.exists(here::here("rdata/dmm.Rds"))) {
  p_lapl <- map(times, function(time_d) {
    
    selector <- sample_data(pseq) %>% 
      as_tibble(rownames = NA) %>% 
      rownames_to_column("sample_id") %>%
      filter(time == time_d) %>%
      .$sample_id
      
    asvs <- prune_samples(selector, pseq) %>% otu_table() %>% t()
    k <- c(1:10)
    fit_ds <- map(k, dmn, count = asvs)
    
    # laplace approximation
    fit_lap_ds <- map_dbl(fit_ds, laplace)
    # plot laplace approximation
    p <- tibble(k = k, laplace = fit_lap_ds) %>%
        ggplot(aes(k, laplace)) +
        geom_point() +
        geom_line() +
        theme_bw()
        
    list(
      "time" = time_d,
      "dmmfit" = fit_ds,
      "lpl_plot" = p
    )
  })

  save(p_lapl, file = here::here("rdata/dmm.Rds"))
 } else {
  load(here::here("rdata/dmm.Rds"))
}

# dmm indicated 2 clusters for each 6 and 10 years samples. Now we explore 
# using samples split by infancy/childhood 

times <- list(
    infacy = c(28, 75, 105),
    childhood = c(2193, 3655)
  )

if (!file.exists(here::here("rdata/dmm_split.Rds"))) {
  
  p_lapl <- map(times, function(time_d) {
    
    selector <- sample_data(pseq) %>% 
      as_tibble(rownames = NA) %>% 
      rownames_to_column("sample_id") %>%
      filter(time %in% time_d) %>%
      .$sample_id
      
    asvs <- prune_samples(selector, pseq) %>% otu_table() %>% t()
    k <- c(1:10)
    fit_ds <- map(k, dmn, count = asvs)
    
    # laplace approximation
    fit_lap_ds <- map_dbl(fit_ds, laplace)
    # plot laplace approximation
    p <- tibble(k = k, laplace = fit_lap_ds) %>%
        ggplot(aes(k, laplace)) +
        geom_point() +
        geom_line() +
        theme_bw()
        
    bestfit <- fit_ds[[which.min(fit_lap_ds)]]
    
    list(
      "time" = time_d,
      "dmmfit" = fit_ds,
      "lpl_plot" = p,
      "bestfit" = bestfit
    )
  })
  

  save(p_lapl, file = here::here("rdata/dmm_split.Rds"))
 } else {
  load(here::here("rdata/dmm_split.Rds"))
}


# now we obtain 3 in infancy and 4 in childhood
# visualize clusters using Aitchison distance for supplementary plots

# extract clusters and add to metadata 
df_k_split <- map(p_lapl, function(dmmlist) {
      df_sample <- sample_data(pseq_clr) %>% 
        as_tibble(rownames = NA) %>% 
        rownames_to_column("sample_id") %>%
        filter(time %in% dmmlist$time) 

      dmm_mod <- dmmlist$bestfit
      df_temp <- mixture(dmm_mod, assign = TRUE) %>% 
        as_tibble(rownames = NA) %>% 
        rownames_to_column("sample_id") %>%
        select(sample_id, dmm = value)

      left_join(df_sample, df_temp, by = "sample_id")
})

# to create the plots I need to add the dmm to pseq object
pseq_temp <- pseq_clr
df_k_all_split <- map_dfr(df_k_split, ~.x) %>% select(time, dmm, sample_id, subject_id)
sample_data(pseq_temp) <- sd_to_df(pseq_temp) %>% 
  select(sample_id) %>% 
  left_join(df_k_all_split, by = "sample_id") %>%
  df_to_sd()


p_dmmsplit_infants <- biplot(
  pseq_temp, 
  filter_samples = df_k_split$infacy$sample_id,
  color = "dmm",
  otu_alpha = 0,
  point_size = 5,
  colors = c("#1b9e77",
               "#d95f02",
               "#7570b3")
)

p_dmmsplit_children <-biplot(
  pseq_temp, 
  filter_samples = df_k_split$childhood$sample_id,
  color = "dmm",
  otu_alpha = 0,
  point_size = 5,
  colors = c("#a6cee3",
               "#1f78b4",
               "#b2df8a",
               "#33a02c")
)

p_dmmsplit_infants <- map(p_dmmsplit_infants, function(x) {
  x + theme_bw(base_size = 20)
})

p_dmmsplit_children <- map(p_dmmsplit_children, function(x) {
  x + theme_bw(base_size = 20)
})

save(p_dmmsplit_children, p_dmmsplit_infants, file = here::here("rdata/dmm_pca.Rds"))




















#####################################################################
##########  5 linear models ef ~ dmm + covariates          ##########
#####################################################################


# for this approach, we first need to identify the trajectories: If there are 
# enough infants/children that are constant in one of the 3/4 clusters across
# sample time points then we can perform this analyses. 

# first for infancy I filter out infants that remain stable in 1 cluster across 
dmminfant_id <- df_k_split[[1]] %>% 
  count(subject_id, dmm) %>% 
  count(subject_id) %>% 
  filter(n == 1) %>% .$subject_id

# looks good, we have equally distributed n and enough children who were only 
# in one cluster. however, we must assume that children with only 1 or 2 samples
# are also stable, which is false. 
df_k_split[[1]] %>% 
  count(subject_id, dmm) %>%
  filter(subject_id %in% dmminfant_id) %>% 
  count(dmm)
# just to double check
df_k_split[[1]] %>% count(subject_id, dmm) %>%
  ungroup() %>%
  pivot_wider(names_from = "dmm", values_from = "n") %>%
  mutate(
    dmm1 = ifelse(is.na(`1`), 0, `1`),
    dmm2 = ifelse(is.na(`2`), 0, `2`),
    dmm3 = ifelse(is.na(`3`), 0, `3`)) %>%
    select(subject_id, dmm1, dmm2, dmm3) %>%
    filter(subject_id %in% dmminfant_id)

dmm_infant <- df_k_split[[1]] %>% 
  filter(subject_id %in% dmminfant_id) %>% 
    select(subject_id, dmm_temp = dmm) %>% 
    arrange(subject_id) %>%
    group_by(subject_id) %>%
    summarise(dmm = max(dmm_temp))


# now the same for childhood:
dmmchild_id <- df_k_split[[2]] %>% 
  count(subject_id, dmm) %>% 
  count(subject_id) %>% 
  filter(n == 1) %>% .$subject_id

# looks good, we have equally distributed n and enough children who were only 
# in one cluster   
df_k_split[[2]] %>% 
  count(subject_id, dmm) %>%
  filter(subject_id %in% dmmchild_id) %>% 
  count(dmm)
# just to double check
df_k_split[[2]] %>% count(subject_id, dmm) %>%
  ungroup() %>%
  pivot_wider(names_from = "dmm", values_from = "n") %>%
  mutate(
    dmm1 = ifelse(is.na(`1`), 0, `1`),
    dmm2 = ifelse(is.na(`2`), 0, `2`),
    dmm3 = ifelse(is.na(`3`), 0, `3`),
    dmm4 = ifelse(is.na(`4`), 0, `4`)) %>%
    select(subject_id, dmm1, dmm2, dmm3, dmm4) %>%
    filter(subject_id %in% dmmchild_id)

dmm_child <- df_k_split[[2]] %>% 
  filter(subject_id %in% 
    dmmchild_id) %>% 
    select(subject_id, dmm_temp = dmm) %>% arrange(subject_id) %>%
    group_by(subject_id) %>%
    summarise(dmm = max(dmm_temp))




# now we can fit the models, starting with the DS outcome
if(!file.exists(here::here("rdata/16s_dmmsplit_models.Rds"))) {
  outcomes <- c("total_lns", "total_bw", "total_f")
  dfis <- c("infancy", "childhood")
  k_models_split <- map(outcomes, function(y) {
    map2(dfis, list(dmm_infant, dmm_child), function(dfi, df) {
      # create filter according to time:
      selector <- if (dfi == "infancy") dmminfant_id else dmmchild_id
      df_lm <- meta %>% filter(time == 28) %>% # just the vars obtained once
        select(
          subject_id,
          sample_id,
          all_of(y),
          edu,
          age10,
          sex,
          bf
        ) %>% 
        mutate(
          subject_id = as.character(subject_id),
          edu = as.integer(as.factor(edu))
        ) %>%
        filter(subject_id %in% selector) %>%
        left_join(df, by = "subject_id") %>%
        select(all_of(y), dmm, edu, age10, sex, bf) %>% 
        na.omit()



      dat_list0 <- list(
        cog = standardize(df_lm[[y]]),
        N = length(df_lm[[y]]),
        age = standardize(df_lm$age10)
      )
      
      dat_list1 <- list(
        cog = standardize(df_lm[[y]]),
        k = as.integer(df_lm[["dmm"]]),
        N = length(df_lm[[y]]),
        age = standardize(df_lm$age10),
        k_vec = ifelse(dfi == "infancy", 3, 4)
      )
      
      dat_list2 <- list(
        cog = standardize(df_lm[[y]]),
        k = df_lm[["dmm"]],
        edu = df_lm[["edu"]],
        alpha = rep(2, 5),
        N = length(df_lm[[y]]),
        age = standardize(df_lm$age10),
        k_vec = ifelse(dfi == "infancy", 3, 4)
      )
      
      dat_list3 <- list(
        cog = standardize(df_lm[[y]]),
        edu = df_lm[["edu"]],
        alpha = rep(2, 5),
        N = length(df_lm[[y]]),
        age = standardize(df_lm$age10)
      )
      
      dat_list4 <- list(
        cog = standardize(df_lm[[y]]),
        edu = df_lm[["edu"]],
        k = df_lm[["dmm"]],
        alpha = rep(2, 5),
        N = length(df_lm[[y]]),
        age = standardize(df_lm$age10),
        sex = as.integer(as.factor(df_lm$sex)),
        k_vec = ifelse(dfi == "infancy", 3, 4)
      )
      
      dat_list5 <- list(
        cog = standardize(df_lm[[y]]),
        edu = df_lm[["edu"]],
        k = df_lm[["dmm"]],
        alpha = rep(2, 5),
        N = length(df_lm[[y]]),
        age = standardize(df_lm$age10),
        sex = as.integer(as.factor(df_lm$sex)),
        k_vec = ifelse(dfi == "infancy", 3, 4),
        bf = df_lm$bf
      )
      
      file0 <- here::here("stan/k_m0.stan")
      mod0 <- cmdstan_model(file0)
      fit0 <- mod0$sample(
        data = dat_list0,
        seed = 123,
        chains = 4,
        parallel_chains = 2,
        refresh = 500
      )
      fit0$summary() %>% filter(rhat > 1.05)
      post_sum0 <-fit0$summary() %>% select(variable, mean, q5, q95)
      loo0 <- loo(fit0$draws("log_lik"))
      
      
      
      file1 <- here::here("stan/k_m1.stan")
      mod1 <- cmdstan_model(file1)
      fit1 <- mod1$sample(
        data = dat_list1,
        seed = 123,
        chains = 4,
        parallel_chains = 2,
        refresh = 500
      )
      
      fit1$summary() %>% filter(rhat > 1.05)
      post_sum1 <-fit1$summary() %>% select(variable, mean, q5, q95)
      loo1 <- loo(fit1$draws("log_lik"))
      
      post1 <- as_draws_df(fit1$draws())
      post1_diff <- post1 %>% mutate(diff = `a[1]` - `a[2]`) %>%
        summarise(
          mean_diff = mean(diff), 
          sd_diff = sd(diff), 
          lower = quantile(diff, 0.025), 
          upper = quantile(diff, 0.975))
      
      file2 <- here::here("stan/k_m2.stan")
      mod2 <- cmdstan_model(file2)
      fit2 <- mod2$sample(
        data = dat_list2,
        seed = 123,
        chains = 4,
        parallel_chains = 2,
        refresh = 500
      )
      
      fit2$summary() %>% filter(rhat > 1.05)
      post_sum2 <-fit2$summary() %>% select(variable, mean, q5, q95)
      loo2 <- loo(fit2$draws("log_lik"))
      
      post2 <- as_draws_df(fit2$draws())
      post2_diff <- post2 %>% mutate(diff = `a[1]` - `a[2]`) %>%
        summarise(
          mean_diff = mean(diff), 
          sd_diff = sd(diff), 
          lower = quantile(diff, 0.025), 
          upper = quantile(diff, 0.975))
      
      file3 <- here::here("stan/k_m3.stan")
      mod3 <- cmdstan_model(file3)
      fit3 <- mod3$sample(
        data = dat_list3,
        seed = 123,
        chains = 4,
        parallel_chains = 2,
        refresh = 500
      )
      
      fit3$summary() %>% filter(rhat > 1.05)
      post_sum3 <-fit3$summary() %>% select(variable, mean, q5, q95)
      loo3 <- loo(fit3$draws("log_lik"))
      
      
      
      file4 <- here::here("stan/k_m4.stan")
      mod4 <- cmdstan_model(file4)
      fit4 <- mod4$sample(
        data = dat_list4,
        seed = 123,
        chains = 4,
        parallel_chains = 2,
        refresh = 500
      )
      
      fit4$summary() %>% filter(rhat > 1.05)
      post_sum4 <-fit4$summary() %>% select(variable, mean, q5, q95)
      loo4 <- loo(fit4$draws("log_lik"))
      post4 <- as_draws_df(fit4$draws())
      
      file5 <- here::here("stan/k_m5.stan")
      mod5 <- cmdstan_model(file5)
      fit5 <- mod5$sample(
        data = dat_list5,
        seed = 123,
        chains = 4,
        parallel_chains = 2,
        refresh = 500
      )
      
      fit5$summary() %>% filter(rhat > 1.05)
      post_sum5 <-fit5$summary() %>% select(variable, mean, q5, q95)
      loo5 <- loo(fit5$draws("log_lik"))
      post5 <- as_draws_df(fit5$draws())
      
      
      loo_comp <- loo_compare(loo0, loo1, loo2, loo3, loo4, loo5)
      
      list(
        "time" = ifelse(dfi == "infancy", "infancy", "childhood"),
        "y" = y,
        "post_sums" = list(post_sum0, post_sum1, post_sum2, post_sum3, post_sum4, post_sum5),
        "loo" = list(loo0, loo1, loo2, loo3, loo4, loo5),
        "model_comp" = loo_comp,
        #"post_diff" = list(post1_diff, post2_diff, post4_diff1, post4_diff2),
        "post" = list(post1, post2, post4, post5)
      )
    })
  })
  save(k_models_split, file = here::here("rdata/16s_dmmsplit_models.Rds"))
 } else {
  load(file = here::here("rdata/16s_dmmsplit_models.Rds"))
}



dmm_k_plots <- map(k_models_split, function(models) {
  map(models, function(model) {
    time_var <- model$time
    y <- ifelse(model$y == "total_lns", "DS Letter Number Sequencing", 
          ifelse(model$y == "total_bw", "DS Backwards", "DS Forwards"))
    
    # please note that the numbers of the clusters are changed because a 
    # colleague of mine first published this DMM approach with the same data 
    # while switching the cluster names (but the clusters are identical).
    # to not confuse readers that read both papers, I stick to her cluster 
    # names.
    if (time_var == "infancy") {
      df <- model$post[[3]] 
      df <- df %>%
        mutate(
          "k1-k2" = `a[1]` - `a[3]`,
          "k1-k3" = `a[1]` - `a[3]`,
          "k2-k3" = `a[2]` - `a[3]`,
        )
    } else {
      df <- model$post[[3]] 
      df <- df %>%
        mutate(
          "k1-k2" = `a[1]` - `a[4]`,
          "k1-k3" = `a[1]` - `a[2]`,
          "k1-k4" = `a[1]` - `a[3]`,
          
          "k2-k3" = `a[4]` - `a[2]`,
          "k2-k4" = `a[4]` - `a[3]`,
          
          "k3-k4" = `a[2]` - `a[3]`,  
        )
    }    
      
      p <- mcmc_areas(df,
        pars = if (time_var == "infancy") c("k1-k2", "k1-k3", "k2-k3")
               else c("k1-k2", "k1-k3","k1-k4", "k2-k3","k2-k4", "k3-k4"),
        prob = 0.95) + ggtitle(y) +
        geom_vline(aes(xintercept = 0), linetype = "dashed") +
        theme_bw(base_size = 20) 
        
      if (model$y != "total_f") {
        p <- p + theme(axis.title.y=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank())
      }
      
      p

  })
})

dmm_k_plots
save(dmm_k_plots, file = here::here("rdata/dmm_k_plots.Rds"))


# same for BRIEF where I can implement again the multilevel structure
if (!file.exists(here::here("rdata/16s_ksplit_models_brief.Rds"))) {
  dfis <- c("infancy", "childhood")
  k_models_brief_split <- map2(dfis, list(dmm_infant, dmm_child), function(dfi, df) {
    
    y <- c("brief_total8_t", "brief_total10_t")
    df_lm <- meta %>% filter(time == 28) %>%
      select(
        subject_id,
        all_of(y),
        edu,
        bf
      ) %>% 
      mutate(
        subject_id = as.character(subject_id),
        edu = as.integer(as.factor(edu))
      ) %>%
      left_join(df, by = "subject_id") %>%
      mutate(dmm = as.integer(dmm)) %>%
      na.omit() %>%
      mutate(subject_id = as.integer(as.factor(subject_id)))


    cog_lm <- df_lm %>% select(subject_id, brief_total8_t, brief_total10_t) %>%
      pivot_longer(
        c(brief_total8_t, brief_total10_t), 
        names_to = "time", 
        values_to = "brief") %>%
        mutate(brief = standardize(brief)) 


    df_temp <- select(df_lm, dmm, edu, bf) 
    


    dat_list0 <- list(
      N_subj = length(df_temp$dmm),
      N_Y_obs_total = dim(cog_lm)[1],
      cog_obs = cog_lm$brief,
      cog_obs_subj = cog_lm$subject_id
    )
    
    
    dat_list1 <- list(
      k = df_temp[["dmm"]],
      N_subj = length(df_temp$dmm),
      N_Y_obs_total = dim(cog_lm)[1],
      cog_obs = cog_lm$brief,
      cog_obs_subj = cog_lm$subject_id,
      k_vec = ifelse(dfi == "infancy", 3, 4)
    )
    
    dat_list2 <- list(
      N_subj = length(df_temp$dmm),
      N_Y_obs_total = dim(cog_lm)[1],
      cog_obs = cog_lm$brief,
      cog_obs_subj = cog_lm$subject_id,
      k = df_temp[["dmm"]],
      edu = df_temp[["edu"]],
      alpha = rep(2, 5),
      k_vec = ifelse(dfi == "infancy", 3, 4)
    )
    
    dat_list3 <- list(
      N_subj = length(df_temp$dmm),
      N_Y_obs_total = dim(cog_lm)[1],
      cog_obs = cog_lm$brief,
      cog_obs_subj = cog_lm$subject_id,
      edu = df_temp[["edu"]],
      alpha = rep(2, 5)
    )
    
    dat_list4 <- list(
      N_subj = length(df_temp$dmm),
      N_Y_obs_total = dim(cog_lm)[1],
      cog_obs = cog_lm$brief,
      cog_obs_subj = cog_lm$subject_id,
      edu = df_temp[["edu"]],
      alpha = rep(2, 5),
      bf = standardize(df_temp[["bf"]])
    )
    
    
    file0 <- here::here("stan/k_brief_m0.stan")
    mod0 <- cmdstan_model(file0)
    fit0 <- mod0$sample(
      data = dat_list0,
      seed = 123,
      chains = 4,
      parallel_chains = 2,
      refresh = 500
    )
    fit0$summary() %>% filter(rhat > 1.05)
    post_sum0 <-fit0$summary() %>% select(variable, mean, q5, q95)
    loo0 <- loo(fit0$draws("log_lik"))
    
    
    
    file1 <- here::here("stan/k_brief_m1.stan")
    mod1 <- cmdstan_model(file1)
    fit1 <- mod1$sample(
      data = dat_list1,
      seed = 123,
      chains = 4,
      parallel_chains = 2,
      refresh = 500
    )
    
    fit1$summary() %>% filter(rhat > 1.05)
    post_sum1 <-fit1$summary() %>% select(variable, mean, q5, q95)
    loo1 <- loo(fit1$draws("log_lik"))
    
    post1 <- as_draws_df(fit1$draws())

    
    file2 <- here::here("stan/k_brief_m2.stan")
    mod2 <- cmdstan_model(file2)
    
    fit2 <- mod2$sample(
      data = dat_list2,
      seed = 123,
      chains = 4,
      parallel_chains = 2,
      refresh = 500
    )
    
    fit2$summary() %>% filter(rhat > 1.05)
    post_sum2 <-fit2$summary() %>% select(variable, mean, q5, q95)
    loo2 <- loo(fit2$draws("log_lik"))
    
    
    post2 <- as_draws_df(fit2$draws())

    
    
    file3 <- here::here("stan/k_brief_m3.stan")
    mod3 <- cmdstan_model(file3)
    
    fit3 <- mod3$sample(
      data = dat_list3,
      seed = 123,
      chains = 4,
      parallel_chains = 2,
      refresh = 500
    )
    
    fit3$summary() %>% filter(rhat > 1.05)
    post_sum3 <-fit3$summary() %>% select(variable, mean, q5, q95)
    loo3 <- loo(fit3$draws("log_lik"))
    
    file4 <- here::here("stan/k_brief_m4.stan")
    mod4 <- cmdstan_model(file4)
    
    fit4 <- mod4$sample(
      data = dat_list4,
      seed = 123,
      chains = 4,
      parallel_chains = 2,
      refresh = 500
    )
    
    fit4$summary() %>% filter(rhat > 1.05)
    post_sum4 <-fit4$summary() %>% select(variable, mean, q5, q95)
    loo4 <- loo(fit4$draws("log_lik"))
    
    
    loo_comp <- loo_compare(loo0, loo1, loo2, loo3, loo4)
    
    list(
      "time" = ifelse(dfi == "infancy", "infancy", "childhood"),
      "post_sums" = list(post_sum0, post_sum1, post_sum2, post_sum3, post_sum4),
      "loo" = list(loo0, loo1, loo2, loo3, loo4),
      "model_comp" = loo_comp,
      "post" = list(post1, post2)
    )
    

  })
  save(k_models_brief_split, file = here::here("rdata/16s_ksplit_models_brief.Rds"))
 } else {
  load(file = here::here("rdata/16s_ksplit_models_brief.Rds"))
}





dmm_brief_plots <- map(k_models_brief_split, function(model) {
    time_var <- model$time
    y <- "BRIEF"
    
    if (time_var == "infancy") {
      df <- model$post[[2]]
      df <- df %>%
        mutate(
          "k1-k2" = `a_cog[1]` - `a_cog[2]`,
          "k1-k3" = `a_cog[1]` - `a_cog[3]`,          
          "k2-k3" = `a_cog[2]` - `a_cog[3]`
        )
    } else {
      df <- model$post[[2]] 
      df <- df %>%
        mutate(
          "k1-k2" = `a_cog[1]` - `a_cog[4]`, # also here naming as above
          "k1-k3" = `a_cog[1]` - `a_cog[2]`,
          "k1-k4" = `a_cog[1]` - `a_cog[3]`,
          
          "k2-k3" = `a_cog[4]` - `a_cog[2]`,
          "k2-k4" = `a_cog[4]` - `a_cog[3]`,
          
          "k3-k4" = `a_cog[2]` - `a_cog[3]`,
        )
    }    


    mcmc_areas(df,
      pars = if (time_var == "infancy") c("k1-k2", "k1-k3", "k2-k3") 
             else c("k1-k2", "k1-k3","k1-k4", "k2-k3","k2-k4", "k3-k4"),
      prob = 0.95) +
      geom_vline(aes(xintercept = 0), linetype = "dashed") +
      theme_bw(base_size = 20) +
      theme(axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank())
      

})


dmm_brief_plots

save(dmm_brief_plots, dmm_k_plots, file = here::here("rdata/dmm_k_plots.Rds"))






















#####################################################################
##########          6 Shannon plots and correlation        ##########
#####################################################################

# correlation 
meta %>% select(subject_id, time, diversity_shannon) %>%
  pivot_wider(names_from = time, values_from = diversity_shannon) %>%
  select(-subject_id) %>%
  cor(use = "pairwise.complete.obs") %>%
  as.data.frame() %>%
  mutate(across(where(is.numeric), round, 2))


# create path plot to visualize AD over time 
# get exact age when stool samples were created:
age_path <- here::here("data/age_per_sample_infancy.xlsx")
inf_age <- readxl::read_excel(age_path) %>%
  select(
    subject_id = ID,  
    day28 = "28 days",
    day75 = "CC -2 days",
    day105  = "CC+28 days"
  ) %>%
  mutate(subject_id = as.character(subject_id)) %>%
  filter(!is.na(subject_id)) %>%
  pivot_longer(contains("day"), names_to = "time", values_to = "age") %>%
  mutate(time = as.numeric(str_replace(time, "day", ""))) %>%
  mutate(age = ifelse(is.na(age), time, age))


# create df needed for that plot
adpdf <- meta %>% 
  select(subject_id, Shannon = diversity_shannon, time) %>%   left_join(inf_age, by = c("subject_id", "time")) %>%
  filter(time %in% c(28, 75, 105)) %>%
  mutate(time = as.factor(time))
# to avoid cluttering we need to split this plot into multiple plots (ids)
df_temp <- tibble(
  subject_id = adpdf$subject_id %>% unique(),
  id = c(rep(seq(1:20), each = 9), rep(21, 5))
)
adpdf <- full_join(adpdf, df_temp, by = "subject_id")
adpdf$subject_id <- as.factor(as.integer(as.factor(adpdf$subject_id)))
adpdf$lower <- quantile(adpdf$Shannon, 0.25)
adpdf$median <- quantile(adpdf$Shannon, 0.5)
adpdf$upper <- quantile(adpdf$Shannon, 0.75)

shannon_stab1 <- map(1:21, function(idd) {
  adpdf %>% 
    filter(id == idd) %>%
    arrange(id, subject_id, time) %>%
    ggplot(aes(age, Shannon, group = subject_id, color = subject_id)) +
      geom_violin(data = adpdf, aes(y = Shannon), group = 0, color = "black") +
      geom_hline(aes(yintercept = median), linetype = "dashed") +
      geom_hline(aes(yintercept = lower), linetype = "dashed") +
      geom_hline(aes(yintercept = upper), linetype = "dashed") +
      geom_point(size = 8) +
      geom_path(size = 0.5) +
      scale_x_continuous(limits = c(10, 170)) +
      scale_y_continuous(limits = c(0.125, 4)) +
      scale_color_manual(
        values = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6')) +
      theme_bw(base_size = 20) +
      theme(legend.position = "none")
})


# same for childhood 
age6 <- readxl::read_excel(here::here("data/age_6.xlsx")) %>%
  select(subject_id = ID, age6 = AGE_M)
age10 <- readxl::read_excel(here::here("data/age_10.xlsx")) %>%
  select(subject_id = ID, age10 = AGE_M)
agedf <- full_join(age6, age10, by = "subject_id") %>%
  pivot_longer(-subject_id, names_to = "time", values_to = "age") %>%
  mutate(time = ifelse(time == "age6", 2193, 3655)) %>%
  mutate(subject_id = as.character(subject_id))

adpdf <- meta %>% 
  select(subject_id, Shannon = diversity_shannon, time, sex) %>%  
  filter(time %in% c(3655, 2193)) %>%
  left_join(agedf, by = c("subject_id", "time")) %>%
  filter(time %in% c(2193, 3655)) %>%
  mutate(time = as.factor(time)) %>%
  mutate(age = ifelse(age >= 1000, NA, age)) # missing values are 99999

df_temp <- tibble(
  subject_id = adpdf$subject_id %>% unique(),
  id = c(rep(seq(1:20), each = 8), rep(21, 5))
)
adpdf <- full_join(adpdf, df_temp, by = "subject_id")
adpdf$subject_id <- as.factor(as.integer(as.factor(adpdf$subject_id)))
adpdf$lower <- quantile(adpdf$Shannon, 0.25)
adpdf$median <- quantile(adpdf$Shannon, 0.5)
adpdf$upper <- quantile(adpdf$Shannon, 0.75)
shannon_stab2 <- map(1:21, function(idd) {
  adpdf %>% 
    filter(id == idd) %>%
    arrange(id, subject_id, time) %>%
    ggplot(aes(age, Shannon, group = subject_id, color = subject_id)) +
      geom_violin(data = adpdf, aes(y = Shannon), group = 0, color = "black") +
      geom_hline(aes(yintercept = median), linetype = "dashed") +
      geom_hline(aes(yintercept = lower), linetype = "dashed") +
      geom_hline(aes(yintercept = upper), linetype = "dashed") +
      geom_point(size = 8) +
      geom_path(size = 0.5) +
      #geom_text(size = 5) +

      scale_x_continuous(limits = c(50, 140)) +
      scale_y_continuous(limits = c(2, 4.5)) +
      scale_color_manual(
        values = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6')) +
      theme_bw(base_size = 20) +
      theme(legend.position = "none")
})

save(shannon_stab1, shannon_stab2, file = here::here("rdata/shannon_stab.Rds"))







#####################################################################
##########                  volatility plots               ##########
#####################################################################
# get all IDs we have to plot
bpids <- meta %>% filter(
    !is.na(time),
    (!is.na(total_f) | 
      !is.na(total_bw) | 
        !is.na(total_lns) |
          !is.na(brief_total8_t) |
            !is.na(brief_total10_t))
  ) %>% distinct(subject_id) %>% .$subject_id


# create a pseq for the biplot function
sd_temp <- meta %>% filter(subject_id %in% bpids) %>%
 mutate(T = 
   ifelse(time == 28, "1", 
    ifelse(time == 75, "2", 
      ifelse(time == 105, "3", 
        ifelse(time == 2193, "4", 
          ifelse(time == 3655, "5",NA))))),
          T = factor(T, levels = c("1", "2", "3", "4", "5"))) %>% 
  right_join(sd_to_df(pseq_clr), by = c("subject_id", "time")) %>%
  filter(!is.na(T)) %>%
  arrange(subject_id, T)

# to split subjects between plots to avoid overplotting
iid <- sd_temp %>% select(subject_id) %>% distinct() %>%
  mutate(
    iid = as.integer(as.factor(subject_id)),
    seq = rep(1:34, each = 5)
  )

sd_temp <- left_join(sd_temp, iid, by = c("subject_id"))
pseq_temp <- prune_samples(sd_temp$sample_id, pseq_clr)
sample_data(pseq_temp) <- df_to_sd(sd_temp)

# set axis limits equal for all plots 
yplus <- 6
yminus <- -10
xminus <- -10
xplus <- 13
beta_vol_plots <- map(1:34, function(seqiid) {
  biplot(
    pseq_temp, 
    filter_samples = filter(sd_temp, seq == seqiid) %>% .$sample_id, 
    connect_series = "time",
    shape = "subject_id",
    otu_alpha = 0,
    path_size = 1,
    label = "T",
    point_size = 15,
    alpha = 0.65,
    # text = TRUE,
    # text_size = 10,
    arrow_size = 1,
    #color = "T",
    #colors = c('#fef0d9','#fdcc8a','#fc8d59','#e34a33','#b30000'),
    path_alpha = 0.5
  )[[1]] +  theme_bw(base_size = 30) + 
            theme(legend.position="none") + 
              xlim(c(xminus, xplus)) + ylim(c(yminus, yplus))

})




save(beta_vol_plots, file = here::here("rdata/beta_vol_plots.Rds"))
 