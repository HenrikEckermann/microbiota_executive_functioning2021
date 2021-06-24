#####################################################################
##########                       Content                   ##########
#####################################################################

#1  import data that has been analyzed using NGTAX2 and stored as pseq object
#2 check how much data we loose when using complete case analyses
#3 save object to use for analyses in the other scripts


library(phyloseq)
library(microbiome)
library(tidyverse)
library(lubridate)
source(here::here("R/mb_helper.R"))

# load 16S data 
if (!file.exists(here::here("rdata/bibo_16s.Rds"))) {
  map_file <- here::here("data/bibo/ngtax1_out_yang/2019_0054 to 2019_0062 20190074 75 82 mapping file.csv")
  biom_file <- here::here("data/bibo/ngtax1_out_yang/Galaxy224-[NG-Tax__BIBO_1m-10y_(20200116)].biom1")
  tree_file <- read_tree(here::here("data/bibo/ngtax1_out_yang/Galaxy226-[Tree_files__NG-Tax__BIBO_1m-10y_(20200116)]/all_otus.tree"))

  bibo <- read_phyloseq(otu.file = biom_file,
                              taxonomy.file = NULL,
                              metadata.file = map_file,
                              type = "biom")
  bibo <- merge_phyloseq(bibo, tree_file)                           

  # the alphabetic letter in each SampleID means the timepoint when stool sample
  # was collected: a, 1m; b, 3m; c, 4m; d, 6y; e, 10y.

  sample_data(bibo) <- bibo %>% sd_to_df() %>%
    mutate(
      subject_id = ifelse(
        grepl("^\\w\\d\\d\\d$", sample_id), 
        str_sub(sample_id, 2, 4), 
        NA),
      time = ifelse(
        grepl("^a\\d\\d\\d$", sample_id), 28, ifelse(
          grepl("^b\\d\\d\\d$", sample_id), 75, ifelse(
            grepl("^c\\d\\d\\d$", sample_id), 105, ifelse(
              grepl("^d\\d\\d\\d$", sample_id), 2193, ifelse(
                grepl("^e\\d\\d\\d$", sample_id), 3655, NA)))))) %>%
    df_to_sd()

  ad <- alpha(bibo)
  sample_data(bibo) <- rownames_to_column(ad, "sample_id") %>%
    left_join(sd_to_df(bibo),
    by = "sample_id") %>%
    df_to_sd()
  save(bibo, file = here::here("rdata/bibo_16s.Rds"))
 } else {
  load(file = here::here("rdata/bibo_16s.Rds"))
}


#####################################################################
##########               combine with ef and covariates    ##########
#####################################################################

load(here::here("data/breastfeeding/anybf.Rds"))
load(here::here("data/medu.Rds"))
ef_data <- readRDS(here::here("data/meta_data/ef_data.RDS")) %>%
    mutate(subject_id = as.character(subject_id)) %>%
    filter(subject_id %in% sd_to_df(bibo)$subject_id) 
df_bf <- rename(df_bf, subject_id = id)
# age ef
age10 <- readxl::read_excel(
  here::here("data/brief/brief_10_years/age_brief10.xlsx"),
  col_types = c("text", "date", "text")) %>% mutate(
    birthdate = dmy(birthdate),
    brief10_age = ymd(brief10_age),
    age10 = brief10_age - birthdate)
age_cog <- readxl::read_excel(
  here::here("data/brief/brief_8_years/age_brief8.xlsx"),
  col_types = c("text", "text", "date")) %>% 
  mutate(frm_date = date(ymd_hms(frm_date))) %>%
  select(subject_id = bibo0002, frm_date) %>%
  full_join(age10, by = "subject_id") %>%
  mutate(age8 = frm_date - birthdate) %>%
  mutate(age8 = as.numeric(age8), age10 = as.numeric(age10)) %>%
  select(subject_id, age8, age10)

meta <- full_join(sd_to_df(bibo), ef_data, by = "subject_id") %>% 
  full_join(maternaledu, by = "subject_id") %>%
  full_join(df_bf, by = "subject_id") %>%
  select(-contains("age")) %>%
  full_join(age_cog, by = "subject_id")

#####################################################################
##########               check missingness                 ##########
#####################################################################


# check how much data we would loose by incl age 
outcomes <- c("total_lns", "total_f", "total_bw")
times <- c(28, 75, 105, 2193, 3655)
map_dfr(outcomes, function(y) {
  map(times, function(time_var) {
      withall <- meta %>% 
        select(
          all_of(y), 
          age10, 
          subject_id, 
          time, 
          edu, 
          diversity_shannon, 
          bf, 
          sex) %>%
        filter(time == time_var) %>%
        na.omit() %>% dim()
      
      without_age <- meta %>% 
        select(
          all_of(y), 
          "subject_id", 
          "time", 
          edu, 
          diversity_shannon, 
          bf, 
          sex) %>%
        filter(time == time_var) %>%
        na.omit() %>% dim()
        
      without_edu <- meta %>% 
        select(
          all_of(y), 
          age10, 
          subject_id, 
          time, 
          diversity_shannon, 
          bf, 
          sex) %>%
        filter(time == time_var) %>%
        na.omit() %>% dim()
      
      without_bf <- meta %>% 
        select(
          all_of(y), 
          age10, 
          subject_id, 
          time, 
          diversity_shannon, 
          edu, 
          sex) %>%
        filter(time == time_var) %>%
        na.omit() %>% dim()
        
      max_possible <- meta %>% 
        select(all_of(y), subject_id, time, diversity_shannon) %>%
        filter(time == time_var) %>%
        na.omit() %>% dim()
        
    list(
      max_possible = max_possible[1],
      withall = withall[1], 
      without_age = without_age[1], 
      without_edu = without_edu[1],
      without_bf = without_bf[1],
      diff = max_possible[1] - withall[1])
  })
})

# adding the covariates causes barely missingness when using cc analyses 
# max is 2.4% 
# store meta data object for further analyses
meta <- meta %>% arrange(subject_id) %>%
  select(
    subject_id, 
    sample_id,
    time, 
    edu,
    sex,
    total_lns, 
    total_f, 
    total_bw, 
    brief_total8_t,
    brief_total10_t,
    diversity_shannon,
    age10,
    bf
  ) %>% filter(!is.na(subject_id), !is.na(time)) # get rid of mock samples 
save(meta, file = here::here("rdata/16s_meta.Rds"))

# for another script, I need the taxtable 
taxtable <- tax_table(bibo) %>% as.data.frame() %>% rownames_to_column("taxid")
save(taxtable, file = here::here("rdata/taxtable.Rds"))
