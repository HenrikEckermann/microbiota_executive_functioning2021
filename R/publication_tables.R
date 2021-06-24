library(tidyverse)
library(huxtable)
library(patchwork)
library(gtsummary)
library(glue)

# 16 S data
load(here::here("rdata/bibo_16s.Rds"))
load(here::here("rdata/16s_meta.Rds"))



#####################################################################
##########               Sample size                       ##########
#####################################################################



Time <- c("28 days", "75 days", "105 days", "6 years", "10 years")
outcomes <- c("DS F", "DS BW", "DS LNS", "BRIEF")
BRIEF <- c(144, 131, 129, 144, 146) %>% as.integer()
DS_F <- c(135, 125, 123, 139, 146) %>% as.integer()
DS_BW <- c(132, 123, 121, 137, 146) %>% as.integer()
DS_LNS <- c(131, 122, 120, 136, 145) %>% as.integer()
sample_size <- tibble(
    Time = Time,
    "DS F" = DS_F,
    "DS BW" = DS_BW,
    "DS LNS" = DS_LNS,
    BRIEF = BRIEF
)

#####################################################################
##########               RF table                          ##########
#####################################################################

load(here::here("rdata/rf_summary_df.Rds"))
summary_df %>% 
  filter(statistic == "pearson1") %>%
  mutate(y = factor(y, levels = c("total_f", "total_bw", "total_lns", "brief"))) %>% 
  arrange(time, y) %>%
  mutate(
    y = 
      ifelse(y == "brief", "BRIEF", 
      ifelse(y == "total_f", "DS Forwards", 
      ifelse(y == "total_bw", "DS Backwards", 
      ifelse(y == "total_lns", "DS LNS", NA)))),
    time = 
      ifelse(time == 28, "28 days",
      ifelse(time == 75, "75 days", 
      ifelse(time == 105, "105 days",
      ifelse(time == 2193, "6 years", "10 years"))))
    ) %>%
  select(Outcome = y, "Age Microbiota Sample" = time, "Median" = median, "Mean" = mean, "SD" = sd, p = pvalue, q) %>%
  mutate(across(where(is.numeric), function(x) format(round(x , 2), nsmall = 2))) %>%
  as_hux() %>%
  set_all_padding(4) %>%
  set_bold(1, everywhere) %>%
  set_bottom_border(1, everywhere) %>%
  add_footnote("Note. SD = standard deviation. DS = Digit Span.") %>%
  set_caption("Correlation between Random Forest prediction and real data.")




#####################################################################
##########              demographic table                  ##########
#####################################################################
times <- c(28, 75, 105, 2193, 3655)
outcomes <- c(
  "total_lns", 
  "total_bw", 
  "total_f", 
  "brief_total8_t", 
  "brief_total10_t"
)
# get ids of involved subjects 
ids <- map_dfr(times, function(time_var) {
  map_dfr(outcomes, function(y) {
    df_temp <- filter(meta, time == time_var) %>%
      select(subject_id, all_of(y), edu, sex, diversity_shannon, age10) %>%
      na.omit() 
  })}) %>%
  .$subject_id %>% unique()



# antibiotics
abx1 <- read_csv(here::here("data/abx/abx_past_year.csv")) %>%
  select(
    sample_id = age.and.ID, 
    time = Age, 
    subject_id = ID, 
    antibiotics = Whether.took.antibiotics.in.the.past.one.year
  )
abx2 <- read_csv(here::here("data/abx/abx_birth_to_sample.csv")) %>%
  select(
    sample_id = age.and.ID, 
    time = Age, 
    subject_id = ID, 
    antibiotics = Whether.took.antibiotics.from.birth.to.stool.sample.collection
  )
bind_rows(abx1, abx2) %>%
  distinct(antibiotics)
abx <- bind_rows(abx1, abx2) %>%
  filter(subject_id %in% ids) %>%
  mutate(time = ifelse(time == "1m", "28 days", 
            ifelse(time == "3m", "75 days", 
              ifelse(time == "4m", "105 days", 
                ifelse(time == "6y", "6 years",
                  ifelse(time == "10y", "10 years", NA))))),
          antibiotics = ifelse(antibiotics %in% c("yes", "Yes"), TRUE, ifelse(antibiotics %in% c("no", "No"), FALSE, NA))) %>%
  select(subject_id, time, antibiotics) %>%
  pivot_wider(names_from = time, values_from = antibiotics) %>%
  mutate(subject_id = as.character(subject_id))

# gather all data for the table 
part1 <- foreign::read.spss(
  here::here("data/meta_data/bibo_confounders.sav"), 
  to.data.frame = T) %>%
  select(
    subject_id = ID,
    Birthweight = BIRTHWEIGTH,
    "Smoking during pregnancy" = preSMOKING,
    "Alcohol during pregnancy" = preALCOHOL,
    # Apgar5min = Apgar5min,
    "Delivery mode" = DELIVERYmode,
    "Gestational length (days)" = gestatlenght,
    Firstborn = firstborn
  ) %>%
  mutate(
    `Smoking during pregnancy` = ifelse(`Smoking during pregnancy` == "nee", FALSE, ifelse(`Smoking during pregnancy` == "ja", TRUE, NA)),
    `Alcohol during pregnancy` = ifelse(`Alcohol during pregnancy` == "nee", FALSE, ifelse(`Alcohol during pregnancy` == "ja", TRUE, NA)),
    `Delivery mode` = ifelse(`Delivery mode` == "natuurlijke bevalling", "vaginal", ifelse(`Delivery mode` == "pomp", "assisted vaginal", ifelse(grepl("keizer", `Delivery mode`), "cesarean section", NA)))
  ) %>%
  mutate(subject_id = as.character(subject_id))



# add age mb samples and ef to demogr:
load(here::here("rdata/ps_upto12.Rds"))
source("https://raw.githubusercontent.com/HenrikEckermann/in_use/master/mb_helper.R")
mbagedays <- sd_to_df(ps) %>% 
  select(subject_id, time, age_d) %>%
  filter(time < 4) %>%
  mutate(time = ifelse(time == 1, "Child age (days) at 28 days FMB collection", ifelse(time == 2, "Child age (days) at 75 days FMB collection", ifelse(time == 3, "Child age (days) at 105 days FMB collection", NA)))) %>%
  pivot_wider(names_from = time, values_from = age_d) %>%
  filter(subject_id %in% ids)

mbageyears <- sd_to_df(ps) %>% 
  select(subject_id, time, age_y) %>%
  filter(time %in% c(4, 5)) %>%
  mutate(time = ifelse(time == 4, "Child age (years) at 6 years FMB collection", ifelse(time == 5, "Child age (years) at 10 years FMB collection", NA))) %>%
  pivot_wider(names_from = time, values_from = age_y) %>%
  filter(subject_id %in% ids)

mbage <- full_join(mbageyears, mbagedays, by = "subject_id") %>%
  select(
    subject_id, 
    "Child age (days) at 28 days FMB collection",
    "Child age (days) at 75 days FMB collection",
    "Child age (days) at 105 days FMB collection",
    "Child age (years) at 6 years FMB collection",
    "Child age (years) at 10 years FMB collection"
  )


# double check median and IQR:
mbage %>% pivot_longer(-subject_id, names_to = "time", values_to = "age") %>%
    na.omit() %>%
    group_by(time) %>%
    summarise(
      med = median(age), 
      lower = quantile(age, 0.25), 
      upper = quantile(age, 0.75))

# add age ef    
load(here::here("rdata/age_ef.Rds"))
age_ef <- select(age_ef, subject_id, "Child age at BRIEF collection 8 years" = age8, "Child age at BRIEF and DS collection" = age10)

age_ef %>% filter(subject_id %in% ids) %>%
  select(-subject_id) %>%
  pivot_longer(everything(), names_to = "time") %>%
  group_by(time) %>%
  summarise(med = median(value), lower = quantile(value, 0.25), upper = quantile(value, 0.75))
demogr <- select(meta, subject_id, edu, Sex = sex) %>% 
  distinct() %>% 
  full_join(part1, by = "subject_id") %>% 
  full_join(mbage, by = "subject_id") %>%
  full_join(age_ef, by = "subject_id") %>%
  full_join(abx, by = "subject_id") %>%
  filter(subject_id %in% ids) %>%
  mutate("Maternal education" = 
        ifelse(edu == 1, "Primary education", 
        ifelse(edu == 2, "Seconday education", 
        ifelse(edu == 3, "Secondary education", 
        ifelse(edu == 4, "Secondary education", 
        ifelse(edu == 5, "Secondary education", 
        ifelse(edu == 6, "Secondary education", 
        ifelse(edu == 7, "College or University", 
        ifelse(edu == 8, "College or University", ifelse(edu == 9, "Other", NA))))))))),
        "Child sex" = ifelse(Sex == 1, "male", ifelse(Sex == 0, "female", NA))
      ) %>%
  select(-edu, -Sex, -subject_id)
demogr$`Maternal education` <- factor(
  demogr$`Maternal education`, 
  levels = c("Secondary education", "College or University")
)


save(demogr, file = here::here("rdata/demogr.Rds"))















#####################################################################
##########               supplementary tables              ##########
#####################################################################

load(file = here::here("rdata/16s_shannon_models.Rds"))

brms_tbl <- function(posterior, rename = NULL) {
  tab <- posterior %>%
    pivot_longer(everything(), names_to = "Parameter", values_to = "Value") %>%
    group_by(Parameter) %>%
    summarise(
      Estimate = median(Value),
      SD = sd(Value),
      lower = quantile(Value, 0.025),
      upper = quantile(Value, 0.975),
      .groups = "drop"
    ) %>%
    mutate(
      across(
        where(is.numeric), 
        function(x) format(round(x, 3), nsmall = 3))
      ) %>%
    mutate("95% CI" = glue("{lower} : {upper}")) %>%
    select(Parameter, Estimate, SD, "95% CI") %>%
    filter(
      str_detect(Parameter, "b[ABES]") |
      str_detect(Parameter, "^a_cog") |
      Parameter == "sigma"
    )
    # the parameter names used in stan are not informative, therefore I rename
    if (!is.null(rename)) {
      for (i in seq_along(rename)) {
        pair <- rename[[i]]
        key <- pair[1]
        value <- pair[2]
        tab <- mutate(tab, Parameter = str_replace(Parameter, key, value))
      }
    }
    tab
}

shannon_tabs1 <- map_dfr(shannon_models, function(lobj) {
  map_dfr(lobj, function(lobj_lowest) {
    tab <- brms_tbl(
      lobj_lowest$posterior, 
      rename = list(
      c("a_cog\\[1\\]", "intercept_female"), 
      c("a_cog\\[2\\]", "intercept_male"),
      c("bA", "age"),
      c("bB", "breastfeeding"),
      c("bE", "education"),
      c("bS", "shannon")
      )
    ) 
    tab$time <- lobj_lowest$time
    tab$y <- lobj_lowest$y
    tab
  }) 
})

# same for brief
load(file = here::here("rdata/16s_shannon_models_brief.Rds"))
shannon_tabs2 <- map_dfr(shannon_models_brief, function(lobj_lowest) {
    tab <- brms_tbl(
      lobj_lowest$posterior, 
      rename = list(
      c("a_cog\\[1\\]", "intercept_female"), 
      c("a_cog\\[2\\]", "intercept_male"),
      c("bA", "age"),
      c("bB", "breastfeeding"),
      c("bE", "education"),
      c("bS", "shannon"),
      c("a_cog", "intercept")
      )
    ) 
    tab$time <- lobj_lowest$time
    tab$y <- "BRIEF"
    tab
})


# merge 
model_tabs <- bind_rows(shannon_tabs1, shannon_tabs2) %>%
  mutate(y = factor(y, levels = c("total_f", "total_bw", "total_lns", "BRIEF"))) %>%
  arrange(time, y) %>%
  group_by(time, y) %>%
  nest()


# print header followed by table for each model 
if (!file.exists(here::here("article/table_supplement.Rmd"))) {
  file.create(here::here("Rmd/publication/model_tables.Rmd"))
  write_lines(
    '---\noutput:\n  word_document:\n    toc: true\n---\n\n', file = here::here("Rmd/publication/model_tables.Rmd"), append = FALSE)
 }

model_tabs$nr <- 1:20
pmap(model_tabs, function(time, y, data, nr) {
  y_pub <- ifelse(y == "BRIEF", "BRIEF", ifelse(y == "total_f", "Digit Span Forwards", ifelse(y == "total_bw", "Digit Span Backwards", ifelse(y == "total_lns", "Digit Span Letter Number Sequencing", NA))))
  time_pub <- ifelse(time == 28, "28 days", ifelse(time == 75, "75 days", ifelse(time == 105, "105 days", ifelse(time == 2193, "6 years", ifelse(time == 3655, "10 years", NA)))))
  write_lines(
    x = glue("\n\n### Table S{nr}: Parameter estimates for the model using microbiota samples obtained at {time_pub} as predictor and {y_pub} as outcome variable.\n \n"),
    file = here::here("Rmd/publication/model_tables.Rmd"), append = TRUE)
  write_lines(
    x = knitr::kable(data),
    file = here::here("Rmd/publication/model_tables.Rmd"), append = TRUE)  
})

rmarkdown::render(here::here("Rmd/publication/model_tables.Rmd"))
