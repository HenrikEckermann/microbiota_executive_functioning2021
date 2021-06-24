library(phyloseq)
library(tidyverse)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(microbiome)
library(qvalue)
library(DirichletMultinomial)

# some helper functions to switch between pseq and dataframe formats 
source(here::here("R/mb_helper.R"))

# import data (bibo is a pseq object)
bibo <- readRDS(file = file.path(here::here("rdata/yang_bibo_meta_orclu_div.rds")))
bibo <- subset_samples(bibo, sample_data(bibo)$Age %in% c("6y", "10y"))
bibo <- prune_taxa(taxa_sums(bibo)>0, bibo)
bibo <- tax_glom(bibo, taxrank = rank_names(bibo)[6])
bibo <- aggregate_taxa(bibo, level = "Genus", verbose = F)
bibo <- microbiome::transform(bibo, transform = "clr")

# extract the dmm clusters from the metadata 
ydmm_ch <- sd_to_df(bibo) %>% 
  select(subject_id = ID, sample_id = "age.and.ID", dmm = Clusters) %>%
  mutate(
    dmm = str_extract(dmm, "\\d"),
    dmm = as.factor(dmm)
  )

# create a df where we have all taxa as columns as well as a dmm column 
otus <- colnames(otu_to_df(bibo))[-1]
df_temp <- sd_to_df(bibo) %>% select(sample_id, Age)
df <- otu_to_df(bibo) %>% 
  full_join(ydmm_ch, by = "sample_id") %>%
  full_join(df_temp, by = "sample_id")

topn <- 50
# load the dmm fit objects
load(here::here("rdata/dmm_split.Rds"))
mbest <- p_lapl$childhood$bestfit
m0 <- p_lapl$childhood$dmmfit[[1]]
p0 <- fitted(m0, scale = TRUE)     # scale by theta
p4 <- fitted(mbest, scale = TRUE)

# the difference between p0 and p4 will reflect how each cluster 
# differs to the whole population average whereas the highest diff values 
# reflect which genera are most different to the average
diff_per_p <- as_tibble(p0, rownames = NA) %>%
  rownames_to_column("taxid") %>%
  full_join(
    as_tibble(p4, rownames = NA) %>% rownames_to_column("taxid"),
    by = "taxid"
  ) %>%
  select(taxid, p0 = V1.x, k1 = V1.y, k2 = V2, k3 = V3, k4 = V4) %>%
  pivot_longer(contains("k"), names_to = "k") %>%
  mutate(diff = abs(value - p0)) %>%
  group_by(k) %>% nest()

toptaxa <- map(diff_per_p$data, function(data) {
  data %>% arrange(desc(diff)) %>%
    head(topn) %>%
    .$taxid
})
topnunion <- Reduce(union, toptaxa)
load(here::here("rdata/taxtable.Rds"))
toptaxa <- taxtable %>% filter(taxid %in% topnunion) %>%
  mutate(Genus = ifelse(
    Genus == "g__", glue("{Domain}_{Phylum}_{Class}_{Order}_{Family}_Unidentified_Genus"), 
      str_replace(Genus, "g__", ""))) %>%
  mutate(Genus = str_replace_all(Genus, "\\w__", ""))
taxtable_y <- tax_table(bibo) %>% as.data.frame()
topnunion <- taxtable_y %>% filter(unique %in% toptaxa$Genus) %>%
  .$unique


# now standardize each taxon and prepare for the plotting
df <- mutate(df, across(all_of(otus), function(x) scale(x)[, 1])) %>%
  as.data.frame() %>%
  column_to_rownames("sample_id") %>%
  arrange(dmm)

# create within cluster sorting of samples 
sample_sort <- map(1:4, function(k) {
  mat_temp <- df %>%
                filter(dmm == k) %>%
                select(-Age, -dmm) %>% 
                as.matrix()
  mat_temp <- mat_temp[, topnunion]
  sample_sort <- neatsort(mat_temp, target = "rows", method = "NMDS", distance = "euclidean")
})
# genus sort for all samples
mat_temp <- select(df, -Age, -dmm) %>% as.matrix()
mat_temp <- mat_temp[, topnunion]
genus_sort <- neatsort(mat_temp, target = "cols", method = "NMDS", distance = "euclidean")
rown_df <- tibble(sampleid = unlist(sample_sort))
df <- rown_df %>% 
  left_join(df %>% rownames_to_column("sampleid"), by = "sampleid") %>%
  column_to_rownames("sampleid")
sample_id <- rownames(df)
clu <- df$dmm
age <- df$Age
df_fin <- df[, genus_sort]
df_fin <- t(df_fin)
# sort only based on small matrix 
colnames(df_fin) <- clu


# fix rownames 
newnames <- str_match(rownames(df_fin), "(.*_)(\\w+_Unidentified_Genus)") %>%
  as_tibble() %>%
  mutate(fullname = rownames(df_fin)) %>%
  mutate(newname = ifelse(is.na(V1), fullname, V3)) %>%
  .$newname
row.names(df_fin) <- newnames

mat <- df_fin

# legend color
nvec <- seq(-1, 1, length.out = 8)
col_fun = colorRamp2(nvec, 
                     c(
                       "#2166ac", 
                       "#4393c3", 
                       "#92c5de", 
                       "#d1e5f0", 
                       "#fddbc7", 
                       "#f4a582", 
                       "#d6604d", 
                       "#b2182b"
                   ))


# annotation
ha <- HeatmapAnnotation(
  Age = factor(age, levels = c("6y", "10y")),
  Clusters = clu,
  col = list(Age = c("6y"="#78c679",
                     "10y"="#006837"),
             Clusters = c("1"="#a6cee3",
                          "2"="#1f78b4",
                          "3"="#b2df8a",
                          "4"="#33a02c")),
  annotation_name_side = "left"
  #show_legend = c("Age" = FALSE, "Clusters" = FALSE)

                      )


# Heatmap
png(here::here("fig/top50union.png"), 
  width = 50, 
  height = 34, 
  units = "cm", 
  res= 900)
  Heatmap(mat,  
          width = unit(32, "cm"),
          height = unit(32, "cm"),
          top_annotation = ha,
          name = "RA", 
          col = col_fun,
          column_title = "",
          column_title_gp = gpar(fontsize = 18, fontface = "bold"),
          row_title = "",
          row_title_side = "right",
          row_title_gp = gpar(fontsize = 18, fontface = "bold"),
          row_names_gp = gpar(fontsize = 15, fontface = "italic"),
          row_names_side = "left",
          show_column_names = F,
          border = T,
          cluster_rows = F,
          cluster_columns = F,
          show_heatmap_legend = TRUE,
          heatmap_legend_param = list(ncol =1)
        )
  dev.off()


