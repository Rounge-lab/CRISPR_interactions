## plasmids, sex, region, age and antibiotic use 

require(tidyverse)
require(broom)
require(knitr)
require(dplyr)

setwd("PATH_TO_MANUS_FOLDER/results/")
metadata_plasmids <- read.delim("PATH_TO_MANUS_FOLDER/datasets/MAGs_AlphaDiv.tsv") 
metadata_plasmids <- metadata_plasmids %>% 
  left_join(metadata %>% rownames_to_column("sample_id") %>% dplyr::select(sample_id, reads_proc)) %>% 
  column_to_rownames("sample_id") %>% 
  mutate(age_cat = factor(age_cat, levels = c("50-59", "60-69", ">=70")))

## need reads_proc for this one

plasmids <- read.delim("PATH_TO_MANUS_FOLDER/datasets/PTUs_relab_renamed.tsv") %>% 
  column_to_rownames("sample_id")

source("PATH_TO_TOOLS/maaslin3-release/R/maaslin3.R")
source("PATH_TO_TOOLS/maaslin3-release/R/fit.R")
source("PATH_TO_TOOLS/maaslin3-release/R/utility_scripts.R")
source("PATH_TO_TOOLS/maaslin3-release/R/viz.R")

# 5) factor reference levels (optional, but wise)
metadata_plasmids$kjonn <- factor(metadata_plasmids$kjonn, levels = c("Female","Male"))
metadata_plasmids$senter <- factor(metadata_plasmids$senter, levels = c("Bærum","Moss"))
metadata_plasmids$age_cat <- factor(metadata_plasmids$age_cat, levels = c("50-59","60-69", ">=70"))
metadata_plasmids$beforeBL <- factor(metadata_plasmids$beforeBL, levels = c("No","Yes"))


### importing MAGS
MAGs_relab <- read.delim("PATH_TO_MANUS_FOLDER/datasets/MAGs_relab.tsv") %>% 
  column_to_rownames("sample_id")

## removing highly correlated features

plasmids
ids_extracted <- intersect(rownames(MAGs_relab), rownames(metadata_plasmids))
MAGs_relab <- MAGs_relab[ids_extracted, ]
plasmids <- plasmids[ids_extracted, ]
all(rownames(MAGs_relab) == rownames(metadata_plasmids))
all(rownames(plasmids) == rownames(metadata_plasmids))

identical(rownames(plasmids), rownames(MAGs_relab))
cor_matrix <- cor(plasmids, MAGs_relab, method = "spearman")

## removing hihgly correlated pathways
high_cor_paths <- colnames(plasmids)[apply(cor_matrix, 1, function(x) any(x > 0.5))]
plasmids_filtered <- plasmids[, !colnames(plasmids) %in% high_cor_paths]
plasmids <- plasmids_filtered

## remove highly correlated features

plasmids_diff_abund_filtered <-
  lapply(c("plasmids"), function(dataset) { 
    if (dataset == "plasmids") tmp_abund <- plasmids
    
    if (!dataset %in% list.dirs("PATH_TO_MANUS_FOLDER/results/diff_abund/", full.names = FALSE)) dir.create(paste0("PATH_TO_MANUS_FOLDER/results/diff_abund/"), recursive = TRUE)
    
    lapply(c("formula_sex", "formula_region", "formula_age", "formula_antibiotics"), function(formula) {
      if (formula == "formula_sex") tmp_formula <- "~ kjonn + reads_proc"
      if (formula == "formula_age") tmp_formula <- "~ age_cat + reads_proc"
      if (formula == "formula_region") tmp_formula <- "~ senter + reads_proc"
      if (formula == "formula_antibiotics") tmp_formula <- "~ beforeBL + reads_proc"
      
      tmp_maaslin <-
        maaslin3(
          input_data = tmp_abund,
          input_metadata = metadata_plasmids,
          transform = "LOG",
          normalization = "NONE",
          output = paste0("PATH_TO_MANUS_FOLDER/results/diff_abund/", "correlation_filtered", tmp_formula),
          formula = c(tmp_formula),
          plot_summary_plot = FALSE,
          plot_associations = FALSE,
          min_prevalence = 0.2
        )
      tmp_maaslin$results %>% 
        tibble() %>% 
        mutate(formula = formula,
               dataset = dataset)
    }) %>%
      bind_rows()
  }) %>% 
  bind_rows()


################# making volcano plot and colouring it by mobility and ARGs in antibiotics group

kjonn <- read.delim("PATH_TO_MANUS_FOLDER/results/diff_abund/correlation_filtered/correlation_filtered~ kjonn + reads_proc/all_results.tsv") %>% 
  mutate(dataset = "kjonn")
senter <- read.delim("PATH_TO_MANUS_FOLDER/results/diff_abund/correlation_filtered/correlation_filtered~ senter + reads_proc/all_results.tsv") %>% 
  mutate(dataset = "senter")
age <- read.delim("PATH_TO_MANUS_FOLDER/results/diff_abund/correlation_filtered/correlation_filtered~ age_cat + reads_proc/all_results.tsv")  %>% 
  mutate(dataset = "age")
antibitoics <- read.delim("PATH_TO_MANUS_FOLDER/results/diff_abund/correlation_filtered/correlation_filtered~ beforeBL + reads_proc/all_results.tsv") %>% 
  mutate(dataset = "beforeBL")
  
all_results <- kjonn %>% rbind(senter) %>% rbind(age) %>% rbind(antibitoics) %>% 
  mutate(feature = str_replace(feature, "CRCbiome.", "CRCbiome-")) %>% 
  dplyr::rename(., PTU = feature) %>% 
  left_join(read.delim("PATH_TO_MANUS_FOLDER/database_to_publish/PTUs_meta_lim1.csv")) %>% 
  filter(metadata != "reads_proc") %>% 
  dplyr::select(PTU, metadata, value, coef, pval_individual, qval_individual, pval_joint, model, predicted_mobility, NumARG_unique, dataset) %>% 
  mutate(negLog10P = -log10(as.numeric(qval_individual)),
    significance = case_when(
      as.numeric(qval_individual) < 0.05 ~ "Significant",
      TRUE ~ "Not Significant")) %>% 
  mutate(NumARG_unique = factor(NumARG_unique))


plot_prevalence <- all_results %>% filter(model == "abundance") %>% 
  #filter(qval_individual < 0.05) %>% 
  ggplot(., aes(
    x = coef,
    y = negLog10P,
    color = significance,
    shape = NumARG_unique
  )) +
  geom_point(size = 3, stroke = 0.7, position = position_jitter(width = 0.02, height = 0.02)) +
  scale_color_manual(values = c("gray", "orange")) +
  scale_shape_manual(values = c(21, 22, 18)) +  # circles and triangles
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    title = "Prevalence",
    x = "log2(Fold Change)",
    y = "-log10(p-value)"
  ) +
  facet_wrap(~ dataset) +
  theme_minimal(base_size = 14)
  

all_results %>% filter(dataset == "beforeBL" & model == "abundance") %>% 
  summarise(qval = sum(qval_individual > 0.05), 
            sign_conj = sum(qval_individual > 0.05 & predicted_mobility == "conjugative"),
                                 sign_mob = sum(qval_individual > 0.05 & predicted_mobility == "mobilizable"),
                                 sign_nonmob = sum(qval_individual > 0.05 & predicted_mobility == "non-mobilizable"),
                                 sign_nARG2 = sum(qval_individual > 0.05 & NumARG_unique == "2"),
                                 sign_nARG1 = sum(qval_individual > 0.05 & NumARG_unique == "1"),
                                 sign_nARG0 = sum(qval_individual > 0.05 & NumARG_unique == "0"))
