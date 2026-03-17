### JACCARD PEMRANOVA with lifestyle and diet
setwd("PATH_TO_MANUS_FOLDER") 

jaccard_interactions <- 
  read.delim("datasets/Individual_interactions_Jaccard_CrisprF.tsv") %>% 
  column_to_rownames("X") %>% 
  as.matrix() %>% 
  as.dist()

jaccard_data <- 
  read.delim("datasets/Individual_interactions_Jaccard_CrisprF.tsv") %>% 
  dplyr::rename("sample_id" = "X")

plasmids <- read.delim("/datasets/PTUs_Jaccard.tsv")

MAGs_AlphaDiv <- read.delim("/datasets/MAGs_AlphaDiv.tsv")

crcbiome_dmm <- read.delim("/datasets/crcbiome_dmm.tsv")

ids_matrix <- jaccard_interactions %>% as.matrix %>% as.data.frame %>% colnames()
metadata_imputed <- metadata[ids_matrix, , drop = FALSE]
metadata_imputed <- metadata_imputed %>% rownames_to_column("sample_id") %>% 
  left_join(MAGs_AlphaDiv %>% dplyr::select(sample_id, beforeBL)) %>% 
  left_join(crcbiome_dmm) %>% 
  column_to_rownames("sample_id")

metadata_imputed_all <- metadata %>% rownames_to_column("sample_id") %>% 
  left_join(MAGs_AlphaDiv %>% dplyr::select(sample_id, beforeBL)) %>% 
  left_join(crcbiome_dmm) %>% 
  column_to_rownames("sample_id")

################################ ANALYSES USED FOR FIGURES IN PAPER, the rest is development

#### metadata is performed similar to the analyses_all_V1 script from the kjonn/scripts folder
### variables to be imputed
group_colors <- c(
  "#9791c6",  
  "#1c6462",  
  "#67A9CF",
  "#FFD55C")  

variables_cat <- 
  c("Smoking", 
    "Snus",
    "Utdanning", 
    "Sivilstatus_cat2",
    "Arbeid_lump",
    "Nasj_cat2",
    "wcrf_index_main", 
    "Antibiotics",
    "Antacids") 

variables_cont <- 
  c("Energi_kcal", 
    "Prot_energi", 
    "Karboh_energi", 
    "Sukker_energi", 
    "Fiber", 
    "Fett_energi", 
    "Mettet_energi", 
    "C_enum_energi",
    "C_flerum_energi", 
    "Trans_u_energi", 
    "Alko",
    "BMI",
    "PhysAct_Score")
## need to impute the data 
table(is.na(metadata_imputed))
metadata_imputed <- metadata_imputed
for (v in variables_cont) {
  metadata_imputed[[v]][is.na(metadata_imputed[[v]])] <- median(metadata_imputed[[v]], na.rm = TRUE)
}
table(is.na(metadata_imputed))
metadata_imputed <- metadata_imputed %>%
  mutate(across(
    all_of(variables_cat),
    ~ replace_na(as.character(.x), "Missing") |> factor()
  ))
table(is.na(metadata_imputed))
#### same for the other metadata table, ideally here one would just subset the metadata table all when needed
table(is.na(metadata_imputed_all))
metadata_imputed_all <- metadata_imputed_all
for (v in variables_cont) {
  metadata_imputed_all[[v]][is.na(metadata_imputed_all[[v]])] <- median(metadata_imputed_all[[v]], na.rm = TRUE)
}
table(is.na(metadata_imputed_all))
metadata_imputed_all <- metadata_imputed_all %>%
  mutate(across(
    all_of(variables_cat),
    ~ replace_na(as.character(.x), "Missing") |> factor()
  ))
table(is.na(metadata_imputed_all))
# Example: edit to match your column names in metadata
blocks <- list(
  diet = c("Energi_kcal", "Prot_energi", "Karboh_energi", "Sukker_energi", "Fiber", 
           "Fett_energi", "Mettet_energi", "C_enum_energi", "C_flerum_energi", "Trans_u_energi", "Alko"),
  lifestyle = c("BMI", "PhysAct_Score", "Smoking", "Snus", "wcrf_index_main", "beforeBL", "Antacids"),
  demographics = c("kjonn", "senter", "age_invitation", "Utdanning", "Sivilstatus_cat2",
                   "Arbeid_lump", "Nasj_cat2"),
  dmm = c("gr"),
  all =  c("Energi_kcal", "Prot_energi", "Karboh_energi", "Sukker_energi", "Fiber", "Fett_energi", "Mettet_energi", "C_enum_energi", "C_flerum_energi", "Trans_u_energi", 
           "Alko", "BMI", "PhysAct_Score", "Smoking", "Snus", "Utdanning", "Sivilstatus_cat2",
           "Arbeid_lump", "Nasj_cat2", "wcrf_index_main", "beforeBL", "Antacids", "kjonn", "senter", "age_invitation"))
#### barplot of all variables. R2 as x axis, and * for singificance, coloured by grouping ie lifestyle, diet etc
diet_vars <- c("Energi_kcal", "Prot_energi", "Karboh_energi", "Sukker_energi", "Fiber", "Fett_energi", "Mettet_energi", "C_enum_energi", "C_flerum_energi", "Trans_u_energi", "Alko")
lifestyle_vars <- c("BMI", "PhysAct_Score", "Smoking", "Snus", "wcrf_index_main", "beforeBL", "Antacids")
demographics_vars <- c("kjonn", "senter", "age_invitation", "Utdanning", "Sivilstatus_cat2", "Arbeid_lump", "Nasj_cat2")
dmm <- c("gr")

###### making barplot for votus and mags as well as interactions side by side
vOTUs_Jaccard <- read.delim("datasets/vOTUs_Jaccard.tsv") %>% 
  column_to_rownames("X") %>% 
  as.matrix() %>% 
  as.dist()

MAGs_Jaccard <- read.delim("datasets/MAGs_Jaccard.tsv") %>% 
  column_to_rownames("X") %>% 
  as.matrix() %>% 
  as.dist()

PTUs_Jaccard <- read.delim("datasets/PTUs_Jaccard.tsv") %>% 
  column_to_rownames("X") %>% 
  as.matrix() %>% 
  as.dist()

metadata_votus <- metadata_imputed_all %>% rownames_to_column("sample_id") %>% 
  mutate(sample_id = str_replace(sample_id, "S-", "S_")) %>% 
  left_join(MAGs_AlphaDiv %>% dplyr::select(sample_id, beforeBL)) %>% 
  left_join(crcbiome_dmm) %>% 
  column_to_rownames("sample_id")

#############################
run_permanova <- function(dist_mat,
                                     metadata,
                                     vars,          # variables to test, character vector
                                     covariates = NULL,
                                     permutations = 999) {
  dist <- dist_mat
  
  res_list <- lapply(vars, function(v) {
    # RHS: this variable + covariates
    rhs_terms <- c(v, covariates)
    rhs <- paste(rhs_terms, collapse = " + ")
    
    form <- as.formula(paste("dist ~", rhs))
    
    fit <- adonis2(
      formula      = form,
      data         = metadata,
      permutations = permutations,
      by           = "margin"   # so v is tested while controlling for covariates
    )
    
    # extract the row corresponding to the variable of interest
    row <- fit[rownames(fit) == v, , drop = FALSE]
    
    data.frame(
      variable = v,
      df       = row$Df,
      R2       = row$R2,
      F        = row$F,
      p        = row$`Pr(>F)`,
      row.names = NULL
    )
  })
  
  results <- do.call(rbind, res_list)
  
  # FDR correction (Benjamini–Hochberg)
  results$p_fdr <- p.adjust(results$p, method = "fdr")
  
  results
}

vars_to_test <- c("Energi_kcal", "Prot_energi", "Karboh_energi", "Sukker_energi", "Fiber", "Fett_energi", "Mettet_energi", "C_enum_energi", 
                  "C_flerum_energi", "Trans_u_energi", "Alko", "BMI", "PhysAct_Score", "Smoking", "Snus", "wcrf_index_main", "beforeBL", 
                  "Antacids", "kjonn", "senter", "age_invitation", "Utdanning", 
                  "Sivilstatus_cat2", "Arbeid_lump", "Nasj_cat2", "gr")

covariates   <- c("Total_Bases_QC_ATLAS", "gr")

results_int <- run_permanova(
  dist_mat   = jaccard_interactions,
  metadata   = metadata_imputed,
  vars       = vars_to_test,
  covariates = covariates,
  permutations = 999
)

results_int

permanova_barplot_int <- results_int %>%
  mutate(group = case_when(
    variable %in% lifestyle_vars  ~ "Lifestyle",
    variable %in% diet_vars       ~ "Diet",
    variable %in% demographics_vars ~ "Demography",
    TRUE                           ~ "Other")) %>% 
  mutate(sig = case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE          ~ ""))%>% 
  mutate(variable = case_when(
    (variable == "Energi_kcal") ~ "Energy (kcal/day)",
    (variable == "Prot_energi") ~ "Proteins (%)", 
    (variable == "Karboh_energi") ~ "Carbohydrates (%)", 
    (variable == "Sukker_energi") ~ "Added sugar (%)",
    (variable == "Fett_energi") ~ "Fats (%)",
    (variable == "Mettet_energi") ~ "SFA (%)",
    (variable == "C_enum_energi") ~ "MUFA (%)", 
    (variable == "C_flerum_energi") ~ "PUFA (%)", 
    (variable == "Trans_u_energi") ~ "TFA (%)", 
    (variable == "Alko") ~ "Alcohol (g/day)",
    (variable == "Fiber") ~ "Fiber (g/day)",
    (variable == "PhysAct_Score") ~ "Physical activity", 
    (variable == "wcrf_index_main") ~ "HLI", 
    (variable == "beforeBL") ~ "Antibiotics", 
    (variable == "kjonn") ~ "Sex", 
    (variable == "senter") ~ "Screening region", 
    (variable == "age_invitation") ~ "Age", 
    (variable == "Utdanning") ~ "Education", 
    (variable == "Sivilstatus_cat2") ~ "Marital status", 
    (variable == "Arbeid_lump") ~ "Employment status", 
    (variable == "Nasj_cat2") ~ "Nationality",
    TRUE ~ variable))

(plot_perm_int <- permanova_barplot_int %>%
    filter(variable != "gr") %>%
    group_by(group) %>%
    mutate(variable = fct_reorder(variable, R2)) %>%
    ggplot(., aes(x = R2, y = variable, fill = group)) +
    geom_col() +
    geom_text(aes(label = sig),
              hjust = -0.2,        # push star a bit to the right of the bar
              size = 4) +
    #facet_wrap(~ group, scales = "free_y") +
    scale_x_continuous(limits = c(0, 0.025), breaks = seq(0, 0.025, by = 0.01),
                       expand = c(0, 0)) +
    scale_fill_manual(values = c(group_colors)) +
    labs(x = "R² int (PERMANOVA)",
         y = NULL,
         fill = "Variable group") +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.text  = element_text(size = 9),
          legend.title = element_text(size = 9),
          axis.title.x = element_text(size = 9)))

results_votus <- run_permanova(
  dist_mat   = vOTUs_Jaccard,
  metadata   = metadata_votus,
  vars       = vars_to_test,
  covariates = covariates,
  permutations = 999
)

results_votus

permanova_barplot_votus <- results_votus %>%
  mutate(group = case_when(
    variable %in% lifestyle_vars  ~ "Lifestyle",
    variable %in% diet_vars       ~ "Diet",
    variable %in% demographics_vars ~ "Demography",
    TRUE                           ~ "Other")) %>% 
  mutate(sig = case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE          ~ ""))%>% 
  mutate(variable = case_when(
    (variable == "Energi_kcal") ~ "Energy (kcal/day)",
    (variable == "Prot_energi") ~ "Proteins (%)", 
    (variable == "Karboh_energi") ~ "Carbohydrates (%)", 
    (variable == "Sukker_energi") ~ "Added sugar (%)",
    (variable == "Fett_energi") ~ "Fats (%)",
    (variable == "Mettet_energi") ~ "SFA (%)",
    (variable == "C_enum_energi") ~ "MUFA (%)", 
    (variable == "C_flerum_energi") ~ "PUFA (%)", 
    (variable == "Trans_u_energi") ~ "TFA (%)", 
    (variable == "Alko") ~ "Alcohol (g/day)",
    (variable == "Fiber") ~ "Fiber (g/day)",
    (variable == "PhysAct_Score") ~ "Physical activity", 
    (variable == "wcrf_index_main") ~ "HLI", 
    (variable == "beforeBL") ~ "Antibiotics", 
    (variable == "kjonn") ~ "Sex", 
    (variable == "senter") ~ "Screening region", 
    (variable == "age_invitation") ~ "Age", 
    (variable == "Utdanning") ~ "Education", 
    (variable == "Sivilstatus_cat2") ~ "Marital status", 
    (variable == "Arbeid_lump") ~ "Employment status", 
    (variable == "Nasj_cat2") ~ "Nationality",
    TRUE ~ variable))

(plot_votus <- permanova_barplot_votus %>%
    filter(variable != "gr") %>%
    group_by(group) %>%
    mutate(variable = fct_reorder(variable, R2)) %>%
    ggplot(., aes(x = R2, y = variable, fill = group)) +
    geom_col() +
    geom_text(aes(label = sig),
              hjust = -0.2,        # push star a bit to the right of the bar
              size = 4) +
    #facet_wrap(~ group, scales = "free_y") +
    scale_x_continuous(limits = c(0, 0.025), breaks = seq(0, 0.025, by = 0.01),
                       expand = c(0, 0)) +
    scale_fill_manual(values = c(group_colors)) +
    labs(x = "R² vOTU (PERMANOVA)",
         y = NULL,
         fill = "Variable group") +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.text  = element_text(size = 9),
          legend.title = element_text(size = 9),
          axis.title.x = element_text(size = 9)))

###############
results_mags <- run_permanova(
  dist_mat   = MAGs_Jaccard,
  metadata   = metadata_votus,
  vars       = vars_to_test,
  covariates = covariates,
  permutations = 999
)


permanova_barplot_mags <- results_mags %>%
  mutate(group = case_when(
    variable %in% lifestyle_vars  ~ "Lifestyle",
    variable %in% diet_vars       ~ "Diet",
    variable %in% demographics_vars ~ "Demography",
    TRUE                           ~ "Other")) %>% 
  mutate(sig = case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE          ~ ""))%>% 
  mutate(variable = case_when(
    (variable == "Energi_kcal") ~ "Energy (kcal/day)",
    (variable == "Prot_energi") ~ "Proteins (%)", 
    (variable == "Karboh_energi") ~ "Carbohydrates (%)", 
    (variable == "Sukker_energi") ~ "Added sugar (%)",
    (variable == "Fett_energi") ~ "Fats (%)",
    (variable == "Mettet_energi") ~ "SFA (%)",
    (variable == "C_enum_energi") ~ "MUFA (%)", 
    (variable == "C_flerum_energi") ~ "PUFA (%)", 
    (variable == "Trans_u_energi") ~ "TFA (%)", 
    (variable == "Alko") ~ "Alcohol (g/day)",
    (variable == "Fiber") ~ "Fiber (g/day)",
    (variable == "PhysAct_Score") ~ "Physical activity", 
    (variable == "wcrf_index_main") ~ "HLI", 
    (variable == "beforeBL") ~ "Antibiotics", 
    (variable == "kjonn") ~ "Sex", 
    (variable == "senter") ~ "Screening region", 
    (variable == "age_invitation") ~ "Age", 
    (variable == "Utdanning") ~ "Education", 
    (variable == "Sivilstatus_cat2") ~ "Marital status", 
    (variable == "Arbeid_lump") ~ "Employment status", 
    (variable == "Nasj_cat2") ~ "Nationality",
    TRUE ~ variable))

(plot_mags <- permanova_barplot_mags %>%
    filter(variable != "gr") %>%
    group_by(group) %>%
    mutate(variable = fct_reorder(variable, R2)) %>%
    ggplot(., aes(x = R2, y = variable, fill = group)) +
    geom_col() +
    geom_text(aes(label = sig),
              hjust = -0.2,        # push star a bit to the right of the bar
              size = 4) +
    #facet_wrap(~ group, scales = "free_y") +
    scale_x_continuous(limits = c(0, 0.025), breaks = seq(0, 0.025, by = 0.01),
      expand = c(0, 0)) +
    scale_fill_manual(values = c(group_colors)) +
    labs(x = "R² MAGS (PERMANOVA)",
         y = NULL,
         fill = "Variable group") +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.text  = element_text(size = 9),
          legend.title = element_text(size = 9),
          axis.title.x = element_text(size = 9)))


###############################


results_PTU <- run_permanova(
  dist_mat   = PTUs_Jaccard,
  metadata   = metadata_votus,
  vars       = vars_to_test,
  covariates = covariates,
  permutations = 999
)

results_PTU

permanova_barplot_PTUs <- results_PTU %>%
  mutate(group = case_when(
    variable %in% lifestyle_vars  ~ "Lifestyle",
    variable %in% diet_vars       ~ "Diet",
    variable %in% demographics_vars ~ "Demography",
    TRUE                           ~ "Other")) %>% 
  mutate(sig = case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE          ~ ""))%>% 
  mutate(variable = case_when(
    (variable == "Energi_kcal") ~ "Energy (kcal/day)",
    (variable == "Prot_energi") ~ "Proteins (%)", 
    (variable == "Karboh_energi") ~ "Carbohydrates (%)", 
    (variable == "Sukker_energi") ~ "Added sugar (%)",
    (variable == "Fett_energi") ~ "Fats (%)",
    (variable == "Mettet_energi") ~ "SFA (%)",
    (variable == "C_enum_energi") ~ "MUFA (%)", 
    (variable == "C_flerum_energi") ~ "PUFA (%)", 
    (variable == "Trans_u_energi") ~ "TFA (%)", 
    (variable == "Alko") ~ "Alcohol (g/day)",
    (variable == "Fiber") ~ "Fiber (g/day)",
    (variable == "PhysAct_Score") ~ "Physical activity", 
    (variable == "wcrf_index_main") ~ "HLI", 
    (variable == "beforeBL") ~ "Antibiotics", 
    (variable == "kjonn") ~ "Sex", 
    (variable == "senter") ~ "Screening region", 
    (variable == "age_invitation") ~ "Age", 
    (variable == "Utdanning") ~ "Education", 
    (variable == "Sivilstatus_cat2") ~ "Marital status", 
    (variable == "Arbeid_lump") ~ "Employment status", 
    (variable == "Nasj_cat2") ~ "Nationality",
    TRUE ~ variable))

(plot_PTUs <- permanova_barplot_PTUs %>%
    filter(variable != "gr") %>%
    group_by(group) %>%
    mutate(variable = fct_reorder(variable, R2)) %>%
    ggplot(., aes(x = R2, y = variable, fill = group)) +
    geom_col() +
    geom_text(aes(label = sig),
              hjust = -0.2,        # push star a bit to the right of the bar
              size = 4) +
    #facet_wrap(~ group, scales = "free_y") +
    scale_x_continuous(limits = c(0, 0.025), breaks = seq(0, 0.025, by = 0.01),
                       expand = c(0, 0)) +
    scale_fill_manual(values = c(group_colors)) +
    labs(x = "R² PTU (PERMANOVA)",
         y = NULL,
         fill = "Variable group") +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.text  = element_text(size = 9),
          legend.title = element_text(size = 9), 
          axis.title.x = element_text(size = 9)))

plot_PTUs 
plot_perm_int
plot_mags
plot_votus


ggsave(filename = "results/permanova_barplot_seqdmmadj_all.pdf",
       plot = ggarrange(plot_perm_int, plot_mags, plot_votus, plot_PTUs, ncol=4, common.legend = T), 
       width = 25, height = 10, unit = "cm")










##################################



run_block_permanova <- function(dist_mat, metadata, blocks, covariates = NULL,
                                permutations = 999) {
  # make sure distance object is called 'dist' for the formula
  dist <- dist_mat
  
  res_list <- lapply(names(blocks), function(block_name) {
    vars <- blocks[[block_name]]
    
    # RHS of formula: block vars + optional covariates
    rhs_terms <- c(vars, if (!is.null(covariates)) covariates)
    rhs <- paste(rhs_terms, collapse = " + ")
    
    form <- as.formula(paste("dist ~", rhs))
    
    fit <- adonis2(
      formula      = form,
      data         = metadata,
      permutations = permutations,
      by           = NULL   # <-- overall model test (whole block)
    )
    
    # For by = NULL, first row = whole model
    data.frame(
      block = block_name,
      df    = fit$Df[1],
      R2    = fit$R2[1],
      F     = fit$F[1],
      p     = fit$`Pr(>F)`[1],
      row.names = NULL
    )
  })
  
  do.call(rbind, res_list)
}


results_covariates <- run_block_permanova(
  dist_mat   = jaccard_interactions,
  metadata   = metadata_imputed,
  blocks     = blocks,
  covariates = covariates,  # or NULL if you don’t want adjustment
  permutations = 999
)




