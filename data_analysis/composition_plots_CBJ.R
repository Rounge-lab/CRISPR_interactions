## making beta diversity plot for Katya, baseline descriptive
## only baseline samples



####### BRAY CURTIS
Combined_BrayCurtis <- read.delim("PATH_TO_MANUS_FOLDER/datasets/Combined_BrayCurtis.tsv") %>% 
  dplyr::rename(., "sample_id" = "X")
View(Combined_BrayCurtis)

MAGs_AlphaDiv <- read.delim("PATH_TO_MANUS_FOLDER/datasets/MAGs_AlphaDiv.tsv") %>% 
  column_to_rownames("sample_id")

set.seed(1234)
bray_dist <- Combined_BrayCurtis %>%
  column_to_rownames("sample_id") %>% 
  as.matrix() %>% 
  as.dist()

p_antibiotics <- 
  bray_dist %>% 
  dist_to_PCoA(group_var = Combined_BrayCurtis %>% 
                 dplyr::select(sample_id) %>% 
                 left_join(MAGs_AlphaDiv %>% rownames_to_column("sample_id") %>% dplyr::select(sample_id, beforeBL)) %>% 
                 pull(beforeBL)) %>% 
  plot_pcoa(dim_1 = "PCoA1", dim_2 = "PCoA2") +
  scale_color_manual(values = c("No" = "#9791c6", "Yes" = "#1c6462")) +
  scale_fill_manual(values = c("No" = "#9791c6", "Yes" = "#1c6462"))

p_sex <- 
  bray_dist %>% 
  dist_to_PCoA(group_var = Combined_BrayCurtis %>% 
                 dplyr::select(sample_id) %>% 
                 left_join(MAGs_AlphaDiv %>% rownames_to_column("sample_id") %>% dplyr::select(sample_id, kjonn)) %>% 
                 pull(kjonn)) %>% 
  plot_pcoa(dim_1 = "PCoA1", dim_2 = "PCoA2") +
  scale_color_manual(values = c("Female" = "#9791c6", "Male" = "#1c6462")) +
  scale_fill_manual(values = c("Female" = "#9791c6", "Male" = "#1c6462")) 

p_region <- 
  bray_dist %>% 
  dist_to_PCoA(group_var = Combined_BrayCurtis %>% 
                 dplyr::select(sample_id) %>% 
                 left_join(MAGs_AlphaDiv %>% rownames_to_column("sample_id") %>% dplyr::select(sample_id, senter)) %>% 
                 pull(senter)) %>% 
  plot_pcoa(dim_1 = "PCoA1", dim_2 = "PCoA2") +
  scale_color_manual(values = c("Bærum" = "#9791c6", "Moss" = "#1c6462")) +
  scale_fill_manual(values = c("Bærum" = "#9791c6", "Moss" = "#1c6462"))

p_age <- 
  bray_dist %>% 
  dist_to_PCoA(group_var = Combined_BrayCurtis %>% 
                 dplyr::select(sample_id) %>% 
                 left_join(MAGs_AlphaDiv %>% rownames_to_column("sample_id") %>% dplyr::select(sample_id, age_cat)) %>% 
                 pull(age_cat)) %>% 
  plot_pcoa(dim_1 = "PCoA1", dim_2 = "PCoA2") +
  scale_color_manual(values = c("50-59" = "#1c6462", "60-69" = "#9791c6", ">=70" = "#5fc0bf")) +
  scale_fill_manual(values = c("50-59" = "#1c6462", "60-69" = "#9791c6", ">=70" = "#5fc0bf"))  

# colours: "#b2b3cd" ace1a5 eac0bf
ggsave(plot = p_antibiotics, filename = "PATH_TO_MANUS_FOLDER/results/composition/comp_plot_antibiotics.pdf",
       unit="cm", width = 10, height = 10)
ggsave(plot = p_sex, filename = "PATH_TO_MANUS_FOLDER/results/composition/comp_plot_sex.pdf",
       unit="cm", width = 10, height = 10)
ggsave(plot = p_region, filename = "PATH_TO_MANUS_FOLDER/results/composition/comp_plot_region.pdf",
       unit="cm", width = 10, height = 10)
ggsave(plot = p_age, filename = "PATH_TO_MANUS_FOLDER/results/composition/comp_plot_age.pdf",
       unit="cm", width = 10, height = 10)



####### JACCARD 
jaccard <- read.delim("PATH_TO_MANUS_FOLDER/datasets/Combined_Jaccard.tsv") %>% 
  dplyr::rename(., "sample_id" = "X")

p_antibiotics_jac <- 
  jaccard %>% 
  dist_to_PCoA(group_var = Combined_Jaccard %>% 
                 dplyr::select(sample_id) %>% 
                 left_join(MAGs_AlphaDiv %>% rownames_to_column("sample_id") %>% dplyr::select(sample_id, beforeBL)) %>% 
                 pull(beforeBL)) %>% 
  plot_pcoa(dim_1 = "PCoA1", dim_2 = "PCoA2") +
  scale_color_manual(values = c("No" = "#9791c6", "Yes" = "#1c6462")) +
  scale_fill_manual(values = c("No" = "#9791c6", "Yes" = "#1c6462")) 

p_sex_jac <- 
  jaccard %>% 
  dist_to_PCoA(group_var = Combined_Jaccard %>% 
                 dplyr::select(sample_id) %>% 
                 left_join(MAGs_AlphaDiv %>% rownames_to_column("sample_id") %>% dplyr::select(sample_id, kjonn)) %>% 
                 pull(kjonn)) %>% 
  plot_pcoa(dim_1 = "PCoA1", dim_2 = "PCoA2") +
  scale_color_manual(values = c("Female" = "#9791c6", "Male" = "#1c6462")) +
  scale_fill_manual(values = c("Female" = "#9791c6", "Male" = "#1c6462")) 

p_region_jac <- 
  jaccard %>% 
  dist_to_PCoA(group_var = Combined_Jaccard %>% 
                 dplyr::select(sample_id) %>% 
                 left_join(MAGs_AlphaDiv %>% rownames_to_column("sample_id") %>% dplyr::select(sample_id, senter)) %>% 
                 pull(senter)) %>% 
  plot_pcoa(dim_1 = "PCoA1", dim_2 = "PCoA2") +
  scale_color_manual(values = c("Bærum" = "#9791c6", "Moss" = "#1c6462")) +
  scale_fill_manual(values = c("Bærum" = "#9791c6", "Moss" = "#1c6462"))

p_age_jac <- 
  jaccard %>% 
  dist_to_PCoA(group_var = Combined_Jaccard %>% 
                 dplyr::select(sample_id) %>% 
                 left_join(MAGs_AlphaDiv %>% rownames_to_column("sample_id") %>% dplyr::select(sample_id, age_cat)) %>% 
                 pull(age_cat)) %>% 
  plot_pcoa(dim_1 = "PCoA1", dim_2 = "PCoA2") +
  scale_color_manual(values = c("50-59" = "#1c6462", "60-69" = "#9791c6", ">=70" = "#5fc0bf")) +
  scale_fill_manual(values = c("50-59" = "#1c6462", "60-69" = "#9791c6", ">=70" = "#5fc0bf")) 

# colours: "#b2b3cd" ace1a5 eac0bf
ggsave(plot = p_antibiotics_jac, filename = "PATH_TO_MANUS_FOLDER/results/composition/comp_plot_antibiotics_jaccard.pdf",
       unit="cm", width = 10, height = 10)
ggsave(plot = p_sex_jac, filename = "PATH_TO_MANUS_FOLDER/results/composition/comp_plot_sex_jaccard.pdf",
       unit="cm", width = 10, height = 10)
ggsave(plot = p_region_jac, filename = "PATH_TO_MANUS_FOLDER/results/composition/comp_plot_region_jaccard.pdf",
       unit="cm", width = 10, height = 10)
ggsave(plot = p_age_jac, filename = "PATH_TO_MANUS_FOLDER/results/composition/comp_plot_age_jaccard.pdf",
       unit="cm", width = 10, height = 10)


