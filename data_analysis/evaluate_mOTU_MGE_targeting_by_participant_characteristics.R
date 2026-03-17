
# install.packages("box")
# library(broom)
box::use(
  readr[ read_csv, read_tsv, cols, read_rds, write_tsv],
  
  dplyr[ mutate, rename, select, distinct,
         left_join, group_by, summarize,
         filter, count, arrange,
         desc, bind_cols, rowwise,
         slice, bind_rows],
  
  tidyr[ pivot_longer, separate_longer_delim,
         pivot_wider, unnest],
  
  tidyselect[ everything],
  
  stringr[ str_remove, str_remove_all],
  
  tibble[ tibble],
  
  magrittr[ `%>%`],
  
  purrr[ map2],
  
  broom[ tidy],
  
  ggplot2[ ggplot, aes, geom_histogram,
           facet_wrap, geom_tile, 
           element_text, theme]
  
)

tmp <-
  rstudioapi::getActiveDocumentContext()$path %>% 
  str_remove("scripts/.*")

if (nchar(tmp) > 0) {
  if (!dir.exists(tmp)) {
    stop("Something's wrong with the paths.")
  } else {
    setwd(tmp)
  }
} else {
  stop("Need to save script")
}

rm(tmp)


## Script for evaluation of associations between mOTU-vOTU/PTU targeting and host factors

## 1. Get targeting data per mOTU

## Relevant vOTUs and PTUs
derep_table_vOTUs <- read_csv("PATH_TO_MANUS_FOLDER/datasets/viral_contigs_list_organized_onlymanus_votus.csv", col_types = cols()) %>% 
  mutate(genome_sample = str_remove(Contig, "_.*")) 
derep_table_PTUs <- read_tsv("PATH_TO_MANUS_FOLDER/datasets/PTU_rename_key.csv", col_types = cols()) %>% 
  mutate(genome_sample = str_remove(PTU, "_.*"))

## List of spacer targets
spacers_manus_table <- read_tsv("PATH_TO_MANUS_FOLDER/datasets/spacers_manus_table.csv", col_types = cols()) %>% 
  rename(sample_id = Sample) %>% 
  mutate(sample_id = str_replace(sample_id, "_", "-")) %>% 
  tibble() %>% 
  select(-Plasmids) %>% 
  distinct() %>% 
  select(-c(participant_id, sample_type, MAG, kjonn, age_invitation, senter, localization,final_result)) %>% 
  pivot_longer(c(vOTUs,PTUs), values_to = "targets", names_to = "target_type") %>% 
  mutate(targets = str_remove(targets, "^\\[") %>% 
           str_remove("\\]$") %>% 
           str_remove_all("\\'")) %>% 
  separate_longer_delim(cols = targets, delim = ", ") %>%
  mutate(manus_target = targets %in% c(derep_table_PTUs$PTU,
                                       derep_table_vOTUs$new_id)) %>% 
  mutate(targets = ifelse(manus_target, targets, NA)) %>% 
  distinct() %>% 
  left_join(read_tsv("PATH_TO_MANUS_FOLDER/datasets/contigs_by_genome_manus.tsv", col_types = cols()) %>% 
              select(Contig = contig, genome, MAG) %>% 
              distinct(), by = "Contig") 

## Get mOTU data
mOTU_data <- 
  read_tsv("PATH_TO_MANUS_FOLDER/datasets/contigs_by_genome_manus.tsv", col_types = cols()) %>% 
  group_by(MAG, genome) %>%
  summarize(n_contigs_per_genome = length(unique(contig))) %>% 
  group_by(MAG) %>% 
  summarize(n_genomes = length(unique(genome)),
            n_contigs = sum(n_contigs_per_genome),
            mean_contigs = mean(n_contigs_per_genome),
            sd_contigs = sd(n_contigs_per_genome),
            .groups = "drop") %>%
  ## Number of genomes with cassettes
  left_join(spacers_manus_table %>% 
              select(genome, MAG, CRISPR) %>% 
              distinct() %>% 
              group_by(MAG) %>% 
              summarize(n_cassettes = length(unique(CRISPR)),
                        n_cassette_genomes = length(unique(genome)), 
                        .groups = "drop"), by = "MAG")%>%
  mutate(across(c(n_cassettes, n_cassette_genomes), .fn = function(x) {
    ifelse(is.na(x), 0, x)
  }))
## Also need number of genomes with cassettes, and number of spacers per genome and in total

tmp <-
  spacers_manus_table %>% 
  filter(!is.na(targets)) %>% 
  filter(!is.na(genome)) %>% 
  filter(!is.na(MAG)) %>% 
  count(targets, genome, MAG) %>% 
  # left_join(mOTU_data %>% select(n_genomes, MAG, n_cassette_genomes, n_cassettes), by = "MAG") %>% 
  group_by(MAG, targets) %>% 
  summarize(n_targeting_genomes = n(),
            n_targeting_spacers = sum(n),
            median_targeting_spacers = median(n),
            min_targeting_spacers = min(n),
            max_targeting_spacers = max(n),
            mean_targeting_spacers = mean(n),
            sd_targeting_spacers = sd(n))

# tmp %>% 
#   filter(n_targeting_genomes > 2) %>% 
#   View()



vir_info <- read_csv(file = "PATH_TO_MANUS_FOLDER/datasets/vOTUs_taxonomy.csv")
MAG_taxonomy <- read_csv("PATH_TO_MANUS_FOLDER/datasets/MAG_taxonomy_full.tsv")


## spacers targeting MGE per mOTU - how many genomes targeting each MGE
mOTU_MGE_targeting <-
  spacers_manus_table %>% 
  filter(!is.na(targets)) %>% 
  select(targets, genome, MAG) %>% 
  distinct() %>% 
  filter(!is.na(genome), !is.na(MAG)) %>% 
  count(targets, MAG, name = "n_targeting_genomes") %>% 
  left_join(mOTU_data %>% select(n_genomes, MAG, n_cassette_genomes, n_cassettes), by = "MAG")

write_table <- FALSE

mOTU_MGE_targeting %>%
  mutate(fraction_genomes_targeting = n_targeting_genomes/n_cassette_genomes) %>% 
  # filter(n_cassette_genomes >= 20) %>%
  arrange(desc(fraction_genomes_targeting)) %>% 
  left_join(MAG_taxonomy %>% select(species, MAG), by = "MAG") %>% 
  left_join(vir_info %>% select(targets = Scaffold, Closer, Family), by = "targets") %>% 
  # summarize(n_single = sum(n_targeting_genomes == 1),
  #           n_tot = n(),
  #           frac_single = n_single/n_tot) %>% 
  (function(x) {
    if (write_table) {
      x %>% 
        write_tsv("PATH_TO_MANUS_FOLDER/datasets/targeting_genomes_per_target.tsv")
    } else {
       x
     }
  })

relevant_comparisons <-
  mOTU_MGE_targeting %>% 
  filter(n_targeting_genomes >= 50)

spacers_targeting_for_comparisons <-
  spacers_manus_table %>%
  filter(MAG %in% relevant_comparisons$MAG) %>% 
  select(sample_id, targets, MAG) %>% 
  distinct() %>% 
  mutate(targets = ifelse(is.na(targets), "none", targets))

participant_data <- 
  read_rds("PATH_TO_MANUS_FOLDER/participant_data/metadata_selected_variables.Rds") %>% 
  select(deltaker_id) %>% 
  bind_cols(read_rds("PATH_TO_MANUS_FOLDER/participant_data/metadata_selected_variables_cat.Rds") %>% select(-deltaker_id)) %>% 
  left_join(read_rds("PATH_TO_MANUS_FOLDER/participant_data/sample_meta.Rds") %>% select(sample_id, deltaker_id), by = "deltaker_id") %>% 
  mutate(sample_id = str_replace(sample_id, "_", "-")) %>% 
  select(sample_id, everything()) %>% 
  select(-deltaker_id)

host_char_targeting_tests <-
  relevant_comparisons %>% 
  rowwise() %>% 
  mutate(participant_associations = map2(.x = MAG, .y = targets, .f = function(testMAG = .x, testtarget = .y) {
    participant_dichotomization <-
      spacers_targeting_for_comparisons %>% 
      filter(MAG %in% testMAG) %>%
      group_by(sample_id) %>% 
      summarize(mOTU_MGE_targeting = any(targets == testtarget), .groups = "drop") %>% 
      left_join(participant_data, by = "sample_id")
    
    lapply(colnames(participant_dichotomization)[-c(1:2)], function(part_var) {
      tmp <- 
        participant_dichotomization %>% 
        select(mOTU_MGE_targeting, p_var = all_of(part_var)) 
      tmp %>%
        lm(as.integer(mOTU_MGE_targeting)~p_var, data = .) %>%
        tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
        mutate(variable = part_var) %>% 
        filter(!str_detect(term, "Intercept")) %>% 
        mutate(term = str_remove(term, "p_var")) %>% 
        left_join(tmp %>% 
                    count(mOTU_MGE_targeting, p_var) %>% 
                    pivot_wider(names_from = mOTU_MGE_targeting, values_from = n, values_fill = 0) %>% 
                    mutate(ref_level = p_var[1], ref_nontarget = `FALSE`[1], ref_target = `TRUE`[1]) %>% 
                    slice(-1) %>% 
                    select(term = p_var, ref_level, ref_nontarget, ref_target, test_nontarget = `FALSE`, test_target = `TRUE`),
                  by = "term")
    }) %>% 
      bind_rows()
  })) %>% 
  left_join(derep_table_PTUs %>% 
              select(targets = PTU, new_name = PTU_new), by = "targets") %>% 
  mutate(targets = ifelse(is.na(new_name), targets, new_name)) %>% 
  select(-new_name)

meta_vars <- read_rds("PATH_TO_MANUS_FOLDER/participant_data/metadata_variables.Rds")

host_char_targeting_tests %>% 
  unnest(participant_associations) %>% 
  mutate(MAG = stringr::str_replace(MAG, "MAG", "mOTU")) %>% 
  select(test_variable = variable, 
         ref_level, 
         test_level = term, 
         mOTU = MAG, 
         target = targets, 
         n_genomes, n_cassette_genomes,
         n_targeting_genomes,
         estimate, conf.low, conf.high, std.error, statistic, p.value) %>% 
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>% 
  mutate(across(c(estimate, std.error, statistic, p.value, conf.low, conf.high, p_adj), .fns = function(x) round(x, 3))) %>% 
  left_join(meta_vars %>% select(test_variable = var_id, var_name)) %>% 
  select(-test_variable) %>% 
  select(test_variable = var_name, everything()) %>% 
  write_tsv("PATH_TO_MANUS_FOLDER/results/mOTU_MGE_targeting_by_participant_characteristics.tsv")



