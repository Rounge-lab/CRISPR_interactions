library(tidyverse)

## get contigs per genome
atlas_path <- "PATH_TO_ATLAS"
samples <- read_tsv("PATH_TO_MANUS_FOLDER/datasets/data_by_samples_crispr.csv", col_types = cols()) %>% 
  select(sample_id) %>% 
  mutate(sample_id = str_replace(sample_id, "_", "-"))


sample_wise_contigs_by_genome <-
  samples %>% 
  rowwise() %>% 
  mutate(contigs_for_genomes = map(.x = sample_id, .f = function(sample_id = .x) {
    read_tsv(glue::glue("{atlas_path}/{sample_id}/binning/DASTool/cluster_attribution.tsv"), 
             col_types = cols(), col_names = FALSE)
  })) %>% 
  unnest(contigs_for_genomes) %>% 
  select(sample_id, contig = X1, genome = X2) %>% 
  left_join(read_tsv(glue::glue("{atlas_path}/genomes/clustering/allbins2genome.tsv"), 
                     col_types = cols()), by = "genome")

sample_wise_contigs_by_genome %>% 
  write_tsv("PATH_TO_MANUS_FOLDER/datasets/contigs_by_genome_manus.tsv")
