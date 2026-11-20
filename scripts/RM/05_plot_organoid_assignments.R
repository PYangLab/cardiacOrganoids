#!/usr/bin/env Rscript


# Run section 1 once; repeat sections 2-4 for SC and SN after scoring each modality.


# Section 1. Setup ----
code_dir <- '.'
source(file.path(code_dir, 'AV_mapping_args.R'), local=TRUE)
# E.g. Rscript 05_plot_organoid_assignments.R --modality SC --t 0.10.
config <- parse_pipeline_args(defaults=list(
  data_dir='../data',
  out_dir='./results',
  gene_file='../data/Cardiomyocyte_atrial_and_ventricular_module_gene_sets.xlsx',
  modality='SC',
  t=0.10
))
print_pipeline_config(config)
source(file.path(code_dir, 'config_AV_mapping.R'), local=TRUE)
source(file.path(code_dir, 'AV_mapping_functions.R'), local=TRUE)

# Section 2. Choose modality and load its saved per-cell scores ----
# Interactive: set config$modality <- 'SN' before rerunning sections 2-4.
modality <- config$modality 
sources <- cfg$sources[[modality]]
datasets <- readLines(file.path(cfg$cache_dir, paste0(modality, '_organoid_datasets.txt')))
scores <- do.call(rbind, lapply(datasets, function(dataset)
  readRDS(file.path(cfg$cache_dir, 'organoid_scores', modality, paste0(dataset, '_scores.rds')))))

# Section 3. Assign A/V/unresolved and calculate proportions ----
# Atrial: A_score - V_score > t AND A_score > 0.
# Ventricular: V_score - A_score > t AND V_score > 0.
cfg$t <- config$t
if (cfg$t < 0) stop('--t must be nonnegative')
scores$assignment <- assign_av(scores$A_score, scores$V_score, cfg$t)
proportions <- organoid_proportions(scores, sources)
threshold_tag <- paste0('t', format(cfg$t, trim=TRUE, scientific=FALSE, nsmall=2))
saveRDS(scores, file.path(cfg$table_dir, paste0(modality, '_organoid_assignments_', threshold_tag, '.rds')))
write.csv(proportions,
          file.path(cfg$table_dir, paste0(modality, '_organoid_proportions_', threshold_tag, '.csv')),
          row.names=FALSE)

# Section 4. Figure 4d (SC) or Figure S4e (SN): stacked proportions ----
figure <- if (modality == 'SC') 'Figure_4d_CM_SC_AV' else 'Figure_S4e_CM_SN_AV'
plot_organoid_proportions(
  proportions, sources,
  output=file.path(cfg$figure_dir, paste0(figure, '_organoid_assignments_', threshold_tag, '.pdf'))
)

