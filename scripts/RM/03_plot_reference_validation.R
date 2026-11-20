#!/usr/bin/env Rscript


# Requires the three completed reference runs from 02_score_references.R.

# Section 1. Setup ----
code_dir <- '.'
source(file.path(code_dir, 'AV_mapping_args.R'), local=TRUE)
config <- parse_pipeline_args(defaults=list(
  data_dir='../data',
  out_dir='./results',
  gene_file='../data/Cardiomyocyte_atrial_and_ventricular_module_gene_sets.xlsx'
))
print_pipeline_config(config)
source(file.path(code_dir, 'config_AV_mapping.R'), local=TRUE)
source(file.path(code_dir, 'AV_mapping_functions.R'), local=TRUE)
score_dir <- file.path(cfg$cache_dir, 'reference_scores')

# Section 2. Figure 4c: single-cell internal and cross-reference validation ----
sc_scores <- do.call(rbind, lapply(cfg$sources$SC, function(ref)
  readRDS(file.path(score_dir, paste0(ref, '_scores.rds')))))
sc_boxes <- plot_reference_validation(
  sc_scores, sources=cfg$sources$SC, targets=cfg$sources$SC,
  output=file.path(cfg$figure_dir, 'Figure_4c_CM_SC_AV_reference_validation.pdf')
)
write.csv(sc_boxes, file.path(cfg$table_dir, 'SC_reference_boxplot_statistics.csv'), row.names=FALSE)

# 3. Figure S4d: single-nucleus reference validation ----
sn_scores <- readRDS(file.path(score_dir, 'Fetal_SN_scores.rds'))
sn_boxes <- plot_reference_validation(
  sn_scores, sources=cfg$sources$SN, targets=cfg$sources$SN,
  output=file.path(cfg$figure_dir, 'Figure_S4d_CM_SN_AV_reference_validation.pdf')
)
write.csv(sn_boxes, file.path(cfg$table_dir, 'SN_reference_boxplot_statistics.csv'), row.names=FALSE)
