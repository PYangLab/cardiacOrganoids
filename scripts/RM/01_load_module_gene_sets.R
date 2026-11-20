#!/usr/bin/env Rscript

# Run all sections once. 
# The Excel workbook is the module-definition input.

# Section 1. Setup ----
code_dir <- '.'
source(file.path(code_dir, 'AV_mapping_args.R'), local=TRUE)
# Edit defaults for an interactive run; command-line arguments override them.
config <- parse_pipeline_args(defaults=list(
  data_dir='../data',
  out_dir='./results',
  gene_file='../data/Cardiomyocyte_atrial_and_ventricular_module_gene_sets.xlsx'
))
print_pipeline_config(config)
source(file.path(code_dir, 'config_AV_mapping.R'), local=TRUE)
source(file.path(code_dir, 'AV_mapping_functions.R'), local=TRUE)

# Section 2. Read the final modules from Excel ----
modules <- read_module_gene_sets()
module_sizes <- do.call(rbind, lapply(names(modules), function(reference)
  data.frame(reference=reference, module=c('Atrial', 'Ventricular'),
             genes=lengths(modules[[reference]]))))
print(module_sizes, row.names=FALSE)

# Section 3. Save the fixed lists for all later scripts ----
saveRDS(modules, file.path(cfg$cache_dir, 'module_gene_sets.rds'))
write.csv(module_sizes, file.path(cfg$table_dir, 'module_gene_counts.csv'), row.names=FALSE)
saveRDS(list(modules=modules, normalization_scale=cfg$normalization_scale,
             scoring=cfg$scoring), file.path(cfg$cache_dir, 'module_scoring_settings.rds'))
