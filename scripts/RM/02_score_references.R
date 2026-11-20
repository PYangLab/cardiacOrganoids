#!/usr/bin/env Rscript


# Run section 1 once; repeat sections 2-4 for Cui, Fetal_SC and Fetal_SN.

# Section 1. Setup ----
# Reads the original reference objects. No DE output or prepared Step 1 file is needed.
code_dir <- '.'
source(file.path(code_dir, 'AV_mapping_args.R'), local=TRUE)
# Edit defaults for an interactive run; e.g. Rscript 02_score_references.R --ref Fetal_SC.
config <- parse_pipeline_args(defaults=list(
  data_dir='../data',
  out_dir='./results',
  gene_file='../data/Cardiomyocyte_atrial_and_ventricular_module_gene_sets.xlsx',
  ref='Cui'
))
print_pipeline_config(config)
source(file.path(code_dir, 'config_AV_mapping.R'), local=TRUE)
source(file.path(code_dir, 'AV_mapping_functions.R'), local=TRUE)
modules <- readRDS(file.path(cfg$cache_dir, 'module_gene_sets.rds'))
score_dir <- file.path(cfg$cache_dir, 'reference_scores')
dir.create(score_dir, recursive=TRUE, showWarnings=FALSE)

# Section 2. Choose one reference and select its labelled cardiomyocytes ----
# Interactive: set config$ref <- 'Fetal_SC' or 'Fetal_SN' before rerunning sections 2-4.
ref <- config$ref 
rc <- cfg$references[[ref]]
obj <- readRDS(rc$object_path)
md <- obj[[]]
if ('cell' %in% names(md)) rownames(md) <- as.character(md$cell)
label <- as.character(md[[rc$label_column]])
keep <- as.character(md[[rc$celltype_column]]) %in% rc$cm_labels &
  label %in% names(rc$zone_map)
md <- md[keep, , drop=FALSE]
md$AV <- unname(rc$zone_map[label[keep]])
table(md$AV)

# Section 3. Prepare and score the full normalized CM RNA matrix ----
counts <- rna_counts(obj)[, rownames(md), drop=FALSE]
rm(obj)
gc()
if (rc$counts_mode == 'recover_log1p_cp10k')
  counts <- recover_hca_counts(counts, as.numeric(as.character(md[[rc$count_library_column]])),
                              rc$recovery_scale)
seu <- normalize_cm_counts(counts, md)
rm(counts)
gc()

scores <- score_av_modules(seu, modules, cfg$sources[[rc$modality]])
scores$target <- ref
scores$modality <- rc$modality
scores$label <- md$AV[match(scores$cell, rownames(md))]

# Section 4. Save per-cell scores ----
saveRDS(scores, file.path(score_dir, paste0(ref, '_scores.rds')))
write.csv(scores, file.path(cfg$table_dir, paste0(ref, '_reference_scores.csv')), row.names=FALSE)
capture.output(sessionInfo(), file=file.path(cfg$log_dir, paste0(ref, '_reference_sessionInfo.txt')))

rm(seu)
gc()

# After all three references finish, run 03_plot_reference_validation.R.
