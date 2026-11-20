# Cardiomyocyte A/V mapping from the supplied, final module gene sets.
# Each entry script supplies config through parse_pipeline_args().
# config contains user arguments. cfg contains the derived paths and analysis
# settings. config$data_dir is input-only; cfg$cache_dir holds generated files.
INPUT_DIR <- path.expand(config$data_dir)
OUTPUT_DIR <- path.expand(config$out_dir)

cfg <- list(
  gene_file = path.expand(config$gene_file),
  gene_sheet = 'gene_sets',

  cache_dir = file.path(OUTPUT_DIR, 'intermediate'),
  table_dir = file.path(OUTPUT_DIR, 'tables'),
  figure_dir = file.path(OUTPUT_DIR, 'figures'),
  log_dir = file.path(OUTPUT_DIR, 'logs'),
  normalization_scale = 10000,
  scoring = list(method='Seurat::AddModuleScore', parameters='Seurat package defaults'),
  # The assignment threshold is supplied by script 05 (--t; default 0.10).
  colors = c(Atrial='#327CB0', Ventricular='#C96727', Unresolved='#999999'),
  label_min_fraction = 0.05,
  source_names = c(Cui='Cui et al.', Fetal_SC='Fetal HCA', Fetal_SN='Fetal HCA'),
  sources = list(SC=c('Cui', 'Fetal_SC'), SN='Fetal_SN'),
  references = list(
    Cui = list(
      object_path=file.path(INPUT_DIR, 'cui_processed.rds'), modality='SC',
      celltype_column='celltype', cm_labels='Cardiomyocyte',
      label_column='Region', zone_map=c(LA='A', RA='A', LV='V', RV='V'),
      counts_mode='integer_counts'
    ),
    Fetal_SC = list(
      object_path=file.path(INPUT_DIR, 'fetal_hca_sc_processed.rds'), modality='SC',
      celltype_column='mid_grain',
      cm_labels=c('AtrialCardiomyocytes', 'VentricularCardiomyocytes'),
      label_column='mid_grain',
      zone_map=c(AtrialCardiomyocytes='A', VentricularCardiomyocytes='V'),
      counts_mode='recover_log1p_cp10k', count_library_column='n_counts',
      recovery_scale=10000
    ),
    Fetal_SN = list(
      object_path=file.path(INPUT_DIR, 'fetal_hca_sn_processed.rds'), modality='SN',
      celltype_column='mid_grain',
      cm_labels=c('AtrialCardiomyocytes', 'VentricularCardiomyocytes'),
      label_column='region',
      zone_map=c(atria='A', 'left ventricle'='V', 'right ventricle'='V'),
      counts_mode='recover_log1p_cp10k', count_library_column='n_counts',
      recovery_scale=10000
    )
  ),
  organoid_files = c(
    SC=file.path(INPUT_DIR, 'sc_organoid_list.rds'),
    SN=file.path(INPUT_DIR, 'sn_organoid_list.rds')
  ),
  organoid_cm_column = 'celltype',
  organoid_cm_labels = 'Cardiomyocyte'
)
