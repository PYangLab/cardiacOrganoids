# Running the analyses

[← Back to the project overview](../README.md)

These notes describe the scripts as uploaded. Input files are not bundled with the code. Relative paths in the command-line examples are resolved from the working directories shown; replace the example absolute paths with your own.

## Input checklist

All organoid lists below are named lists of Seurat objects. Gene names, cell names, RNA layers and metadata must match the assumptions in the relevant scripts.

| Module | Expected inputs | Important details |
| :--- | :--- | :--- |
| Identity and maturity | `HCA_v1_sc_processed_filtered.rds`, `HCA_v1_sn_processed_filtered.rds`, `sc_organoids.rds`, `sn_organoids.rds` | Stored under `~/Cardiac_organoids/processed_data/` in the original scripts. Organoids need `celltype`; identity also uses reference `sangerID`, `region`, `donor`, `gender`. |
| Regional mapping (RM) | `cui_processed.rds`, `fetal_hca_sc_processed.rds`, `fetal_hca_sn_processed.rds`, `sc_organoid_list.rds`, `sn_organoid_list.rds`, `Cardiomyocyte_atrial_and_ventricular_module_gene_sets.xlsx` | Filenames are configured in [`config_AV_mapping.R`](../scripts/RM/config_AV_mapping.R). Only cardiomyocytes are scored. |
| Off-target annotation (OTA) | `sc_organoid_list.rds`, `hca_sc_processed.rds`, `fetal_hca_sc_processed.rds`; a Cao fetal Azimuth reference directory | Requires exactly 22 uniquely named organoid datasets. The Azimuth reference directory must contain `ref.Rds` and `idx.annoy`. |

**Filename differences require preparation.** The processing scripts use `sc_organoids.rds` / `sn_organoids.rds`, while RM and OTA expect `sc_organoid_list.rds` / `sn_organoid_list.rds`. OTA expects the adult reference as `hca_sc_processed.rds`. Check that the objects meet each module's requirements before saving under those names; renaming alone does not establish compatibility.

The current upload does not provide a complete data-download manifest or preparation scripts for every external reference and resource. The filenames above identify expected local inputs.

## Dataset processing

Work through the numbered files in [`scripts/dataset_processing/`](../scripts/dataset_processing) in order, preserving the objects needed by later scripts.

1. **Adult heart reference:** read the HDF/H5AD data, align genes and cells, filter by modality, run dimensionality reduction and scReClassify, and harmonise labels. Replace the `~filepath/` output paths.
2. **Organoid datasets:** load the study-specific matrices and metadata, then process the single-cell and single-nucleus datasets. Review each study's input paths and format.
3. **Cell-type annotation:** train modality-specific scClassify models and obtain scClassify, Seurat and SingleR predictions. A consensus requires at least two matching labels. Atrial and ventricular cardiomyocyte annotations are merged to `Cardiomyocyte` for downstream organoid assessment.
4. **Save objects:** split by time point, chamber or condition. The supplied list contains 22 single-cell and 4 single-nucleus subsets. Review the rare-cell-type handling and enable or adapt the commented `saveRDS()` calls.

## Cell identity

From the repository root, load the helper functions:

```r
source("scripts/identity/01_cepo_functions.R")
```

Then work through [`02_cepo_analysis.R`](../scripts/identity/02_cepo_analysis.R) interactively, selecting **either** its single-cell setup block **or** its single-nucleus setup block for each run. Running both setup blocks consecutively replaces the single-cell variables with single-nucleus variables.

The script filters cell types and common genes, calculates and averages reference Cepo scores, and compares them with organoid scores using Pearson correlation. Organoid subsampling uses 10 repetitions with 80% of each cell type retained. The final long-format result is `bootstrap_correlations_df`; save it explicitly if required.

## Maturity

[`01_cytotrace.R`](../scripts/maturity/01_cytotrace.R) runs CytoTRACE2 on combined organoid/reference counts separately for single cells and single nuclei. It creates filtered results in `cytotrace2_filtered_sc` and `cytotrace2_filtered_sn`.

For entropy, load the local helper file before working through the analysis:

```r
source("scripts/maturity/02_entropy_functions.R")
```

[`03_entropy.R`](../scripts/maturity/03_entropy.R) also loads the external file `clean_nodatasets_060720.RData`, which is **not included**. Update that path and check the resource contents; the human entropy calculation requires `G_list_human`, including its `symbol` and `use_for_entropy` fields. Outputs are held in `combined_entropy_df_sc` and `combined_entropy_df_sn`.

CytoTRACE2 scores and expression entropy are complementary transcriptomic measures. Their interpretation as organoid maturity requires biological context.

## Regional mapping

Run from `scripts/RM/`. Use the same output directory across all steps.

The gene workbook must have a `gene_sets` sheet, two header rows, and six gene columns in this order: Cui atrial, Cui ventricular, fetal HCA single-cell atrial, fetal HCA single-cell ventricular, fetal HCA single-nucleus atrial, fetal HCA single-nucleus ventricular.

Review [`config_AV_mapping.R`](../scripts/RM/config_AV_mapping.R) before running:

- Organoid cardiomyocytes are identified by `celltype == "Cardiomyocyte"`.
- The Cui reference uses `celltype == "Cardiomyocyte"` and `Region` values `LA`, `RA`, `LV`, `RV`.
- Fetal references use `mid_grain` labels `AtrialCardiomyocytes` / `VentricularCardiomyocytes`; the single-nucleus reference also uses `region`.
- The fetal HCA reference inputs are expected to contain log1p-CP10K values in the RNA counts layer, with library totals in `n_counts`. The reconstruction helper checks the round trip. Do not substitute a differently encoded matrix without adapting the configuration and preparation.

The following shell commands run the supplied sequence after these inputs are prepared:

```bash
cd scripts/RM

DATA_DIR="/absolute/path/to/prepared_data"
OUT_DIR="/absolute/path/to/results/regional_mapping"
GENE_FILE="$DATA_DIR/Cardiomyocyte_atrial_and_ventricular_module_gene_sets.xlsx"

Rscript 01_load_module_gene_sets.R \
  --data-dir "$DATA_DIR" --out-dir "$OUT_DIR" --gene-file "$GENE_FILE"

for ref in Cui Fetal_SC Fetal_SN; do
  Rscript 02_score_references.R \
    --data-dir "$DATA_DIR" --out-dir "$OUT_DIR" \
    --gene-file "$GENE_FILE" --ref "$ref"
done

Rscript 03_plot_reference_validation.R \
  --data-dir "$DATA_DIR" --out-dir "$OUT_DIR" --gene-file "$GENE_FILE"

for modality in SC SN; do
  Rscript 04_score_organoids.R \
    --data-dir "$DATA_DIR" --out-dir "$OUT_DIR" \
    --gene-file "$GENE_FILE" --modality "$modality"
  Rscript 05_plot_organoid_assignments.R \
    --data-dir "$DATA_DIR" --out-dir "$OUT_DIR" \
    --gene-file "$GENE_FILE" --modality "$modality" --t 0.10
done
```

The assignment threshold defaults to `t = 0.10`. A cell is assigned atrial when `A_score - V_score > t` and `A_score > 0`; the ventricular rule is symmetric. Otherwise it is unresolved. Single-cell organoids are scored with Cui and fetal HCA single-cell modules; single-nucleus organoids use fetal HCA single-nucleus modules.

Results are written under `intermediate/`, `tables/`, `figures/` and `logs/` within the chosen output directory, including reference validation plots, per-cell scores, assignment/proportion tables and organoid assignment plots.

## Off-target annotation

Run from `scripts/OTA/`, keeping the same output directory across all steps.

Before running, check the input metadata:

- Organoids: `celltype` plus a usable normalised RNA `data` layer.
- Adult reference: `scReClassify_label`, `age`, `donor`.
- Fetal reference: `coarse_grain`, `mid_grain`, `fine_grain`, `age`, `week`, `donor`, `cell_or_nuclei`. Every entry must be marked `cell`; the reference must also contain at least 30 PCA dimensions for label transfer.

```bash
cd scripts/OTA

DATA_DIR="/absolute/path/to/prepared_data"
OUT_DIR="/absolute/path/to/results/offtarget"
CAO_REFERENCE_DIR="/absolute/path/to/cao_fetal_azimuth"

Rscript 01_validate_inputs.R \
  --data-dir "$DATA_DIR" --out-dir "$OUT_DIR"

Rscript 02_classify_adult_unassigned_with_fetal_hca.R \
  --data-dir "$DATA_DIR" --out-dir "$OUT_DIR" --workers 4

Rscript 03_map_doubly_unassigned_with_cao_fetal_azimuth.R \
  --out-dir "$OUT_DIR" --cao-reference-dir "$CAO_REFERENCE_DIR"
```

Step 1 validates the prepared files and saves per-dataset queries. Step 2 classifies only organoid cells labelled `unassigned`, using fetal-heart scClassify, Seurat and SingleR predictions. Step 3 maps cells still unassigned to the Cao fetal multi-organ Azimuth reference.

The stages write input summaries, prediction tables, checkpoints and annotated organoid objects. The final list is saved at:

```text
<OUT_DIR>/03_cao_fetal_azimuth/sc_organoid_list_with_fetal_hca_and_cao_labels.rds
```

For visualisation, open [`04_extra-cardiac_visualization.Rmd`](../scripts/OTA/04_extra-cardiac_visualization.Rmd), set `results_dir` to the same output directory, and review the chunks before rendering with R Markdown.

Use a new output directory when inputs or analysis settings change: the scripts reuse existing query files, predictions and checkpoints.
