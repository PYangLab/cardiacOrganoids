# Running the analyses

[Back to the project overview](../README.md)

Run each command-line example from the directory shown and replace the example absolute paths with your own. Relative paths are resolved from that working directory. You'll need to obtain the input files separately.

## Input checklist

Each organoid list must be a named list of Seurat objects. Check that its gene names, cell names, RNA layers and metadata match what the scripts expect.

| Module | Expected inputs | Important details |
| :--- | :--- | :--- |
| Identity and maturity | `HCA_v1_sc_processed_filtered.rds`, `HCA_v1_sn_processed_filtered.rds`, `sc_organoids.rds`, `sn_organoids.rds` | The scripts read these from `~/Cardiac_organoids/processed_data/`. Organoids need `celltype`. Identity analysis also uses reference metadata `sangerID`, `region`, `donor`, `gender`. |
| Regional mapping (RM) | `cui_processed.rds`, `fetal_hca_sc_processed.rds`, `fetal_hca_sn_processed.rds`, `sc_organoid_list.rds`, `sn_organoid_list.rds`, `Cardiomyocyte_atrial_and_ventricular_module_gene_sets.xlsx` | Set filenames in [`config_AV_mapping.R`](../scripts/RM/config_AV_mapping.R). Scoring uses cardiomyocytes only. |
| Off-target annotation (OTA) | `sc_organoid_list.rds`, `hca_sc_processed.rds`, `fetal_hca_sc_processed.rds` and a Cao fetal Azimuth reference directory | The organoid list must have exactly 22 uniquely named datasets. The Azimuth reference directory needs `ref.Rds` and `idx.annoy`. |

The processing scripts use `sc_organoids.rds` / `sn_organoids.rds`, while RM and OTA expect `sc_organoid_list.rds` / `sn_organoid_list.rds`. OTA reads the adult reference from `hca_sc_processed.rds`. Before saving objects under these names, check their contents against the requirements above. A matching filename alone does not make an object compatible.

The repository does not yet include a complete list of data downloads or preparation scripts for every external reference and resource. The filenames above are the local files the scripts expect.

## Dataset processing

Work through the numbered scripts in [`scripts/dataset_processing/`](../scripts/dataset_processing) in order. Keep the R objects needed by later steps in your session.

1. Read the adult heart HDF/H5AD data, align genes and cells, and filter each modality. Run dimensionality reduction and scReClassify, then harmonise the labels. Change the `~filepath/` output paths.
2. Load the organoid matrices and metadata for each study, checking its input paths and file format. Process single-cell and single-nucleus data separately.
3. Train scClassify models for each modality and obtain predictions from scClassify, Seurat and SingleR. At least two labels must match for a consensus. Merge atrial and ventricular cardiomyocyte labels to `Cardiomyocyte` for the downstream organoid analyses.
4. Split the organoids by time point, chamber or condition to create 22 single-cell and 4 single-nucleus subsets. Check how rare cell types are handled, then enable or adapt the commented `saveRDS()` calls.

## Cell identity

From the repository root, load the helper functions:

```r
source("scripts/identity/01_cepo_functions.R")
```

Work through [`02_cepo_analysis.R`](../scripts/identity/02_cepo_analysis.R) interactively for one modality at a time. Choose the single-cell setup block for a single-cell run, or the single-nucleus block for a single-nucleus run. Running both consecutively overwrites the single-cell variables.

The script filters cell types and finds common genes before calculating and averaging the reference Cepo scores. It then compares organoid and reference scores with Pearson correlation. For each organoid, it repeats subsampling 10 times, keeping 80% of each cell type in each subsample. Save `bootstrap_correlations_df` if you need the final long-format results on disk.

## Maturity

[`01_cytotrace.R`](../scripts/maturity/01_cytotrace.R) combines organoid and reference counts, then runs CytoTRACE2 separately for single cells and single nuclei. The filtered results are `cytotrace2_filtered_sc` and `cytotrace2_filtered_sn`.

Load the entropy helper functions before running the entropy analysis:

```r
source("scripts/maturity/02_entropy_functions.R")
```

[`03_entropy.R`](../scripts/maturity/03_entropy.R) also needs `clean_nodatasets_060720.RData`, an external file that is not included here. Update its path and check that it provides `G_list_human` with the `symbol` and `use_for_entropy` fields needed for human entropy calculations. The results are `combined_entropy_df_sc` and `combined_entropy_df_sn`.

CytoTRACE2 scores and expression entropy measure different aspects of the transcriptome. Interpret both alongside biological evidence when assessing organoid maturity.

## Regional mapping

Run these scripts from `scripts/RM/` and keep the same output directory throughout.

The gene workbook must have a `gene_sets` sheet, two header rows, and six gene columns in this order: Cui atrial, Cui ventricular, fetal HCA single-cell atrial, fetal HCA single-cell ventricular, fetal HCA single-nucleus atrial, fetal HCA single-nucleus ventricular.

Check the input settings in [`config_AV_mapping.R`](../scripts/RM/config_AV_mapping.R):

- Organoid cardiomyocytes must have `celltype == "Cardiomyocyte"`.
- The Cui reference uses `celltype == "Cardiomyocyte"` and `Region` values `LA`, `RA`, `LV`, `RV`.
- Fetal references use `mid_grain` labels `AtrialCardiomyocytes` / `VentricularCardiomyocytes`. The single-nucleus reference also needs `region`.
- The fetal HCA references need log1p-CP10K values in the RNA counts layer and library totals in `n_counts`. The helper reconstructs counts and checks the round trip. If your matrix uses a different encoding, adapt the configuration and input preparation first.

Once the inputs are ready, run:

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

The default assignment threshold is `t = 0.10`. The script assigns an atrial label when `A_score - V_score > t` and `A_score > 0`. The ventricular rule is symmetric, and cells meeting neither rule are unresolved. Single-cell organoids use Cui and fetal HCA single-cell modules. Single-nucleus organoids use fetal HCA single-nucleus modules.

The scripts save results under `intermediate/`, `tables/`, `figures/` and `logs/` in your output directory. These include reference validation plots, per-cell scores, assignment and proportion tables, and organoid assignment plots.

## Off-target annotation

Run these scripts from `scripts/OTA/` with the same output directory for every step.

Before running, check the input metadata:

- Organoids: `celltype` plus a usable normalised RNA `data` layer.
- Adult reference: `scReClassify_label`, `age`, `donor`.
- Fetal reference: `coarse_grain`, `mid_grain`, `fine_grain`, `age`, `week`, `donor`, `cell_or_nuclei`. Every entry must be marked `cell`. Label transfer also needs at least 30 PCA dimensions in the reference.

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

The first script checks the input files and saves queries for each dataset. The second takes organoid cells labelled `unassigned` and classifies them against the fetal heart reference with scClassify, Seurat and SingleR. Cells still unassigned then pass to the third script for mapping to the Cao fetal multi-organ Azimuth reference.

The scripts save input summaries, prediction tables, checkpoints and annotated organoid objects. The final organoid list is at:

```text
<OUT_DIR>/03_cao_fetal_azimuth/sc_organoid_list_with_fetal_hca_and_cao_labels.rds
```

To plot the results, open [`04_extra-cardiac_visualization.Rmd`](../scripts/OTA/04_extra-cardiac_visualization.Rmd) and set `results_dir` to your output directory. Check the chunks before rendering with R Markdown.

Choose a new output directory if you change the inputs or analysis settings. Otherwise, the scripts reuse existing query files, predictions and checkpoints.
