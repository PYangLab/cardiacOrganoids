get_rna_layer <- function(object, layer = "data") {
  if (!"RNA" %in% names(object@assays)) {
    stop("Object has no RNA assay")
  }
  
  # Seurat v5 path.
  if (exists("LayerData", mode = "function")) {
    value <- tryCatch(
      LayerData(object = object, assay = "RNA", layer = layer),
      error = function(e) NULL
    )
    if (!is.null(value)) {
      return(value)
    }
  }
  
  # Legacy Seurat assay path.
  assay <- object@assays$RNA
  value <- tryCatch(
    GetAssayData(object = object, assay = "RNA", slot = layer),
    error = function(e) NULL
  )
  if (!is.null(value) && length(dim(value)) == 2L && all(dim(value) > 0L)) {
    return(value)
  }
  
  # Read-only fallback for inspecting Seurat v5 Assay5 objects in an older
  # Seurat installation. Classification itself should be run with Seurat v5.
  layers <- attr(assay, "layers")
  if (is.null(layers) || !layer %in% names(layers)) {
    stop("RNA layer is unavailable: ", layer)
  }
  value <- layers[[layer]]
  feature_map <- attr(assay, "features")
  cell_map <- attr(assay, "cells")
  feature_names <- rownames(feature_map)
  cell_names <- rownames(cell_map)
  
  if (nrow(value) != length(feature_names) || ncol(value) != length(cell_names)) {
    stop(
      "The fallback Assay5 reader only supports layers containing all features ",
      "and cells; use Seurat v5 for this object"
    )
  }
  dimnames(value) <- list(feature_names, cell_names)
  value
}

get_object_cells <- function(object) {
  cells <- tryCatch(Cells(object), error = function(e) NULL)
  if (is.null(cells) || length(cells) == 0L) {
    cells <- rownames(object@meta.data)
  }
  cells
}


write_csv_gz <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.csv(x, con, row.names = FALSE, quote = TRUE)
  invisible(path)
}

exact_consensus <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  votes <- table(x)
  if (length(votes) > 0L && max(votes) >= 2L) {
    return(names(votes)[which.max(votes)])
  }
  "unassigned"
}


standardize_adult_hca_labels <- function(x) {
  x <- as.character(x)
  replacements <- c(
    "Mast cell" = "Mast.cell",
    "Neural cell" = "Neural.cell",
    "Mural cell" = "Mural.cell",
    "Endothelial cell" = "Endothelial.cell",
    "Ventricular Cardiomyocyte" = "Ventricular.cardiomyocyte",
    "Atrial Cardiomyocyte" = "Atrial.cardiomyocyte",
    "Lymphatic Endothelial cell" = "Lymphatic.endothelial.cell",
    "Mesothelial cell" = "Mesothelial.cell"
  )
  for (from in names(replacements)) {
    x[grepl(from, x, fixed = TRUE)] <- replacements[[from]]
  }
  x
}


assert_unique_keys <- function(x, keys, label) {
  duplicated_key <- duplicated(x[, keys, drop = FALSE])
  if (any(duplicated_key)) {
    stop(label, " contains duplicated keys: ", paste(keys, collapse = " + "))
  }
  invisible(TRUE)
}



parse_pipeline_args <- function(
    defaults,
    args = commandArgs(trailingOnly = TRUE)
) {
  if (!is.list(defaults) ||
      is.null(names(defaults)) ||
      any(!nzchar(names(defaults)))) {
    stop("defaults must be a named list")
  }
  
  if (!length(args)) {
    return(defaults)
  }
  
  supplied <- list()
  i <- 1L
  
  while (i <= length(args)) {
    token <- args[[i]]
    
    if (!startsWith(token, "--")) {
      stop("Expected an argument beginning with '--'; received: ", token)
    }
    
    token <- substring(token, 3L)
    
    if (grepl("=", token, fixed = TRUE)) {
      separator <- regexpr("=", token, fixed = TRUE)[1]
      key <- substring(token, 1L, separator - 1L)
      value <- substring(token, separator + 1L)
      i <- i + 1L
    } else {
      key <- token
      
      if (i == length(args)) {
        stop("Missing value for argument --", key)
      }
      
      value <- args[[i + 1L]]
      i <- i + 2L
    }
    
    # Translate command-line names such as --data-dir into data_dir.
    key <- gsub("-", "_", key, fixed = TRUE)
    
    if (!key %in% names(defaults)) {
      stop(
        "Unknown argument --",
        gsub("_", "-", key, fixed = TRUE),
        ". Allowed arguments are: ",
        paste(
          paste0("--", gsub("_", "-", names(defaults), fixed = TRUE)),
          collapse = ", "
        )
      )
    }
    
    if (key %in% names(supplied)) {
      stop("Argument supplied more than once: --", gsub("_", "-", key))
    }
    
    template <- defaults[[key]]
    
    if (is.integer(template)) {
      converted <- suppressWarnings(as.integer(value))
      if (length(converted) != 1L || is.na(converted)) {
        stop("--", gsub("_", "-", key), " must be an integer")
      }
    } else if (is.numeric(template)) {
      converted <- suppressWarnings(as.numeric(value))
      if (length(converted) != 1L || !is.finite(converted)) {
        stop("--", gsub("_", "-", key), " must be numeric")
      }
    } else if (is.logical(template)) {
      logical_value <- tolower(value)
      
      if (logical_value %in% c("true", "t", "1", "yes")) {
        converted <- TRUE
      } else if (logical_value %in% c("false", "f", "0", "no")) {
        converted <- FALSE
      } else {
        stop("--", gsub("_", "-", key), " must be TRUE or FALSE")
      }
    } else {
      converted <- as.character(value)
    }
    
    supplied[[key]] <- converted
  }
  
  defaults[names(supplied)] <- supplied
  defaults
}
