# Shared input preparation, Seurat scoring and plotting functions.

suppressPackageStartupMessages({library(Seurat); library(Matrix); library(ggplot2)})
invisible(Sys.setlocale('LC_COLLATE', 'C'))
RNGkind('Mersenne-Twister', 'Inversion', 'Rejection')
future::plan('sequential')
for (path in c(cfg$cache_dir, cfg$table_dir, cfg$figure_dir, cfg$log_dir))
  dir.create(path, recursive=TRUE, showWarnings=FALSE)

read_module_gene_sets <- function() {
  # Row 1 contains reference names, row 2 contains module names and counts.
  # Columns A:B = Cui, C:D = fetal HCA sc, E:F = fetal HCA sn.
  table <- readxl::read_excel(cfg$gene_file, sheet=cfg$gene_sheet, skip=2,
                            col_names=FALSE, col_types=rep('text', 6),
                            .name_repair='minimal')
  columns <- lapply(table, function(x) {
    x <- trimws(x[!is.na(x)])

    sort(unique(x[nzchar(x)]))
  })
  list(Cui=list(A=columns[[1]], V=columns[[2]]),
       Fetal_SC=list(A=columns[[3]], V=columns[[4]]),
       Fetal_SN=list(A=columns[[5]], V=columns[[6]]))
}

rna_data <- function(obj, layer='data') {
  if (utils::packageVersion('SeuratObject') >= '5.0.0')
    SeuratObject::GetAssayData(obj, assay='RNA', layer=layer)
  else
    SeuratObject::GetAssayData(obj, assay='RNA', slot=layer)
}

rna_counts <- function(obj) {
  if (inherits(obj[['RNA']], 'Assay5')) {
    layers <- grep('^counts($|[.])', SeuratObject::Layers(obj[['RNA']], search=NA), value=TRUE)
    if (length(layers) > 1L)
      obj <- SeuratObject::JoinLayers(obj, assay='RNA', layers=layers, new='counts')
  }
  as(rna_data(obj, 'counts'), 'dgCMatrix')
}

recover_hca_counts <- function(x, library_total, scale=10000) {

  stored_values <- x@x
  recovered_values <- numeric(length(stored_values))
  pointers <- x@p
  for (j in seq_len(ncol(x))) {
    lo <- pointers[j] + 1L
    hi <- pointers[j + 1L]
    if (lo > hi) stop('Empty HCA library: ', colnames(x)[j])
    stored <- stored_values[lo:hi]
    inverse <- expm1(stored) * library_total[j] / scale
    error <- max(abs(inverse - round(inverse)))
    if (!is.finite(error) || error >= 0.02)
      stop('HCA count reconstruction failed for ', colnames(x)[j])
    counts <- round(inverse)
    log_error <- max(abs(log1p(counts * scale / library_total[j]) - stored))
    if (sum(counts) != library_total[j] || !is.finite(log_error) || log_error >= 1e-6)
      stop('HCA count reconstruction round trip failed for ', colnames(x)[j])
    recovered_values[lo:hi] <- counts
    if (j %% 1000L == 0L || j == ncol(x))
      message('Recovered ', j, '/', ncol(x), ' HCA cells')
  }
  x@x <- recovered_values
  x
}

normalize_cm_counts <- function(counts, metadata) {
  # Retain every RNA feature for library-size normalization and the control pool.
  obj <- Seurat::CreateSeuratObject(counts=counts, meta.data=metadata, assay='RNA',
                                    min.cells=0, min.features=0)
  obj <- Seurat::NormalizeData(obj, assay='RNA', normalization.method='LogNormalize',
                                scale.factor=cfg$normalization_scale, verbose=FALSE)
  assay <- SeuratObject::CreateAssayObject(data=rna_data(obj))
  SeuratObject::Key(assay) <- 'rna_'
  SeuratObject::CreateSeuratObject(counts=assay, assay='RNA', meta.data=metadata)
}

score_av_modules <- function(obj, modules, sources) {
  parts <- lapply(sources, function(source) {
    missing <- setdiff(unlist(modules[[source]], use.names=FALSE), rownames(obj))
    if (length(missing)) stop(source, ': supplied module genes absent from RNA: ',
                              paste(missing, collapse=', '))
    scored <- Seurat::AddModuleScore(
      object=obj, features=modules[[source]][c('A', 'V')],
      assay='RNA', name='AV_'
    )
    data.frame(cell=colnames(scored), source=source,
               A_score=scored[[]]$AV_1, V_score=scored[[]]$AV_2,
               stringsAsFactors=FALSE)
  })
  do.call(rbind, parts)
}

box_statistics <- function(values) {
  quartiles <- unname(quantile(values, c(0.25, 0.50, 0.75), type=7))
  iqr <- quartiles[3] - quartiles[1]
  inside <- values[values >= quartiles[1] - 1.5*iqr & values <= quartiles[3] + 1.5*iqr]
  data.frame(n=length(values), q1=quartiles[1], median=quartiles[2], q3=quartiles[3],
             whisker_low=min(inside), whisker_high=max(inside))
}

reference_box_statistics <- function(scores) {
  keys <- unique(scores[, c('source', 'target', 'label')])
  parts <- lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, , drop=FALSE]
    z <- scores[scores$source == key$source & scores$target == key$target &
                  scores$label == key$label, , drop=FALSE]
    do.call(rbind, lapply(c('Atrial', 'Ventricular'), function(module) {
      values <- z[[if (module == 'Atrial') 'A_score' else 'V_score']]
      cbind(key, module=module, box_statistics(values))
    }))
  })
  do.call(rbind, parts)
}

plot_reference_validation <- function(scores, sources, targets, output) {
  boxes <- reference_box_statistics(scores)
  rows <- as.vector(t(outer(targets, c('A', 'V'), paste, sep='__')))
  boxes$row <- paste(boxes$target, boxes$label, sep='__')
  boxes$y <- match(boxes$row, rev(rows)) + ifelse(boxes$module == 'Atrial', 0.16, -0.16)
  boxes$source <- factor(boxes$source, levels=sources)
  bands <- data.frame(y=seq_along(rows),
                      label=sub('.*__', '', rev(rows)))
  axis_labels <- vapply(rev(rows), function(row) {
    key <- strsplit(row, '__', fixed=TRUE)[[1]]
    paste0(cfg$source_names[key[1]], '\n', if (key[2] == 'A') 'atrial cells' else 'ventricular cells')
  }, character(1))
  strips <- setNames(paste0(cfg$source_names[sources], '\ndefined cardiomyocyte\ngene modules'), sources)
  limits <- range(c(boxes$whisker_low, boxes$whisker_high, 0))
  padding <- max(0.025, diff(limits)*0.05)
  p <- ggplot(boxes) +
    geom_rect(data=bands[bands$label == 'A', , drop=FALSE],
              aes(xmin=-Inf, xmax=Inf, ymin=y-0.5, ymax=y+0.5),
              fill='#DEEBF5', inherit.aes=FALSE) +
    geom_rect(data=bands[bands$label == 'V', , drop=FALSE],
              aes(xmin=-Inf, xmax=Inf, ymin=y-0.5, ymax=y+0.5),
              fill='#FBE5D6', inherit.aes=FALSE) +
    geom_vline(xintercept=0, linetype=2, linewidth=0.3, color='grey55') +
    geom_segment(aes(x=whisker_low, xend=whisker_high, y=y, yend=y),
                 linewidth=0.35, color='grey45') +
    geom_rect(aes(xmin=q1, xmax=q3, ymin=y-0.115, ymax=y+0.115, fill=module),
              color='grey30', linewidth=0.3) +
    geom_segment(aes(x=median, xend=median, y=y-0.115, yend=y+0.115), linewidth=0.45) +
    facet_grid(~source, labeller=labeller(source=as_labeller(strips))) +
    scale_fill_manual(values=cfg$colors[c('Atrial', 'Ventricular')],
                      breaks=c('Atrial', 'Ventricular')) +
    scale_y_continuous(breaks=seq_along(rows), labels=axis_labels,
                       limits=c(0.5, length(rows)+0.5), expand=expansion(mult=0)) +
    coord_cartesian(xlim=limits+c(-padding, padding)) +
    labs(x='Gene module score', y=NULL, fill=NULL) +
    theme_bw(base_size=10) +
    theme(panel.grid=element_blank(), panel.spacing=grid::unit(0.45, 'cm'),
          strip.background=element_rect(fill='grey85', color='grey40'),
          strip.text=element_text(size=9), legend.position='bottom',
          legend.key.size=grid::unit(0.3, 'cm'), axis.text.y=element_text(color='black'),
          plot.margin=margin(6, 16, 4, 6))
  ggsave(output, p, width=if (length(sources) == 2L) 6.6 else 4.2,
         height=if (length(targets) == 2L) 3.5 else 2.65, useDingbats=FALSE)
  invisible(boxes)
}

assign_av <- function(A_score, V_score, t) {
  label <- rep('Unresolved', length(A_score))
  label[A_score - V_score > t & A_score > 0] <- 'Atrial'
  label[V_score - A_score > t & V_score > 0] <- 'Ventricular'
  factor(label, levels=c('Atrial', 'Ventricular', 'Unresolved'))
}

organoid_proportions <- function(scores, sources) {
  parts <- lapply(sort(unique(scores$dataset)), function(dataset) {
    do.call(rbind, lapply(sources, function(source) {
      z <- scores[scores$dataset == dataset & scores$source == source, , drop=FALSE]
      counts <- table(z$assignment)
      data.frame(dataset=dataset, source=source, assignment=names(counts),
                 n_cells=as.integer(counts), total_cells=nrow(z),
                 fraction=as.numeric(counts)/nrow(z), threshold=cfg$t)
    }))
  })
  proportions <- do.call(rbind, parts)
  if (length(sources) == 2L) {
    combined <- aggregate(fraction ~ dataset + assignment + threshold,
                          data=proportions, FUN=mean)
    combined$source <- 'Combined'
    combined$n_cells <- NA_integer_
    sizes <- proportions[!duplicated(proportions$dataset), c('dataset', 'total_cells')]
    combined$total_cells <- sizes$total_cells[match(combined$dataset, sizes$dataset)]
    proportions <- rbind(proportions, combined[, names(proportions)])
  }

  proportions
}

plot_organoid_proportions <- function(proportions, sources, output) {
  datasets <- sort(unique(proportions$dataset))
  columns <- c(sources, if (length(sources) == 2L) 'Combined')
  q <- proportions
  q$source <- factor(q$source, levels=columns)
  q$assignment <- factor(q$assignment, levels=c('Atrial', 'Ventricular', 'Unresolved'))
  q <- q[order(q$source, match(q$dataset, datasets), q$assignment), ]
  key <- interaction(q$source, q$dataset, drop=TRUE)
  q$upper <- ave(q$fraction, key, FUN=cumsum)
  q$lower <- q$upper - q$fraction
  q$y <- match(q$dataset, rev(datasets))
  q$text <- ifelse(q$fraction >= cfg$label_min_fraction,
                   sprintf('%.0f%%', 100*q$fraction), '')
  strips <- c(cfg$source_names[sources], Combined='Combined')
  p <- ggplot(q, aes(fill=assignment)) +
    geom_rect(aes(xmin=lower, xmax=upper, ymin=y-0.4, ymax=y+0.4)) +
    geom_text(aes(x=(lower+upper)/2, y=y, label=text), size=2.7, color='white') +
    facet_grid(~source, labeller=labeller(source=as_labeller(strips))) +
    scale_fill_manual(values=cfg$colors, breaks=c('Atrial', 'Ventricular', 'Unresolved'), drop=FALSE) +
    scale_x_continuous(breaks=seq(0, 1, 0.25), labels=c('0%', '25%', '50%', '75%', '100%'),
                       expand=expansion(mult=0)) +
    scale_y_continuous(breaks=seq_along(datasets), labels=rev(datasets),
                       limits=c(0.5, length(datasets)+0.5), expand=expansion(mult=0)) +
    coord_cartesian(xlim=c(0, 1)) +
    labs(x='Proportion of cardiomyocytes', y=NULL, fill=NULL) +
    theme_bw(base_size=10) +
    theme(panel.grid=element_blank(), panel.spacing=grid::unit(0.8, 'cm'),
          strip.background=element_rect(fill='grey85', color='grey40'),
          strip.text=element_text(size=10), legend.position='bottom',
          legend.key.size=grid::unit(0.3, 'cm'), axis.text.y=element_text(color='black'),
          plot.margin=margin(6, 18, 4, 6))
  ggsave(output, p, width=if (length(sources) == 2L) 10 else 4.8,
         height=if (length(sources) == 2L) 8 else 2.8, useDingbats=FALSE)
  invisible(p)
}
