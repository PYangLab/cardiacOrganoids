
# Accept --name value and --name=value. With no arguments, use defaults.
parse_pipeline_args <- function(defaults, args=commandArgs(trailingOnly=TRUE)) {
  if (!length(args)) return(defaults)
  supplied <- list()
  i <- 1L
  while (i <= length(args)) {
    token <- args[[i]]
    if (!startsWith(token, '--'))
      stop("Expected an argument beginning with '--'; received: ", token)
    token <- substring(token, 3L)
    if (grepl('=', token, fixed=TRUE)) {
      separator <- regexpr('=', token, fixed=TRUE)[1]
      key <- substring(token, 1L, separator-1L)
      value <- substring(token, separator+1L)
      i <- i+1L
    } else {
      key <- token
      if (i == length(args)) stop('Missing value for argument --', key)
      value <- args[[i+1L]]
      i <- i+2L
    }
    key <- gsub('-', '_', key, fixed=TRUE)
    if (!key %in% names(defaults))
      stop('Unknown argument --', gsub('_', '-', key, fixed=TRUE),
           '. Allowed arguments: ',
           paste(paste0('--', gsub('_', '-', names(defaults), fixed=TRUE)), collapse=', '))
    if (key %in% names(supplied)) stop('Argument supplied more than once: --', key)
    template <- defaults[[key]]
    if (is.integer(template)) {
      converted <- suppressWarnings(as.integer(value))
      if (length(converted) != 1L || is.na(converted)) stop('--', key, ' must be an integer')
    } else if (is.numeric(template)) {
      converted <- suppressWarnings(as.numeric(value))
      if (length(converted) != 1L || !is.finite(converted)) stop('--', key, ' must be numeric')
    } else if (is.logical(template)) {
      logical_value <- tolower(value)
      if (logical_value %in% c('true', 't', '1', 'yes')) converted <- TRUE
      else if (logical_value %in% c('false', 'f', '0', 'no')) converted <- FALSE
      else stop('--', key, ' must be TRUE or FALSE')
    } else converted <- as.character(value)
    supplied[[key]] <- converted
  }
  defaults[names(supplied)] <- supplied
  defaults
}

print_pipeline_config <- function(config) {
  message('A/V mapping configuration:')
  for (name in names(config)) message('  ', name, ': ', config[[name]])
}
