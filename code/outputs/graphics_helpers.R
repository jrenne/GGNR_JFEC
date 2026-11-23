# Portable PDF output with embedded fonts. Prefer a working Cairo device;
# otherwise use base PDF followed by Poppler's pdftocairo.
if (!exists("ggnr_graphics", inherits = FALSE)) {
  ggnr_graphics <- new.env(parent = emptyenv())
  ggnr_graphics$devices <- list()
  ggnr_graphics$cairo <- NULL
}

ggnr_pdf_converter <- function() {
  converter <- Sys.getenv("GGNR_PDF_CONVERTER", "")
  if (!nzchar(converter)) converter <- unname(Sys.which("pdftocairo"))
  if (nzchar(converter) && !file.exists(converter)) {
    stop("GGNR_PDF_CONVERTER does not name an existing executable.")
  }
  converter
}

ggnr_has_cairo <- function() {
  if (!is.null(ggnr_graphics$cairo)) return(ggnr_graphics$cairo)
  test_file <- tempfile(fileext = ".pdf")
  old_device <- dev.cur()
  ok <- suppressWarnings(tryCatch({
    grDevices::cairo_pdf(test_file, width = 2, height = 2)
    if (dev.cur() == old_device) stop("Cairo did not open a device.")
    plot.new(); text(0.5, 0.5, "Font test")
    grDevices::dev.off()
    file.exists(test_file) && file.info(test_file)$size > 100
  }, error = function(e) FALSE))
  if (dev.cur() != old_device) grDevices::dev.off()
  unlink(test_file)
  ggnr_graphics$cairo <- ok
  ok
}

ggnr_pdf <- function(file, width = 7, height = 7, ...) {
  args <- list(...)
  cairo <- ggnr_has_cairo()
  converter <- if (cairo) "" else ggnr_pdf_converter()
  if (!cairo && !nzchar(converter)) {
    stop(paste("PDF figures need working Cairo support or pdftocairo (Poppler).",
               "Install Poppler or set GGNR_PDF_CONVERTER to its executable."))
  }
  if (cairo) {
    args$useDingbats <- args$paper <- NULL
    do.call(grDevices::cairo_pdf,
            c(list(filename = file, width = width, height = height), args))
  } else {
    args$useDingbats <- FALSE
    do.call(grDevices::pdf,
            c(list(file = file, width = width, height = height), args))
  }
  ggnr_graphics$devices[[as.character(dev.cur())]] <-
    list(file = file, converter = converter)
  invisible(dev.cur())
}

ggnr_dev_off <- function() {
  device <- as.character(dev.cur())
  item <- ggnr_graphics$devices[[device]]
  value <- grDevices::dev.off()
  ggnr_graphics$devices[[device]] <- NULL
  if (!is.null(item) && nzchar(item$converter)) {
    embedded <- tempfile(tmpdir = dirname(item$file), fileext = ".pdf")
    log <- tempfile()
    status <- system2(item$converter,
      c("-pdf", shQuote(item$file), shQuote(embedded)), stdout = log, stderr = log)
    if (status != 0 || !file.exists(embedded) || file.info(embedded)$size < 100) {
      stop("Font embedding failed: ", paste(readLines(log, warn = FALSE), collapse = " "))
    }
    if (!file.copy(embedded, item$file, overwrite = TRUE)) stop("Cannot save the embedded PDF.")
    unlink(c(embedded, log))
  }
  invisible(value)
}
