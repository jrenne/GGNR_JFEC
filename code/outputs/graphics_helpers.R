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

# Reserve a separate top margin for legends, outside the plotting rectangle.
ggnr_panel <- function(fun, args) {
  main <- args$main; args$main <- NULL
  mar <- par("mar"); mar[3] <- max(mar[3], if (is.null(main)) 3.4 else 6.4)
  par(mar = mar)
  value <- do.call(fun, args)
  if (!is.null(main)) graphics::title(main = main, line = mar[3] - 1.6)
  invisible(value)
}
ggnr_plot <- function(...) ggnr_panel(graphics::plot, list(...))
ggnr_matplot <- function(...) ggnr_panel(graphics::matplot, list(...))
ggnr_legend <- function(x, ..., y = NULL) {
  args <- list(...); args$x <- args$y <- args$inset <- args$plot <- NULL
  args$xpd <- NA; args$bty <- "n"; args$bg <- NULL; args$box.col <- NULL
  n <- length(args$legend)
  if (is.null(args$legend) && length(args)) {
    args$legend <- args[[1]]; args[[1]] <- NULL; n <- length(args$legend)
  }
  if (is.null(args$cex)) args$cex <- 1
  usr <- par("usr"); width <- diff(usr[1:2])
  requested <- if (is.null(args$ncol)) min(3L, n) else args$ncol
  for (cols in seq.int(requested, 1L)) {
    args$ncol <- cols
    box <- do.call(graphics::legend, c(list(x = "top", plot = FALSE), args))
    if (box$rect$w <= .98 * width) break
  }
  if (box$rect$w > .98 * width) {
    args$cex <- args$cex * .98 * width / box$rect$w
    box <- do.call(graphics::legend, c(list(x = "top", plot = FALSE), args))
  }
  gap <- diff(usr[3:4]) * .10 / par("pin")[2]
  do.call(graphics::legend, c(list(x = mean(usr[1:2]),
    y = usr[4] + box$rect$h + gap, xjust = .5, yjust = 1), args))
}
