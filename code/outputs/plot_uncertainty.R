# Diagnostic plots for Hamilton-style filtering-plus-parameter uncertainty.

input_directory <- Sys.getenv(
  "GGNR_HAMILTON_DIR",
  "outputs/diagnostics/uncertainty"
)
output_directory <- Sys.getenv("GGNR_FIGURE_DIR", "outputs/figures")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
states <- read.csv(file.path(input_directory, "state_bands.csv"))
premia <- read.csv(file.path(input_directory, "premium_bands.csv"))
states$date <- as.Date(states$date)
premia$date <- as.Date(premia$date)

draw_band <- function(data, prefix, title, zero = FALSE) {
  y95 <- c(data[[paste0(prefix, "_p025")]],
           data[[paste0(prefix, "_p975")]])
  limits <- range(y95[is.finite(y95)])
  plot(data$date, data[[paste0(prefix, "_estimate")]], type = "n",
       xlab = "", ylab = "percentage points", main = title, ylim = limits)
  polygon(c(data$date, rev(data$date)),
          c(data[[paste0(prefix, "_p025")]],
            rev(data[[paste0(prefix, "_p975")]])),
          col = "grey85", border = NA)
  polygon(c(data$date, rev(data$date)),
          c(data[[paste0(prefix, "_p16")]],
            rev(data[[paste0(prefix, "_p84")]])),
          col = "grey65", border = NA)
  if (zero) abline(h = 0, col = "grey45", lty = 3, lwd = 1.2)
  lines(data$date, data[[paste0(prefix, "_estimate")]], lwd = 2.5,
        col = "black")
}

draw_state_figure <- function() {
  par(mfrow = c(2, 1), mar = c(3.4, 5.1, 3.0, 1.0),
      cex.axis = 1.12, cex.lab = 1.18, cex.main = 1.28,
      mgp = c(2.8, 0.8, 0), tcl = -0.3)
  draw_band(states, "r_star", expression(paste("(a) Real-rate trend ", r*"*")),
            zero = TRUE)
  legend("topright", c("Estimate", "68% interval", "95% interval"),
         lty = c(1, NA, NA), lwd = c(2.5, NA, NA),
         pch = c(NA, 15, 15), pt.cex = c(NA, 2.0, 2.0),
         col = c("black", "grey65", "grey85"), bty = "n", cex = 1.05)
  draw_band(states, "pi_star",
            expression(paste("(b) Inflation trend ", pi*"*")), zero = TRUE)
}

draw_premium_figure <- function() {
  par(mfrow = c(3, 1), mar = c(3.4, 5.1, 3.0, 1.0),
      cex.axis = 1.12, cex.lab = 1.18, cex.main = 1.28,
      mgp = c(2.8, 0.8, 0), tcl = -0.3)
  par(cex = 1) # Restore full-size text after the multi-panel layout.
  draw_band(premia, "nominal_10y", "(a) 10-year nominal term premium", TRUE)
  legend("topright", c("Estimate", "68% interval", "95% interval"),
         lty = c(1, NA, NA), lwd = c(2.5, NA, NA),
         pch = c(NA, 15, 15), pt.cex = c(NA, 2.0, 2.0),
         col = c("black", "grey65", "grey85"), bty = "n", cex = 1.05)
  draw_band(premia, "real_10y", "(b) 10-year real term premium", TRUE)
  draw_band(premia, "inflation_10y", "(c) 10-year inflation risk premium", TRUE)
}

png(file.path(output_directory, "hamilton_state_bands.png"),
    width = 1800, height = 1100, res = 180)
draw_state_figure()
dev.off()
pdf(file.path(output_directory, "hamilton_state_bands.pdf"),
    width = 10, height = 6.1, useDingbats = FALSE)
draw_state_figure()
dev.off()

png(file.path(output_directory, "hamilton_premium_bands.png"),
    width = 1800, height = 1500, res = 180)
draw_premium_figure()
dev.off()
pdf(file.path(output_directory, "hamilton_premium_bands.pdf"),
    width = 10, height = 8.3, useDingbats = FALSE)
draw_premium_figure()
dev.off()
