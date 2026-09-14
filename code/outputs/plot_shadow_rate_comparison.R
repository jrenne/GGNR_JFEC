# Compare the retained shadow-rate estimate with two external model-based series.
# The figure has no overall title because its interpretation is carried by the
# self-contained LaTeX caption.

state_file <- "outputs/figures/states_and_shadow_rate.csv"
comparison_file <- "data/external/shadow_rate_comparisons.rds"
output_file <- "outputs/figures/shadow_rate_benchmarks.pdf"

model <- read.csv(state_file)
model$date <- as.Date(model$date)

comparison <- readRDS(comparison_file)
comparison$date <- as.Date(comparison$date)

series <- Reduce(
  function(x, y) merge(x, y, by = "date", all = TRUE),
  list(
    model[c("date", "shadow_rate", "effective_federal_funds_rate")],
    comparison
  )
)

draw_panel <- function(from, to, panel_title, show_legend = FALSE) {
  keep <- series$date >= as.Date(from) & series$date <= as.Date(to)
  z <- series[keep, ]
  limits <- range(
    z[c("shadow_rate", "wu_xia", "krippner", "effective_federal_funds_rate")],
    na.rm = TRUE
  )
  plot(
    z$date, z$shadow_rate, type = "l", col = "black", lwd = 2.8,
    xlab = "", ylab = "Percent", main = panel_title, ylim = limits
  )
  lines(z$date, z$wu_xia, col = "grey35", lwd = 2.2, lty = 2)
  lines(z$date, z$krippner, col = "grey60", lwd = 2.2, lty = 3)
  lines(
    z$date, z$effective_federal_funds_rate,
    col = "grey75", lwd = 1.8, lty = 4
  )
  abline(h = 0, col = "grey78", lwd = 0.8)
  grid(col = "grey90", lty = 1)
  if (show_legend) {
    legend(
      "bottomleft",
      c("Model", "Wu--Xia", "Krippner", "Effective federal funds rate"),
      col = c("black", "grey35", "grey60", "grey75"),
      lty = c(1, 2, 3, 4), lwd = c(2.8, 2.2, 2.2, 1.8),
      cex = 1.02, bg = adjustcolor("white", alpha.f = 0.90)
    )
  }
}

pdf(output_file, width = 10, height = 7.5, useDingbats = FALSE)
par(
  mfrow = c(2, 1), mar = c(3.5, 4.8, 2.8, 0.8),
  cex.axis = 1.15, cex.lab = 1.20, cex.main = 1.18,
  font.main = 2, mgp = c(2.6, 0.8, 0), tcl = -0.35
)
draw_panel("1995-01-01", "2023-06-01", "(a) Common sample", TRUE)
draw_panel("2008-01-01", "2022-03-01", "(b) Effective-lower-bound periods")
dev.off()

# Compact diagnostics used to keep the accompanying discussion quantitative.
month_key <- function(x) format(as.Date(x), "%Y-%m")
comparison <- merge(
  data.frame(month = month_key(model$date), model = model$shadow_rate),
  merge(
    data.frame(month = month_key(comparison$date),
               wu_xia = comparison$wu_xia),
    data.frame(month = month_key(comparison$date),
               krippner = comparison$krippner),
    by = "month"
  ),
  by = "month"
)
diagnostics <- function(z, alternative) {
  keep <- is.finite(z$model) & is.finite(z[[alternative]])
  model <- z$model[keep]
  comparison <- z[[alternative]][keep]
  difference <- model - comparison
  c(
    correlation = cor(model, comparison),
    mean_absolute_difference = mean(abs(difference)),
    root_mean_squared_difference = sqrt(mean(difference^2))
  )
}
print(rbind(
  wu_xia = diagnostics(comparison, "wu_xia"),
  krippner = diagnostics(comparison, "krippner")
))
cat("Saved", output_file, "\n")
