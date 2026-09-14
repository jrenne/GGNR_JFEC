# Optional reconstruction of the estimation database
#
# Open this file in RStudio, choose the last month of the updated sample, and
# click Source. The frozen data used in the paper are never modified.

sample_start <- as.Date("1968-10-01")
sample_end <- as.Date("2026-08-01")

source("code/data/build_updated_data.R")

updated_data <- build_updated_data(sample_start, sample_end)
saveRDS(updated_data, "data/updated/updated_data.rds", version = 3)

message("Updated database written to data/updated/updated_data.rds")
message("The frozen paper data in data/input/ were not changed.")
