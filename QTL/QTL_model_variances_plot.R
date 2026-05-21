# QTL_model_variances_plot.R
# Reads *QTL_model_output_251203 summary (fitqtl) text files and plots
# % variance explained by component: unexplained, QTL, other genotypic, environment.
# Run from project root or with setwd to QTL/.

library(dplyr)
library(ggplot2)

# Working directory: use QTL/ if files live there
if (!file.exists("GR_QTL_model_output_251203") &&
    file.exists(file.path("QTL", "GR_QTL_model_output_251203"))) {
  setwd("QTL")
  message("Set working directory to QTL/")
}

# Find all model output files
pattern <- "*_QTL_model_output_251203"
model_files <- list.files(pattern = glob2rx(pattern), full.names = FALSE)
if (length(model_files) == 0) {
  stop("No files matching ", pattern, " found in ", getwd())
}
message("Found ", length(model_files), " model output file(s).")

utils_path <- if (file.exists("fitqtl_summary_utils.R")) {
  "fitqtl_summary_utils.R"
} else {
  file.path("QTL", "fitqtl_summary_utils.R")
}
source(utils_path, local = TRUE)
classify_term <- classify_fitqtl_term

# Process each file and build long data for plotting
plot_data_list <- list()
for (f in model_files) {
  trait <- sub("_QTL_model_output_251203$", "", f)
  parsed <- parse_fitqtl_summary_file(f)
  full_pct <- parsed$full_model_pctvar
  unexplained <- 100 - full_pct
  if (is.na(full_pct)) {
    warning("Could not parse full model %var from ", f)
    next
  }

  # Sum drop-one %var by component
  do <- parsed$drop_one
  do$component <- vapply(do$term, classify_term, character(1))
  do <- do[!is.na(do$component), ]

  sums <- do %>%
    group_by(component) %>%
    summarise(pct_var = sum(pctvar, na.rm = TRUE), .groups = "drop")

  qtl_var <- if (any(sums$component == "QTL")) sums$pct_var[sums$component == "QTL"] else 0
  other_var <- if (any(sums$component == "other_genotypic")) sums$pct_var[sums$component == "other_genotypic"] else 0
  env_var <- if (any(sums$component == "environment")) sums$pct_var[sums$component == "environment"] else 0

  # One row per (trait, component); include unexplained
  plot_data_list[[trait]] <- data.frame(
    trait = trait,
    component = c("Unexplained", "QTL", "Family and/or hybrid index", "Environment"),
    pct_var = c(unexplained, qtl_var, other_var, env_var),
    stringsAsFactors = FALSE
  )
}

plot_data <- bind_rows(plot_data_list)

# Order traits (use order of appearance in files, or alphabetically)
trait_levels <- unique(plot_data$trait)
plot_data$trait <- factor(plot_data$trait, levels = trait_levels)
plot_data$component <- factor(plot_data$component,
  levels = c("Unexplained", "QTL", "Family and/or hybrid index", "Environment"))

# Display labels for x-axis (abbreviations with periods)
trait_labels <- setNames(as.character(trait_levels), trait_levels)
trait_labels["Chlreabs"] <- "Chl.reabs"
trait_labels["Nreabs"]   <- "N.reabs"

# Plot: trait on x, % variance on y. Draw Environment (squares) first, then Family (triangles) on top.
# All points black; Unexplained = open circle (shape 1), others filled.
p <- ggplot(plot_data, aes(x = trait, y = pct_var, shape = component)) +
  geom_point(data = filter(plot_data, component != "Family and/or hybrid index"), size = 3.5, colour = "black") +
  geom_point(data = filter(plot_data, component == "Family and/or hybrid index"), size = 3.5, colour = "black") +
  scale_shape_manual(
    values = c(
      "Unexplained"               = 1,
      "QTL"                       = 16,
      "Family and/or hybrid index" = 17,
      "Environment"               = 15
    ),
    name = "Variance component"
  ) +
  labs(
    x = NULL,
    y = "% phenotypic variance"
  ) +
  scale_x_discrete(labels = trait_labels) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 0),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) +
  ylim(0, NA)

print(p)

# Save
out_pdf <- "QTL_model_variances_plot.pdf"
out_png <- "QTL_model_variances_plot.png"
ggsave(out_pdf, p, width = 5, height = 2.5)
ggsave(out_png, p, width = 5, height = 2.5, dpi = 150)
message("Saved: ", out_pdf, ", ", out_png)

# Optional: save summary table
write.csv(plot_data, "QTL_model_variances_by_trait.csv", row.names = FALSE)
message("Saved: QTL_model_variances_by_trait.csv")
