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

#' Parse a single fitqtl summary text file.
#' Returns list with full_model_pctvar and drop_one (data frame: term, pctvar).
parse_fitqtl_summary_file <- function(path) {
  lines <- readLines(path, warn = FALSE)
  full_pctvar <- NA_real_
  drop_one <- data.frame(term = character(), pctvar = numeric(), stringsAsFactors = FALSE)

  # Full model: find "Full model result" then the "Model" ANOVA table row (Model + df + SS + ...)
  in_full <- FALSE
  for (i in seq_along(lines)) {
    if (grepl("Full model result", lines[i], fixed = TRUE)) in_full <- TRUE
    if (in_full && grepl("^Model\\s+", lines[i])) {
      parts <- strsplit(trimws(lines[i]), "\\s+")[[1]]
      # Table row has Model, df (integer), SS, MS, LOD, %var; skip "Model formula:" and "Model:  normal"
      if (length(parts) >= 6 && !is.na(suppressWarnings(as.numeric(parts[2]))))
        full_pctvar <- as.numeric(parts[6])
      if (!is.na(full_pctvar)) break
    }
  }

  # Drop-one table: between "Drop one QTL at a time" and the short "---" line (not "----...")
  in_drop <- FALSE
  for (i in seq_along(lines)) {
    if (grepl("Drop one QTL at a time", lines[i], fixed = TRUE)) {
      in_drop <- TRUE
      next
    }
    if (in_drop && trimws(lines[i]) == "---") break
    if (in_drop && nzchar(trimws(lines[i]))) {
      parts <- strsplit(trimws(lines[i]), "\\s+")[[1]]
      # Skip header line
      if (identical(parts[1], "df") || identical(parts[1], "Model")) next
      if (length(parts) >= 5 && !is.na(suppressWarnings(as.numeric(parts[2])))) {
        term <- parts[1]
        pctvar <- as.numeric(parts[5])  # %var is 5th column in drop-one table
        if (!is.na(pctvar)) drop_one <- rbind(drop_one, data.frame(term = term, pctvar = pctvar, stringsAsFactors = FALSE))
      }
    }
  }

  list(full_model_pctvar = full_pctvar, drop_one = drop_one)
}

#' Classify drop-one term into component: QTL, other_genotypic, environment.
classify_term <- function(term) {
  term <- trimws(term)
  if (grepl("^[0-9]+@[0-9.]+$", term)) return("QTL")           # e.g. 3@59.2, 1@4.1
  if (term == "hybrid_index") return("other_genotypic")
  if (grepl("^F21_", term)) return("other_genotypic")
  if (grepl("^coh[0-9]+$", term) || grepl("^s\\.coh", term)) return("environment")
  NA_character_
}

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
plot_data$trait <- factor(plot_data$trait, levels = unique(plot_data$trait))
plot_data$component <- factor(plot_data$component,
  levels = c("Unexplained", "QTL", "Family and/or hybrid index", "Environment"))

# Plot: trait on x, % variance on y, one point per (trait, component), colored by component
# Colors: Wong/Nature Methods colorblind-friendly palette
p <- ggplot(plot_data, aes(x = trait, y = pct_var, colour = component, shape = component)) +
  geom_point(size = 3.5) +
  scale_colour_manual(
    values = c(
      "Unexplained"               = "#666666",
      "QTL"                       = "#E69F00",
      "Family and/or hybrid index" = "#56B4E9",
      "Environment"               = "#009E73"
    ),
    name = "Variance component"
  ) +
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
    x = "Trait",
    y = "% variance"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) +
  ylim(0, NA)

print(p)

# Save
out_pdf <- "QTL_model_variances_plot.pdf"
out_png <- "QTL_model_variances_plot.png"
ggsave(out_pdf, p, width = 5, height = 3.5)
ggsave(out_png, p, width = 5, height = 3.5, dpi = 150)
message("Saved: ", out_pdf, ", ", out_png)

# Optional: save summary table
write.csv(plot_data, "QTL_model_variances_by_trait.csv", row.names = FALSE)
message("Saved: QTL_model_variances_by_trait.csv")
