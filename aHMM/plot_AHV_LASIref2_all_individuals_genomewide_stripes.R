# Genome-wide aHMM ancestry: all individuals in one PDF, alternating
# black / dark grey points by chromosome. Requires CRAN package ggh4x (strip styling).
#
# From project root (directory that contains `aHMM/`):
#   Rscript aHMM/plot_AHV_LASIref2_all_individuals_genomewide_stripes.R
#
# When sourced from `aHMM/popgen_aHMM.qmd`, the working directory is usually `aHMM/`;
# paths are set automatically.

library(readr)
library(ggplot2)
library(dplyr)
library(grid)
library(ggh4x)

if (dir.exists("aHMM/AHV_LASIref2_260210")) {
  ahv_base <- "aHMM/AHV_LASIref2_260210"
  out_pdf <- "aHMM/sup_wild_aHMM.pdf"
} else if (dir.exists("AHV_LASIref2_260210")) {
  ahv_base <- "AHV_LASIref2_260210"
  out_pdf <- "AHV_LASIref2_260210_all_individuals_genomewide.pdf"
} else {
  stop(
    "Cannot find AHV_LASIref2_260210. Run from project root, or render/sourced with wd = aHMM/"
  )
}
chr_subdirs <- paste0("chr", 1:9)

read_ahv_posterior <- function(path) {
  x <- read_tsv(path, show_col_types = FALSE)
  names(x) <- c("chrom", "position", "ref_homo", "het", "alt_homo")
  x
}

posterior_to_plot <- function(x) {
  x %>%
    filter(ref_homo > 0.9 | alt_homo > 0.9 | het > 0.9) %>%
    mutate(ancestry = ref_homo + het * 0.5)
}

ind_files <- list.files(
  file.path(ahv_base, "chr1"),
  pattern = "\\.posterior$",
  full.names = FALSE
)
individuals <- sort(sub("\\.posterior$", "", ind_files))
# Strip labels: digits and dots only (narrower facet strip than e.g. 19.516_HYB)
id_strip_label <- gsub("[^0-9.]", "", individuals)

missing <- NULL
for (id in individuals) {
  for (sd in chr_subdirs) {
    p <- file.path(ahv_base, sd, paste0(id, ".posterior"))
    if (!file.exists(p)) missing <- c(missing, p)
  }
}
if (length(missing)) {
  stop("Missing posterior files:\n", paste(missing, collapse = "\n"))
}

ref_id <- individuals[[1]]
chr_starts <- list()
cum <- 0
chr_band_tbl <- NULL
for (i in seq_along(chr_subdirs)) {
  p <- file.path(ahv_base, chr_subdirs[[i]], paste0(ref_id, ".posterior"))
  d0 <- read_ahv_posterior(p)
  chrom <- as.character(unique(d0$chrom)[[1]])
  mx <- max(d0$position)
  point_col <- if (i %% 2L == 1L) "grey35" else "black"
  chr_band_tbl <- rbind(
    chr_band_tbl,
    data.frame(
      chrom = chrom,
      xmin = cum,
      xmax = cum + mx,
      point_col = point_col,
      stringsAsFactors = FALSE
    )
  )
  chr_starts[[chrom]] <- cum
  cum <- cum + mx
}

chr_start_vec <- unlist(chr_starts)

all_rows <- list()
for (id in individuals) {
  for (sd in chr_subdirs) {
    p <- file.path(ahv_base, sd, paste0(id, ".posterior"))
    d <- read_ahv_posterior(p) %>% posterior_to_plot()
    chrom <- as.character(d$chrom[[1]])
    d$genome_pos <- d$position + chr_start_vec[[chrom]]
    d$individual <- id
    all_rows[[length(all_rows) + 1L]] <- d
  }
}
plot_df <- bind_rows(all_rows)
id_to_strip <- stats::setNames(id_strip_label, individuals)
plot_df$individual <- factor(
  id_to_strip[as.character(plot_df$individual)],
  levels = id_strip_label
)
plot_df <- plot_df %>%
  left_join(chr_band_tbl %>% select(chrom, point_col), by = "chrom")

chr_band_tbl$chr_num <- seq_len(nrow(chr_band_tbl))
chr_band_tbl$mid <- (chr_band_tbl$xmin + chr_band_tbl$xmax) / 2
chr_boundaries <- chr_band_tbl$xmin[-1L]

row_bg <- data.frame(
  individual = factor(id_strip_label, levels = id_strip_label),
  bg = ifelse(seq_along(id_strip_label) %% 2L == 1L, "grey90", "white"),
  stringsAsFactors = FALSE
)

p <- ggplot(plot_df, aes(genome_pos, ancestry)) +
  geom_rect(
    data = row_bg,
    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = bg),
    inherit.aes = FALSE
  ) +
  scale_fill_identity() +
  geom_vline(
    xintercept = chr_boundaries,
    color = "grey55",
    linewidth = 0.35,
    inherit.aes = FALSE
  ) +
  geom_point(size = 1, stroke = 0, aes(color = point_col), alpha = 0.85) +
  scale_color_identity(guide = "none") +
  facet_grid2(
    rows = vars(individual),
    scales = "fixed",
    switch = "y",
    strip = strip_themed(
      background_y = list(
        element_rect(fill = "grey90", colour = NA),
        element_rect(fill = "white", colour = NA)
      )
    )
  ) +
  scale_x_continuous(
    breaks = chr_band_tbl$mid,
    labels = chr_band_tbl$chr_num,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    breaks = c(0, 0.5, 1),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Chromosome",
    y = "Ancestry",
    title = "wild hybrid aHMM posteriors"
  ) +
  theme_bw() +
  theme(
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    panel.spacing.y = unit(0.45, "lines"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = NA, colour = NA),
    plot.background = element_rect(fill = "white", colour = NA)
  )

# Single-column page width; height scales with n but stays within a typical full page.
fig_w_cm <- 17.8
n_ind <- length(individuals)
fig_h_cm <- min(26, max(14, n_ind * 1.9))

ggsave(
  out_pdf,
  plot = p,
  width = fig_w_cm,
  height = fig_h_cm,
  units = "cm",
  limitsize = FALSE
)
message("Wrote ", normalizePath(out_pdf, mustWork = FALSE))
