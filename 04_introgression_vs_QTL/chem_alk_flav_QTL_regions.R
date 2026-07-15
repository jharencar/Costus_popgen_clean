# Shared helper: alkaloid / flavonoid 10% FDR QTL regions (15 cM clustering).
# Source from chemistry overlay scripts under 04_introgression_vs_QTL/.

load_chem_alk_flav_regions <- function(project_root, cM_window = 15) {
  chem_path <- file.path(project_root, "03_QTL_mapping", "data", "Chem_QTL_mass_WATERS.csv")
  sig_path <- file.path(project_root, "QTL", "Chem_260703", "WatersMass_chemistry_significant_10pct.csv")
  summ_path <- file.path(project_root, "QTL", "Chem_260703", "WatersMass_chemistry_scanone_summary.csv")
  if (!file.exists(summ_path)) {
    summ_path <- file.path(project_root, "03_QTL_mapping", "data", "WatersMass_chemistry_scanone_summary.csv")
  }
  stopifnot(file.exists(chem_path), file.exists(sig_path), file.exists(summ_path))

  clean_npc <- function(x) {
    x <- trimws(as.character(x))
    x[x %in% c("#N/A", "nan", "0", "None", "NA", "")] <- ""
    x
  }
  parse_marker_bp <- function(marker) {
    as.numeric(sub("^Chrom[0-9]+_", "", as.character(marker)))
  }

  chem_raw <- read.csv(chem_path, check.names = FALSE)
  chem_raw <- chem_raw[, nzchar(names(chem_raw)) & !is.na(names(chem_raw)), drop = FALSE]
  chem <- chem_raw %>%
    dplyr::mutate(
      pathway = clean_npc(.data[["NPC#pathway"]]),
      superclass = clean_npc(.data[["NPC#superclass"]]),
      row_ID = as.integer(.data[["row ID"]]),
      trait = paste0("row", row_ID),
      group = dplyr::case_when(
        pathway == "Alkaloids" ~ "Alkaloid",
        superclass == "Flavonoids" ~ "Flavonoid",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(group))

  sig <- read.csv(sig_path, check.names = FALSE)
  summ <- read.csv(summ_path, check.names = FALSE)

  peaks <- summ %>%
    dplyr::filter(trait %in% chem$trait, trait %in% sig$trait) %>%
    dplyr::left_join(chem %>% dplyr::select(trait, group), by = "trait") %>%
    dplyr::mutate(
      peak_bp = parse_marker_bp(marker),
      chrN = as.integer(as.character(chr))
    ) %>%
    dplyr::filter(is.finite(peak_bp), !is.na(group))

  cluster_group <- function(df, cM = cM_window) {
    df <- df %>% dplyr::arrange(chrN, pos)
    clusters <- list()
    cur <- NULL
    for (i in seq_len(nrow(df))) {
      row <- df[i, ]
      if (is.null(cur) || row$chrN != cur$chrN || abs(row$pos - cur$pos_center) > cM) {
        if (!is.null(cur)) clusters[[length(clusters) + 1]] <- cur
        cur <- list(
          chrN = row$chrN,
          group = row$group,
          traits = row$trait,
          pos_min = row$pos,
          pos_max = row$pos,
          pos_center = row$pos,
          max_lod = row$lod,
          peak_bp = row$peak_bp,
          best = row$trait
        )
      } else {
        cur$traits <- c(cur$traits, row$trait)
        cur$pos_min <- min(cur$pos_min, row$pos)
        cur$pos_max <- max(cur$pos_max, row$pos)
        cur$pos_center <- (cur$pos_min + cur$pos_max) / 2
        if (row$lod > cur$max_lod) {
          cur$max_lod <- row$lod
          cur$peak_bp <- row$peak_bp
          cur$best <- row$trait
        }
      }
    }
    if (!is.null(cur)) clusters[[length(clusters) + 1]] <- cur
    dplyr::bind_rows(lapply(clusters, function(c) {
      tibble::tibble(
        chrN = c$chrN,
        group = c$group,
        n_traits = length(c$traits),
        sharing = if (length(c$traits) >= 2) "shared" else "independent",
        pos_min_cM = c$pos_min,
        pos_max_cM = c$pos_max,
        max_lod = c$max_lod,
        peak_bp = c$peak_bp,
        best_trait = c$best,
        traits = paste(c$traits, collapse = ", ")
      )
    }))
  }

  dplyr::bind_rows(
    cluster_group(peaks %>% dplyr::filter(group == "Alkaloid")),
    cluster_group(peaks %>% dplyr::filter(group == "Flavonoid"))
  ) %>%
    dplyr::mutate(
      group = factor(group, levels = c("Alkaloid", "Flavonoid")),
      sharing = factor(sharing, levels = c("shared", "independent"))
    )
}

chem_colors_alk_flav <- c(Alkaloid = "#D55E00", Flavonoid = "#0072B2")

add_chem_vlines_bp <- function(base_plot, chem_df) {
  # Per-chromosome plots use bp on x; hard-code colors so scales do not clash.
  p <- base_plot
  for (g in names(chem_colors_alk_flav)) {
    sub <- chem_df %>% dplyr::filter(as.character(group) == g)
    if (nrow(sub) == 0) next
    p <- p +
      ggplot2::geom_vline(
        data = sub,
        ggplot2::aes(xintercept = peak_bp, linetype = sharing),
        color = chem_colors_alk_flav[[g]],
        linewidth = 0.9,
        inherit.aes = FALSE
      )
  }
  p +
    ggplot2::scale_linetype_manual(
      values = c(shared = "solid", independent = "dashed"),
      drop = FALSE,
      guide = "none"
    )
}
