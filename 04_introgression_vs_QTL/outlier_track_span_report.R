# =============================================================================
# HTML report: Hudson Fst and Dxy outlier-track spans (95th and 99th).
#
# Two track definitions per percentile:
#   adjacent: consecutive 10 kb outlier tiles (starts 10 kb apart).
#   gap10:    outliers may be bridged by up to 10 intervening 10 kb slots that
#             are missing (absent from the pixy table or non-finite). A finite
#             below-threshold window breaks the track.
#
# Chromosome 4 is excluded from counts. Thresholds use the genome-wide
# distribution including chromosome 4. QTL columns use the seven retained
# 95% CIs (Nreabs chr2 and GR chr5 dropped).
#
# Output: 04_introgression_vs_QTL/outlier_track_span_report.html
# =============================================================================

suppressPackageStartupMessages(library(tidyverse))

for (p in c(
  file.path("scripts", "project_paths.R"),
  file.path("..", "scripts", "project_paths.R")
)) {
  if (file.exists(p)) {
    source(p)
    break
  }
}
if (!exists("project_paths")) {
  stop("Could not source scripts/project_paths.R. Run from the repo root.")
}

paths <- project_paths()
project_root <- paths$project_root
pixy_folder <- paths$pixy_folder
qtl_ci_path <- paths$qtl_ci_path
if (is.na(qtl_ci_path) || !nzchar(qtl_ci_path) || !file.exists(qtl_ci_path)) {
  alt <- file.path(project_root, "QTL", "all_QTL_CIs_bp_cm.251130.csv")
  if (file.exists(alt)) qtl_ci_path <- alt
}
if (!file.exists(pixy_folder)) stop("Pixy folder not found: ", pixy_folder)
if (!file.exists(qtl_ci_path)) stop("QTL CI table not found: ", qtl_ci_path)

out_html <- file.path(
  project_root, "04_introgression_vs_QTL", "outlier_track_span_report.html"
)
WINDOW_BP <- 10000L
MAX_MISSING_BRIDGE <- 10L

qtl <- read.csv(qtl_ci_path, stringsAsFactors = FALSE) %>%
  filter(!(chr == 2 & trait == "Nreabs")) %>%
  filter(!(chr == 5 & trait == "GR")) %>%
  filter(as.integer(chr) != 4L) %>%
  mutate(
    chrN = as.integer(chr),
    label = paste0(trait, " chr", chrN)
  )

read_pixy_stat <- function(pattern, value_col) {
  files <- list.files(pixy_folder, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) stop("No files matching ", pattern, " in ", pixy_folder)
  bind_rows(lapply(files, read.delim, sep = "\t", stringsAsFactors = FALSE)) %>%
    mutate(
      chrN = as.integer(gsub("[^0-9]", "", chromosome)),
      value = suppressWarnings(as.numeric(.data[[value_col]])),
      window_pos_1 = as.numeric(window_pos_1),
      window_pos_2 = as.numeric(window_pos_2)
    ) %>%
    select(chrN, window_pos_1, window_pos_2, value)
}

# Build tracks on one chromosome.
# mode = "adjacent": only merge when the next outlier starts exactly 10 kb later.
# mode = "gap10": merge when every intervening 10 kb slot is missing and the
# number of intervening slots is between 1 and MAX_MISSING_BRIDGE (0 is adjacent).
tracks_one_chr <- function(df_chr, thr, mode = c("adjacent", "gap10")) {
  mode <- match.arg(mode)
  df_chr <- df_chr %>% arrange(window_pos_1)
  outs <- which(is.finite(df_chr$value) & df_chr$value >= thr)
  if (length(outs) == 0L) {
    return(tibble(
      start = numeric(), end = numeric(), n_outlier_windows = integer(),
      n_slots_spanned = integer()
    ))
  }

  # For each pair of successive outlier rows in the table, count intervening
  # 10 kb slots and whether any finite below-threshold window sits between them.
  starts <- df_chr$window_pos_1[outs]
  ends <- df_chr$window_pos_2[outs]
  n_out <- length(outs)
  can_join <- logical(max(0L, n_out - 1L))
  if (n_out > 1L) {
    for (i in seq_len(n_out - 1L)) {
      a <- starts[i]
      b <- starts[i + 1L]
      slots_between <- as.integer((b - a) / WINDOW_BP) - 1L
      if (slots_between < 0L || (b - a) %% WINDOW_BP != 0) {
        can_join[i] <- FALSE
        next
      }
      if (mode == "adjacent") {
        can_join[i] <- slots_between == 0L
        next
      }
      # gap10: 0..MAX_MISSING_BRIDGE missing slots; no finite non-outlier between.
      if (slots_between > MAX_MISSING_BRIDGE) {
        can_join[i] <- FALSE
        next
      }
      if (slots_between == 0L) {
        can_join[i] <- TRUE
        next
      }
      between_idx <- seq.int(outs[i] + 1L, outs[i + 1L] - 1L)
      between_vals <- df_chr$value[between_idx]
      # Any finite window between the two outliers is a real observation and
      # must be below threshold (otherwise it would itself be an outlier). A
      # finite below-threshold window breaks the bridge. Completely absent
      # slots are inferred from the bp gap when the table skips them.
      present_slots <- length(between_idx)
      absent_slots <- slots_between - present_slots
      if (absent_slots < 0L) {
        can_join[i] <- FALSE
        next
      }
      has_finite_non_outlier <- any(is.finite(between_vals))
      can_join[i] <- !has_finite_non_outlier
    }
  }

  run_id <- cumsum(c(TRUE, !can_join))
  tibble(
    start = starts,
    end = ends,
    n_outlier_windows = 1L,
    run_id = run_id
  ) %>%
    group_by(run_id) %>%
    summarise(
      start = min(start),
      end = max(end),
      n_outlier_windows = sum(n_outlier_windows),
      n_slots_spanned = as.integer((max(end) - min(start) + 1L) / WINDOW_BP),
      .groups = "drop"
    ) %>%
    select(-run_id)
}

# Completely absent tiles: expand each chromosome to a regular 10 kb grid so
# gap counting uses missing slots, not only rows present as NA.
expand_chr_grid <- function(df_chr) {
  df_chr <- df_chr %>% arrange(window_pos_1)
  if (nrow(df_chr) == 0L) return(df_chr)
  min_start <- min(df_chr$window_pos_1)
  max_end <- max(df_chr$window_pos_2)
  grid_starts <- seq(min_start, max_end - WINDOW_BP + 1L, by = WINDOW_BP)
  # Prefer starts that align to the observed tiling (usually 1, 10001, ...).
  # If the chromosome does not begin at 1, keep the observed minimum.
  tibble(
    chrN = df_chr$chrN[1],
    window_pos_1 = grid_starts,
    window_pos_2 = grid_starts + WINDOW_BP - 1L
  ) %>%
    left_join(
      df_chr %>% select(window_pos_1, value),
      by = "window_pos_1"
    )
}

build_tracks <- function(win, thr, mode) {
  bind_rows(lapply(split(win, win$chrN), function(df_chr) {
    if (mode == "gap10") df_chr <- expand_chr_grid(df_chr)
    tr <- tracks_one_chr(df_chr, thr, mode)
    if (nrow(tr) == 0L) {
      return(tr %>% mutate(chrN = integer()))
    }
    tr %>% mutate(chrN = df_chr$chrN[1])
  }))
}

annotate_qtl <- function(tracks) {
  if (nrow(tracks) == 0L) {
    return(tracks %>% mutate(in_qtl = FALSE, qtls = NA_character_))
  }
  tracks %>%
    rowwise() %>%
    mutate(
      hit = {
        labs <- qtl$label[
          qtl$chrN == chrN &
            start <= qtl$upper95CI_bp &
            end >= qtl$lower95CI_bp
        ]
        if (length(labs) == 0L) NA_character_ else paste(sort(unique(labs)), collapse = ", ")
      }
    ) %>%
    ungroup() %>%
    mutate(
      in_qtl = !is.na(hit),
      qtls = hit
    ) %>%
    select(-hit)
}

summarise_tracks <- function(tracks) {
  if (nrow(tracks) == 0L) {
    return(tibble(
      n_slots_spanned = integer(),
      span_kb = integer(),
      n_tracks = integer(),
      n_tracks_in_qtl = integer(),
      qtl_cis = character()
    ))
  }
  tracks %>%
    group_by(n_slots_spanned) %>%
    summarise(
      n_tracks = n(),
      n_tracks_in_qtl = sum(in_qtl),
      qtl_cis = {
        labs <- sort(unique(unlist(strsplit(na.omit(qtls), ", "))))
        if (length(labs) == 0L) "—" else paste(labs, collapse = ", ")
      },
      .groups = "drop"
    ) %>%
    mutate(span_kb = as.integer(n_slots_spanned) * 10L) %>%
    relocate(n_slots_spanned, span_kb, n_tracks, n_tracks_in_qtl, qtl_cis) %>%
    arrange(n_slots_spanned)
}

html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

table_html <- function(tab, caption) {
  if (nrow(tab) == 0L) {
    return(paste0(
      "<table><caption>", html_escape(caption),
      "</caption><tbody><tr><td>No outlier tracks</td></tr></tbody></table>"
    ))
  }
  rows <- character(nrow(tab))
  for (i in seq_len(nrow(tab))) {
    rows[i] <- paste0(
      "<tr>",
      "<td>", tab$n_slots_spanned[i], "</td>",
      "<td>", tab$span_kb[i], "</td>",
      "<td>", tab$n_tracks[i], "</td>",
      "<td>", tab$n_tracks_in_qtl[i], "</td>",
      "<td>", html_escape(tab$qtl_cis[i]), "</td>",
      "</tr>"
    )
  }
  paste0(
    "<table>\n",
    "<caption>", html_escape(caption), "</caption>\n",
    "<thead><tr>",
    "<th>Spanned 10 kb slots</th>",
    "<th>Span (kb)</th>",
    "<th>Tracks</th>",
    "<th>Tracks in a QTL CI</th>",
    "<th>QTL CIs</th>",
    "</tr></thead>\n<tbody>\n",
    paste(rows, collapse = "\n"),
    "\n</tbody></table>"
  )
}

analyse_stat <- function(win_all, label, value_name) {
  thr99 <- as.numeric(quantile(win_all$value, 0.99, na.rm = TRUE))
  thr95 <- as.numeric(quantile(win_all$value, 0.95, na.rm = TRUE))
  win <- win_all %>% filter(chrN != 4L)
  specs <- list(
    list(prob = 0.99, thr = thr99, mode = "adjacent", title = paste0(label, " 99th percentile — adjacent tracks")),
    list(prob = 0.99, thr = thr99, mode = "gap10", title = paste0(label, " 99th percentile — bridge up to 10 missing windows")),
    list(prob = 0.95, thr = thr95, mode = "adjacent", title = paste0(label, " 95th percentile — adjacent tracks")),
    list(prob = 0.95, thr = thr95, mode = "gap10", title = paste0(label, " 95th percentile — bridge up to 10 missing windows"))
  )
  tables <- lapply(specs, function(sp) {
    tr <- build_tracks(win, sp$thr, sp$mode) %>% annotate_qtl()
    list(
      title = sp$title,
      thr = sp$thr,
      mode = sp$mode,
      prob = sp$prob,
      n_outlier_windows = sum(is.finite(win$value) & win$value >= sp$thr),
      n_tracks = nrow(tr),
      summary = summarise_tracks(tr),
      tracks = tr
    )
  })
  list(
    label = label,
    value_name = value_name,
    thr99 = thr99,
    thr95 = thr95,
    tables = tables
  )
}

cat("Reading pixy Fst and Dxy...\n")
fst_win <- read_pixy_stat("_fst\\.txt$", "avg_hudson_fst")
dxy_win <- read_pixy_stat("_dxy\\.txt$", "avg_dxy")

cat("Building Fst track tables...\n")
fst_res <- analyse_stat(fst_win, "Hudson Fst", "avg_hudson_fst")
cat("Building Dxy track tables...\n")
dxy_res <- analyse_stat(dxy_win, "Dxy", "avg_dxy")

page_html <- function(res) {
  blocks <- vapply(res$tables, function(tb) {
    note <- if (tb$mode == "adjacent") {
      "Tracks are runs of outlier 10 kb windows that start exactly 10 kb apart."
    } else {
      paste0(
        "Tracks also merge outliers separated by 1–", MAX_MISSING_BRIDGE,
        " intervening 10 kb slots that are missing (absent or non-finite). ",
        "A finite below-threshold window breaks the track. Span is the full ",
        "base-pair distance from the first to last window in the track."
      )
    }
    paste0(
      "<section class=\"table-block\">\n",
      "<h2>", html_escape(tb$title), "</h2>\n",
      "<p class=\"meta\">Cutoff = ", signif(tb$thr, 4),
      " (genome-wide ", tb$prob * 100, "th percentile of ",
      html_escape(res$value_name), ", including chr 4). ",
      "Outlier windows on chr ≠ 4: ", tb$n_outlier_windows,
      ". Tracks: ", tb$n_tracks, ". ", note, "</p>\n",
      table_html(tb$summary, tb$title),
      "\n</section>\n"
    )
  }, character(1))
  paste(blocks, collapse = "\n")
}

html <- paste0(
  "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n",
  "<title>Outlier track spans — Fst and Dxy</title>\n",
  "<style>\n",
  ":root { --ink:#1a1a1a; --muted:#555; --line:#ccc; --bg:#fafafa; --accent:#245b4a; }\n",
  "body { margin:0; font-family: 'Iowan Old Style', 'Palatino Linotype', Palatino, Georgia, serif;",
  " color:var(--ink); background:var(--bg); line-height:1.45; }\n",
  "header { background:#fff; border-bottom:1px solid var(--line); padding:1.2rem 1.5rem 0.8rem; }\n",
  "header h1 { margin:0 0 0.35rem; font-size:1.55rem; font-weight:600; }\n",
  "header p { margin:0; color:var(--muted); font-size:0.95rem; max-width:52rem; }\n",
  "nav { display:flex; gap:0.5rem; margin-top:0.9rem; }\n",
  "nav button { font:inherit; border:1px solid var(--line); background:#fff; padding:0.4rem 0.9rem;",
  " border-radius:4px; cursor:pointer; color:var(--ink); }\n",
  "nav button[aria-selected='true'] { background:var(--accent); color:#fff; border-color:var(--accent); }\n",
  "main { padding:1.25rem 1.5rem 2.5rem; max-width:70rem; }\n",
  ".page { display:none; }\n",
  ".page.active { display:block; }\n",
  "h2 { font-size:1.15rem; margin:1.6rem 0 0.4rem; }\n",
  ".meta { color:var(--muted); font-size:0.92rem; margin:0 0 0.7rem; max-width:52rem; }\n",
  "table { border-collapse:collapse; width:100%; background:#fff; margin-bottom:0.5rem;",
  " font-variant-numeric:tabular-nums; }\n",
  "caption { text-align:left; font-weight:600; padding:0.4rem 0; caption-side:top; }\n",
  "th, td { border:1px solid var(--line); padding:0.4rem 0.55rem; text-align:left; vertical-align:top; }\n",
  "th { background:#f0f3f1; font-weight:600; }\n",
  "td:nth-child(1), td:nth-child(2), td:nth-child(3), td:nth-child(4) { text-align:right; width:7rem; }\n",
  "footer { color:var(--muted); font-size:0.85rem; padding:0 1.5rem 2rem; max-width:70rem; }\n",
  "@media print { nav { display:none; } .page { display:block !important; page-break-after:always; } }\n",
  "</style>\n</head>\n<body>\n",
  "<header>\n",
  "<h1>Outlier track spans</h1>\n",
  "<p>Hudson <i>F</i><sub>st</sub> and <i>D</i><sub>xy</sub> 10 kb outlier tracks at the 95th and 99th percentiles. ",
  "Chromosome 4 is excluded from track counts. QTL columns use the seven retained 95% credible intervals ",
  "(Nreabs chr2 and GR chr5 dropped). A track in a QTL CI means the track interval intersects a CI.</p>\n",
  "<nav role=\"tablist\" aria-label=\"Statistic\">\n",
  "<button role=\"tab\" id=\"tab-fst\" aria-selected=\"true\" aria-controls=\"page-fst\">Hudson Fst</button>\n",
  "<button role=\"tab\" id=\"tab-dxy\" aria-selected=\"false\" aria-controls=\"page-dxy\">Dxy</button>\n",
  "</nav>\n</header>\n",
  "<main>\n",
  "<div class=\"page active\" id=\"page-fst\" role=\"tabpanel\" aria-labelledby=\"tab-fst\">\n",
  page_html(fst_res),
  "</div>\n",
  "<div class=\"page\" id=\"page-dxy\" role=\"tabpanel\" aria-labelledby=\"tab-dxy\">\n",
  page_html(dxy_res),
  "</div>\n",
  "</main>\n",
  "<footer>Generated by <code>04_introgression_vs_QTL/outlier_track_span_report.R</code>. ",
  "Adjacent mode requires successive outlier windows to start 10 kb apart. ",
  "Gap mode allows 1–10 missing intervening 10 kb slots; finite below-threshold windows break tracks.</footer>\n",
  "<script>\n",
  "const tabs = [...document.querySelectorAll('nav button')];\n",
  "const pages = { 'tab-fst': 'page-fst', 'tab-dxy': 'page-dxy' };\n",
  "tabs.forEach(btn => btn.addEventListener('click', () => {\n",
  "  tabs.forEach(b => b.setAttribute('aria-selected', b === btn ? 'true' : 'false'));\n",
  "  Object.values(pages).forEach(id => document.getElementById(id).classList.remove('active'));\n",
  "  document.getElementById(pages[btn.id]).classList.add('active');\n",
  "}));\n",
  "</script>\n</body>\n</html>\n"
)

writeLines(html, out_html)
cat("Wrote ", out_html, "\n", sep = "")

# Also print a compact console summary for the adjacent 99th tables.
for (res in list(fst_res, dxy_res)) {
  cat("\n=== ", res$label, " ===\n", sep = "")
  cat("99th thr ", signif(res$thr99, 4), "; 95th thr ", signif(res$thr95, 4), "\n", sep = "")
  for (tb in res$tables) {
    cat("\n", tb$title, " | tracks ", tb$n_tracks, "\n", sep = "")
    print(as.data.frame(tb$summary), row.names = FALSE)
  }
}
