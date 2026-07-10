# Code to plot lepmap output as genetic map

for (p in c(
  file.path("scripts", "project_paths.R"),
  file.path("..", "..", "scripts", "project_paths.R")
)) {
  if (file.exists(p)) { source(p); break }
}
paths <- project_paths()
lepmap_wd <- Sys.getenv("LEPMAP_WD", "")
if (!nzchar(lepmap_wd)) {
  lepmap_wd <- if (!is.na(paths$lepmap_dir)) {
    dirname(paths$lepmap_dir)
  } else {
    file.path(Sys.getenv("COSTUS_ARCHIVE", "/Volumes/dissert/Costus_archive"), "04_lepmap3", "lepmap3")
  }
}
if (!dir.exists(lepmap_wd)) {
  stop("lepmap3 directory not found: ", lepmap_wd, "\nSet LEPMAP_WD or COSTUS_ARCHIVE.")
}
setwd(normalizePath(lepmap_wd, winslash = "/", mustWork = FALSE))
# Clear workspace
rm(list = ls())

# import .call file to correspond marker numbers with reference chrom/pos:
callfile <- "LASIref_minDif0.7_data_0.05covCO_08192025.call"
  # "LASIref_minDif0.7_data_0.05covCO_08192025.call"

# Skip header lines (assume first non-header line has numeric CHR)
call_lines <- readLines(callfile)
data_start <- which(grepl("^\\d", call_lines))[1]

# Read only the lines with marker data
call_data <- read.table(callfile, header = FALSE, skip = data_start - 1)

# Store the CHR column (original reference chromosome) as vector
ref_chr <- call_data$V1  # V1 = CHR

# set color palette
unique_chrs <- sort(unique(ref_chr))
#chr_colors <- setNames(RColorBrewer::brewer.pal(13, "Set3"), unique_chrs)
okabe_ito <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#000000")

# Add extra distinguishable colors
extra_colors <- c("#999999", "#A65628", "#FF69B4", "#00FFFF", "#8A2BE2")

chr_colors <- setNames(
  c(okabe_ito, extra_colors)[seq_along(unique_chrs)],
  unique_chrs
)

# Output PDF file name
output_pdf_name <- "./best_iter_ordermarkers/genetic_map_best_iters_LodLim11_theta0.001_YES_DistLod_allF2s_250923_0.05covCO.pdf"
#output_pdf_name <- "genetic_map_9LGs_LodLim6_theta0.03_NODistLod_250729.pdf"

# Set up output PDF
pdf(output_pdf_name, height = 20, width = 10)
par(mfrow = c(5, 2), mar = c(2, 4, 2, 2))  # set margin sizes

for (i in 1:9) {
  filename <- Sys.glob(
    paste0("./best_iter_ordermarkers/order_LG", i, "_rep*LodLim11_theta0.001_allF2s_0.5covCO.txt")
  )[1]
  lepmap <- read.table(filename, header = FALSE, comment.char = "#")
  
  marker_nums <- lepmap$V1
  cM <- lepmap$V3
  
  # Lookup reference chromosome per marker
  marker_chr <- ref_chr[marker_nums]
  colors <- chr_colors[as.character(marker_chr)]
  
  # Plot
  plot(cM ~ marker_nums,
       main = paste("Linkage Group", i),
       xlab = "Marker number (LepMap)",
       ylab = "cM (LepMap)",
       col = colors,
       pch = 16,
       cex.axis = 1.2,
       cex.lab = 1.2)
  
  # Optional: add legend to only the first plot
  if (i == 1) {
    legend("topleft", legend = names(chr_colors), col = chr_colors,
           pch = 16, cex = 1.5, title = "Ref Chr", ncol = 2)
  }
}

# Close PDF
dev.off()

### old plotting code without color by chromosome:
# # Loop over 9 chromosomes
# for (i in 1:9) {
#   # Construct file name order_chr9_LodLim6_theta0.03
#   filename <- paste0("order_chr", i, "_LodLim5_theta0.03.txt")
#   
#   # Read in data
#   lepmap <- read.table(filename, header = FALSE)
#   
#   # Extract position (V1) and cM (V3)
#   positions <- as.numeric(lepmap$V1)
#   cM <- lepmap$V3
#   
#   # Plot
#   plot(cM ~ positions,
#        main = paste("Linkage Group", i),
#        xlab = "marker number",
#        ylab = "cM (lepmap)",
#        cex.axis = 1.2,
#        cex.lab = 1.2,
#        pch = 16, col = "blue")
# }

######## Lepmap per LG files with chromosome number and position ########

# # import .call file to correspond marker numbers with reference chrom/pos:
# callfile <- "LASIref_minDif0.7_data_0.05covCO_08192025.call"
# 
# # Skip header lines (assume first non-header line has numeric CHR)
# call_lines <- readLines(callfile)
# data_start <- which(grepl("^\\d", call_lines))[1]
# 
# # Read only the lines with marker data
# call_data <- read.table(callfile, header = FALSE, skip = data_start - 1)

# Store the CHR and POS columns (original reference chromosome and bp position)
ref_chr <- call_data$V1  # V1 = CHR
ref_pos <- call_data$V2  # V2 = POS

for (i in 1:9) {
  filenames <- Sys.glob(
    paste0("./best_iter_ordermarkers/order_chr", i, "_rep*LodLim11_theta0.001_allF2s_0.5covCO.txt")
  )
  
  for (filename in filenames) {
    lepmap <- read.table(filename, header = FALSE, comment.char = "#")
    
    marker_nums <- lepmap$V1
    
    # Lookup reference chromosome and position per marker
    marker_chr <- ref_chr[marker_nums]
    marker_pos <- ref_pos[marker_nums]
    
    # Keep only marker number, LepMap position (V2), refCHR, refPOS
    new_df <- data.frame(
      marker = marker_nums,
      lepmap_pos = lepmap$V2,
      refCHR = marker_chr,
      refPOS = marker_pos
    )
    
    # Write new version of file
    outname <- sub("\\.txt$", "_withRef.txt", filename)
    write.table(new_df, outname, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  }
}

#### fix LG4/chr7 ####
# read in data
LG4 <- read.table("./best_iter_ordermarkers/order_LG4_rep93_LodLim11_theta0.001_allF2s_0.5covCO_withRef.txt", header = FALSE)

# remove piece of chr5 (only keep chr 7, V3==7) and bits of 7 that aren't placed correctly (much lower cM than rest of chrom, keep V2>20)
chr7 <- LG4 %>% 
  filter(V3==7, V2>20)

# 4297 markers after all filtering and removals (down from 4838)

# export
write.table(chr7, "./best_iter_ordermarkers/order_chr7_clean_LodLim11_theta0.001_allF2s_0.5covCO_withRef.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# identify region of chrom5 that needs to be masked in the reference since it belongs on 7:
bad_chr5 <- LG4 %>% 
  filter(V3==5)
min(bad_chr5$V4) # 85117780
max(bad_chr5$V4) # 92012694
# 92012694-85117780=6894914, so ~7Mb removed

# identify misplaced markers on chrom7 (really low cM compared to rest)
bad_chr7 <- LG4 %>% 
  filter(V2<20)
min(bad_chr7$V4) # 33052598
max(bad_chr7$V4) # 36425584
# 36425584-33052598=3372986, so ~3Mb removed

#### fix LG6/chr9 ####
# read in data
LG6 <- read.table("./best_iter_ordermarkers/order_LG6_rep56_LodLim11_theta0.001_allF2s_0.5covCO_withRef.txt", header = FALSE)

# remove piece of chr6 that is not actually on chr6 (much higher cM than rest of chrom)
chr9 <- LG6 %>% 
  filter(V2 < 50)

# 2421 markers after filtering

# export
write.table(chr9, "./best_iter_ordermarkers/order_chr9_clean_LodLim11_theta0.001_allF2s_0.5covCO_withRef.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# identify misplaced markers on chrom9 (really high cM compared to rest)
bad_chr9 <- LG6 %>% 
  filter(V2 > 50)
min(bad_chr9$V4) # 75245142
max(bad_chr9$V4) # 78197448
# 78197448-75245142=2952306, so ~3Mb removed

#### marker nums ######
# get number of markers on each chromosome and rename to match cleaned chr 7 and 9:

# DUMB - can just rename and get marker counts from map file:

# LG1/chr5
# 7616 markers

# LG2/chr3
# 7108 markers

# LG3/chr1
# 7021 markers

# LG5/chr8
# 4387 markers

# LG7/chr2
# 1875 markers

# LG8/chr6
# 805 markers

# LG9/chr4
# 4305 markers (after merging LGs 9, 11, and 12, and calling all chr4 singles LG9 too)