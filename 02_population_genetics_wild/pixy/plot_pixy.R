### Code to plot pixy output: pi, Dxy, and Fst
# Julia Harencar
# 03/16/2024

# load packages
library(tidyverse)
library(gridExtra)

# function from the pixy website (https://pixy.readthedocs.io/en/latest/plotting.html)
pixy_to_long <- function(pixy_files){
  pixy_df <- list()
  for(i in 1:length(pixy_files)){
    stat_file_type <- gsub(".*_|.txt", "", pixy_files[i])
    if(stat_file_type == "pi"){
      df <- read_delim(pixy_files[i], delim = "\t")
      df <- df %>%
        gather(-pop, -window_pos_1, -window_pos_2, -chromosome,
               key = "statistic", value = "value") %>%
        rename(pop1 = pop) %>%
        mutate(pop2 = NA)
      
      pixy_df[[i]] <- df
    } else{
      df <- read_delim(pixy_files[i], delim = "\t")
      df <- df %>%
        gather(-pop1, -pop2, -window_pos_1, -window_pos_2, -chromosome,
               key = "statistic", value = "value")
      pixy_df[[i]] <- df
    }
  }
  bind_rows(pixy_df) %>%
    arrange(pop1, pop2, chromosome, window_pos_1, statistic)
}

# Import my data
for (p in c(
  file.path("scripts", "project_paths.R"),
  file.path("..", "..", "scripts", "project_paths.R")
)) {
  if (file.exists(p)) { source(p); break }
}
paths <- project_paths()
pixy_folder <- paths$pixy_folder
if (is.na(pixy_folder) || !nzchar(pixy_folder) || !dir.exists(pixy_folder)) {
  stop("pixy output not found: ", pixy_folder, "\nSet PIXY_FOLDER to the Hudson 10 kb run.")
}
pixy_files <- list.files(pixy_folder, pattern = "_(pi|dxy|fst)\\.txt$", full.names = TRUE)
if (length(pixy_files) == 0) {
  stop("No pixy *_pi.txt / *_dxy.txt / *_fst.txt files in ", pixy_folder)
}
pixy_df <- pixy_to_long(pixy_files)
# Fst and Dxy use the same 10 kb floors as the genome-wide figures. Pi is unchanged.
window_keys <- c("pop1", "pop2", "chromosome", "window_pos_1", "window_pos_2")
fst_keep <- pixy_df %>%
  filter(statistic == "no_snps", suppressWarnings(as.numeric(value)) >= 10) %>%
  distinct(across(all_of(window_keys)))
dxy_keep <- pixy_df %>%
  filter(statistic == "no_sites", suppressWarnings(as.numeric(value)) >= 1000) %>%
  distinct(across(all_of(window_keys)))
pixy_df <- bind_rows(
  pixy_df %>% filter(statistic == "avg_pi"),
  pixy_df %>% filter(statistic == "avg_dxy") %>% inner_join(dxy_keep, by = window_keys),
  pixy_df %>% filter(statistic == "avg_hudson_fst") %>% inner_join(fst_keep, by = window_keys)
)
pixy_stats <- unique(pixy_df$statistic)
if (!"avg_hudson_fst" %in% pixy_stats) {
  if ("avg_wc_fst" %in% pixy_stats) {
    stop("Pixy folder is Weir–Cockerham (avg_wc_fst). Use the Hudson run, not output_251124_LASIref.")
  }
  stop("Pixy folder is missing avg_hudson_fst: ", pixy_folder)
}

# custom labeller for special characters in pi/dxy/fst
pixy_labeller <- as_labeller(c(avg_pi = "pi",
                               avg_dxy = "D[XY]",
                               avg_hudson_fst = "F[ST]"),
                             default = label_parsed)

# plotting summary statistics along a single chromosome
plot_list = list()
for (i in 1:9) {
  chr_num <- paste0("Chrom", i)
p <- pixy_df %>%
  filter(chromosome == chr_num) %>% 
  filter(statistic %in% c("avg_pi", "avg_dxy", "avg_hudson_fst")) %>%
  mutate(chr_position = ((window_pos_1 + window_pos_2)/2)/1000000) %>%
  ggplot(aes(x = chr_position, y = value, color = statistic))+
  geom_line(size = 0.25)+
  facet_grid(statistic ~ .,
             scales = "free_y", switch = "x", space = "free_x",
             labeller = labeller(statistic = pixy_labeller,
                                 value = label_value))+
  xlab(paste("Position on Chromosome", i,"(Mb)")) +
  ylab("")+
  theme_bw()+
  theme(panel.spacing = unit(0.1, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")+
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  scale_color_brewer(palette = "Set1")
plot_list[[i]] = p
}

plot_list[[1]] + plot_list[[2]] + plot_list[[3]] + plot_list[[4]] + 
  plot_list[[5]] + plot_list[[6]] + plot_list[[7]] + plot_list[[8]] + 
  plot_list[[9]] + plot_layout(ncol = 3)

# not working:
#ggpubr::ggarrange(plotlist = plot_list, nrow = 3, ncol = 3, byrow = TRUE)


pixy_df %>%
  mutate(chrom_color_group = case_when(as.numeric(chromosome) %% 2 != 0 ~ "even",
                                       TRUE ~ "odd" )) %>%
  mutate(chromosome = factor(chromosome, levels = c(1:22, "X", "Y"))) %>%
  filter(statistic %in% c("avg_pi", "avg_dxy", "avg_hudson_fst")) %>%
  ggplot(aes(x = (window_pos_1 + window_pos_2)/2, y = value, color = chrom_color_group))+
  geom_point(size = 0.5, alpha = 0.5, stroke = 0)+
  facet_grid(statistic ~ chromosome,
             scales = "free_y", switch = "x", space = "free_x",
             labeller = labeller(statistic = pixy_labeller,
                                 value = label_value))+
  xlab("Chromsome")+
  ylab("Statistic Value")+
  scale_color_manual(values = c("grey50", "black"))+
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position ="none")+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,NA))
