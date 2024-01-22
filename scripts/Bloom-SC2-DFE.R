# -----------------------------------------------------------------------------
# Analysis of Bloom and Neher 2023 preprint data for approximate DFE shape
# Author: Chase W. Nelson, chase.nelson@nih.gov
# Cite: Nelson et al.; https://github.com/chasewnelson/Soni_response
# Source: https://academic.oup.com/ve/article/9/2/vead055/7265011?login=false
# Source: https://raw.githubusercontent.com/jbloomlab/SARS2-mut-fitness/main/results/aa_fitness/aamut_fitness_all.csv
# -----------------------------------------------------------------------------


# Import libraries
library(tidyverse)
library(MASS)


# set working directory (top-level)
setwd('Soni_response')  # change to wherever this repository is located on your machine


# -----------------------------------------------------------------------------
# IMPORT mutation data
# Obtain this file from Bloom & Neher, https://raw.githubusercontent.com/jbloomlab/SARS2-mut-fitness/main/results/aa_fitness/aamut_fitness_all.csv
(aamut_fitness_all <- read_csv("data/aamut_fitness_all.csv"))

# if we REMOVE 'clade_founder_aa' column, ensure all rows remain distinct
nrow(aamut_fitness_all)  # 121915
nrow(distinct(dplyr::select(aamut_fitness_all, -clade_founder_aa)))  # 121915 QED


# -----------------------------------------------------------------------------
# CLASSIFY mutation type

# extract WT and MUT aa
(aamut_fitness_all$aa1 <- str_replace(string = aamut_fitness_all$aa_mutation, pattern = "^([\\w*])\\d+([\\w*])$", replacement = "\\1"))
(aamut_fitness_all$aa2 <- str_replace(string = aamut_fitness_all$aa_mutation, pattern = "^([\\w*])\\d+([\\w*])$", replacement = "\\2"))
sort(unique(aamut_fitness_all$aa1))  # "*" "A" "C" "D" "E" "F" "G" "H" "I" "K" "L" "M" "N" "P" "Q" "R" "S" "T" "V" "W" "Y" QED
sort(unique(aamut_fitness_all$aa2))  # "*" "A" "C" "D" "E" "F" "G" "H" "I" "K" "L" "M" "N" "P" "Q" "R" "S" "T" "V" "W" "Y" QED

# verify all lengths 1
all(lengths(aamut_fitness_all$aa1) == 1)  # TRUE QED
all(lengths(aamut_fitness_all$aa2) == 1)  # TRUE QED

# mut_type - order of execution matters
aamut_fitness_all$mut_type <- as.character(NA)
aamut_fitness_all[aamut_fitness_all$aa1 == aamut_fitness_all$aa2, ]$mut_type <- 'syn'
aamut_fitness_all[aamut_fitness_all$aa1 != aamut_fitness_all$aa2, ]$mut_type <- 'nsyn'
aamut_fitness_all[grepl("\\*", aamut_fitness_all$aa1) | grepl("\\*", aamut_fitness_all$aa2), ]$mut_type <- 'stop'


# -----------------------------------------------------------------------------
# FILTER out redundant mutations

# filter to nr = exclude duplicates from ORF1ab; form tibble
(aamut_fitness_all_nr <- as_tibble(filter(aamut_fitness_all, subset_of_ORF1ab == FALSE)))

# INITIALIZE gene list in 5' to 3' order
genes_sorted <- c("ORF1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF9b", "ORF10")

# summarize by mut_type
(aamut_fitness_all_nr %>%
  group_by(mut_type) %>%
  summarise(
   n = n(),
   prop_del = sum(delta_fitness < 0) / n,
   prop_ben = sum(delta_fitness > 0) / n
  ))
# mut_type     n prop_del prop_ben
# 1 nsyn     58324    0.869   0.131 <== too much
# 2 stop      3175    0.985   0.0148
# 3 syn       9557    0.634   0.366 


# -----------------------------------------------------------------------------
# RECREATE & VERIFY Bloom & Neher 2023, Fig. 3A - DFE by type

# factor mut_type
aamut_fitness_all_nr$mut_type <- factor(aamut_fitness_all_nr$mut_type, levels = c('syn', 'nsyn', 'stop'))

# PLOT
(aamut_fitness_all_nr_DFE_PLOT <- ggplot(data = aamut_fitness_all_nr, mapping = aes(x = delta_fitness)) +
  
  geom_histogram(boundary = 0, binwidth = 0.5, mapping = aes(fill = mut_type)) +
  
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey') +
  
  facet_wrap(. ~ mut_type, scales = 'free_y') +
  
  xlab('Fitness effect (ln[O/E])') + ylab('Number of mutations') +
  
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.025))) +
  scale_fill_manual(values = c('#3498DB', '#E74C3C', '#2C3E50')))

# matches, QED


# -----------------------------------------------------------------------------
# DFE: determine % beneficial, etc.
# view(filter(aamut_fitness_all_nr, mut_type == 'syn'))

# (1) NEUTRAL: the middle 95% of synonymous mutations
(NEUTRAL_lo <- as.numeric(quantile(filter(aamut_fitness_all_nr, mut_type == 'syn')$delta_fitness, 0.025)))  # -1.7512 
(NEUTRAL_hi <- as.numeric(quantile(filter(aamut_fitness_all_nr, mut_type == 'syn')$delta_fitness, 0.975)))  # 1.19863 

# (2) LETHAL: anything with median STOP effect or lower
(LETHAL_lo <- min(filter(aamut_fitness_all_nr, mut_type == 'nsyn')$delta_fitness))  # -7.1397
(LETHAL_hi <- median(filter(aamut_fitness_all_nr, mut_type == 'stop')$delta_fitness))  # -3.946 median stop mutation fitness - LETHAL

# (3) DELETERIOUS: greater than lethal but less than 2.5%ile of synonymous
(DEL_lo <- LETHAL_hi)  # -3.946
(DEL_hi <- NEUTRAL_lo)  # -1.7512

# (4) BENEFICIAL: greater than 2.5%ile of synonymous
(BEN_lo <- NEUTRAL_hi)  # 1.19863
(BEN_hi <- max(filter(aamut_fitness_all_nr, mut_type == 'nsyn')$delta_fitness))  # 6.1665

# CLASSIFY NSYN MUTATIONS
# (1) NEUTRAL are <= NEUTRAL_hi, >= NEUTRAL_lo
# (2) LETHAL are <= LETHAL_hi
# (3) DEL are >= LETHAL_hi, <= NETURAL_lo
# (4) BEN are >= NEUTRAL_hi
# CONVENTION is inclusive-exclusive, i.e., [a,b) intervals
aamut_fitness_all_nr$effect_class <- NA
aamut_fitness_all_nr[aamut_fitness_all_nr$delta_fitness <= LETHAL_hi, ]$effect_class <- 'lethal'
aamut_fitness_all_nr[aamut_fitness_all_nr$delta_fitness >= DEL_lo & aamut_fitness_all_nr$delta_fitness < DEL_hi, ]$effect_class <- 'del'
aamut_fitness_all_nr[aamut_fitness_all_nr$delta_fitness >= NEUTRAL_lo & aamut_fitness_all_nr$delta_fitness < NEUTRAL_hi, ]$effect_class <- 'neutral'
aamut_fitness_all_nr[aamut_fitness_all_nr$delta_fitness >= BEN_lo, ]$effect_class <- 'ben'
aamut_fitness_all_nr[aamut_fitness_all_nr$mut_type != 'nsyn', ]$effect_class <- NA

# factor
aamut_fitness_all_nr$effect_class <- factor(aamut_fitness_all_nr$effect_class, levels = c('lethal', 'del', 'neutral', 'ben', NA))


# ------------------------------------------------------------------------------
# WHOLE GENOME SUMMARY
(aamut_fitness_all_nr_NYSN_SUMMARY <- filter(aamut_fitness_all_nr, mut_type == 'nsyn') %>%
      group_by(effect_class) %>%  # mut_type, 
      summarise(
         count = n()
      ))

# add total
aamut_fitness_all_nr_NYSN_SUMMARY$total <- sum(aamut_fitness_all_nr_NYSN_SUMMARY$count)
aamut_fitness_all_nr_NYSN_SUMMARY$total_nonneutral <- sum(filter(aamut_fitness_all_nr_NYSN_SUMMARY, effect_class != 'neutral')$count)

# proportion of all mutations, nonneutral mutations
aamut_fitness_all_nr_NYSN_SUMMARY$prop <- aamut_fitness_all_nr_NYSN_SUMMARY$count / aamut_fitness_all_nr_NYSN_SUMMARY$total
aamut_fitness_all_nr_NYSN_SUMMARY$prop_nonneutral <- aamut_fitness_all_nr_NYSN_SUMMARY$count / aamut_fitness_all_nr_NYSN_SUMMARY$total_nonneutral
aamut_fitness_all_nr_NYSN_SUMMARY[aamut_fitness_all_nr_NYSN_SUMMARY$effect_class == 'neutral', ]$prop_nonneutral <- NA

# ADD 95% CI
aamut_fitness_all_nr_NYSN_SUMMARY$prop_95CI_lo <- as.numeric(NA)
aamut_fitness_all_nr_NYSN_SUMMARY$prop_95CI_hi <- as.numeric(NA)
aamut_fitness_all_nr_NYSN_SUMMARY$prop_nonneutral_95CI_lo <- as.numeric(NA)
aamut_fitness_all_nr_NYSN_SUMMARY$prop_nonneutral_95CI_hi <- as.numeric(NA)

for (i in 1:nrow(aamut_fitness_all_nr_NYSN_SUMMARY)) {
   # i <- 1
   prop_test_result <- prop.test(x = aamut_fitness_all_nr_NYSN_SUMMARY[i, ]$count, n = aamut_fitness_all_nr_NYSN_SUMMARY[i, ]$total, conf.level = 0.95)
   aamut_fitness_all_nr_NYSN_SUMMARY[i, ]$prop_95CI_lo <- prop_test_result$conf.int[1]
   aamut_fitness_all_nr_NYSN_SUMMARY[i, ]$prop_95CI_hi <- prop_test_result$conf.int[2]
   
   if(aamut_fitness_all_nr_NYSN_SUMMARY[i, ]$effect_class != 'neutral') {
      prop_test_nonneutral_result <- prop.test(x = aamut_fitness_all_nr_NYSN_SUMMARY[i, ]$count, n = aamut_fitness_all_nr_NYSN_SUMMARY[i, ]$total_nonneutral, conf.level = 0.95)
      aamut_fitness_all_nr_NYSN_SUMMARY[i, ]$prop_nonneutral_95CI_lo <- prop_test_nonneutral_result$conf.int[1]
      aamut_fitness_all_nr_NYSN_SUMMARY[i, ]$prop_nonneutral_95CI_hi <- prop_test_nonneutral_result$conf.int[2]
   }
}

# view
aamut_fitness_all_nr_NYSN_SUMMARY
# 1.49% of all mutations (2.84% of nonneutral mutations) beneficial <== more reasonable


# ----------
# ORF-BY-ORF SUMMARY
(aamut_fitness_all_nr_NYSN_SUMMARY_byORF <- filter(aamut_fitness_all_nr, mut_type == 'nsyn') %>%
    group_by(gene, effect_class) %>%
    summarise(
       count = n()
    ))  
# N.B.: some classes are missing for some ORFs

# add missing classes
for (this_effect in c('lethal', 'del', 'neutral', 'ben')) {
   for (this_gene in genes_sorted) {
      if(nrow(filter(aamut_fitness_all_nr_NYSN_SUMMARY_byORF, gene == this_gene, effect_class == this_effect)) == 0) {
         cat(paste0("gene=", this_gene, " | effect=", this_effect, '\n'))
         
         # add
         aamut_fitness_all_nr_NYSN_SUMMARY_byORF <- rbind(aamut_fitness_all_nr_NYSN_SUMMARY_byORF,
                                                          tibble(gene = this_gene, effect_class = this_effect, count = 0))
      }
   }
}

# ORF-BY-ORF SUMS
(aamut_fitness_all_nr_NYSN_ORFsums <- filter(aamut_fitness_all_nr, mut_type == 'nsyn') %>%
      group_by(gene) %>%
      summarise(
         total = n()
      ))

# ORF-BY-ORF SUMS NONNEUTRAL
(aamut_fitness_all_nr_NYSN_ORFsumsNonneutral <- filter(aamut_fitness_all_nr, mut_type == 'nsyn', effect_class != 'neutral') %>%
      group_by(gene) %>%
      summarise(
         total_nonneutral = n()
      ))

# JOIN SUMS
aamut_fitness_all_nr_NYSN_SUMMARY_byORF <- left_join(aamut_fitness_all_nr_NYSN_SUMMARY_byORF, aamut_fitness_all_nr_NYSN_ORFsums, by = 'gene')
aamut_fitness_all_nr_NYSN_SUMMARY_byORF <- left_join(aamut_fitness_all_nr_NYSN_SUMMARY_byORF, aamut_fitness_all_nr_NYSN_ORFsumsNonneutral, by = 'gene')

# proportion of all mutations, nonneutral mutations
aamut_fitness_all_nr_NYSN_SUMMARY_byORF$prop <- aamut_fitness_all_nr_NYSN_SUMMARY_byORF$count / aamut_fitness_all_nr_NYSN_SUMMARY_byORF$total
aamut_fitness_all_nr_NYSN_SUMMARY_byORF$prop_nonneutral <- 
   aamut_fitness_all_nr_NYSN_SUMMARY_byORF$count / aamut_fitness_all_nr_NYSN_SUMMARY_byORF$total_nonneutral
aamut_fitness_all_nr_NYSN_SUMMARY_byORF[aamut_fitness_all_nr_NYSN_SUMMARY_byORF$effect_class == 'neutral', ]$prop_nonneutral <- NA

# ADD 95% CI
aamut_fitness_all_nr_NYSN_SUMMARY_byORF$prop_95CI_lo <- as.numeric(NA)
aamut_fitness_all_nr_NYSN_SUMMARY_byORF$prop_95CI_hi <- as.numeric(NA)
aamut_fitness_all_nr_NYSN_SUMMARY_byORF$prop_nonneutral_95CI_lo <- as.numeric(NA)
aamut_fitness_all_nr_NYSN_SUMMARY_byORF$prop_nonneutral_95CI_hi <- as.numeric(NA)

for (i in 1:nrow(aamut_fitness_all_nr_NYSN_SUMMARY_byORF)) {
   # i <- 11
   # cat(i, ' ')
   prop_test_result <- prop.test(x = aamut_fitness_all_nr_NYSN_SUMMARY_byORF[i, ]$count, n = aamut_fitness_all_nr_NYSN_SUMMARY_byORF[i, ]$total, conf.level = 0.95)
   aamut_fitness_all_nr_NYSN_SUMMARY_byORF[i, ]$prop_95CI_lo <- prop_test_result$conf.int[1]
   aamut_fitness_all_nr_NYSN_SUMMARY_byORF[i, ]$prop_95CI_hi <- prop_test_result$conf.int[2]
   
   if(aamut_fitness_all_nr_NYSN_SUMMARY_byORF[i, ]$effect_class != 'neutral') {
      prop_test_nonneutral_result <- prop.test(x = aamut_fitness_all_nr_NYSN_SUMMARY_byORF[i, ]$count, n = aamut_fitness_all_nr_NYSN_SUMMARY_byORF[i, ]$total_nonneutral, conf.level = 0.95)
      aamut_fitness_all_nr_NYSN_SUMMARY_byORF[i, ]$prop_nonneutral_95CI_lo <- prop_test_nonneutral_result$conf.int[1]
      aamut_fitness_all_nr_NYSN_SUMMARY_byORF[i, ]$prop_nonneutral_95CI_hi <- prop_test_nonneutral_result$conf.int[2]
   }
}


# ----- 
# COMBINE

# add column and rearrange
aamut_fitness_all_nr_NYSN_SUMMARY$gene <- 'Genome'
aamut_fitness_all_nr_NYSN_SUMMARY <- dplyr::select(aamut_fitness_all_nr_NYSN_SUMMARY, gene, everything())

# combine
aamut_fitness_all_nr_NYSN_SUMMARY <- rbind(aamut_fitness_all_nr_NYSN_SUMMARY, aamut_fitness_all_nr_NYSN_SUMMARY_byORF)

# factor gene
aamut_fitness_all_nr_NYSN_SUMMARY$gene <- factor(aamut_fitness_all_nr_NYSN_SUMMARY$gene,
                                                      levels = c(genes_sorted, 'Genome'))
aamut_fitness_all_nr_NYSN_SUMMARY <- arrange(aamut_fitness_all_nr_NYSN_SUMMARY, gene)  # sorts by factor levels

# SAVE (already in repository for you) or RELOAD
# write_tsv(aamut_fitness_all_nr_NYSN_SUMMARY, "data/aamut_fitness_all_nr_NYSN_SUMMARY.tsv")
# aamut_fitness_all_nr_NYSN_SUMMARY <- read_tsv("data/aamut_fitness_all_nr_NYSN_SUMMARY.tsv")

# -----
# TO WIDE

# pivot WIDE
aamut_fitness_all_nr_NYSN_SUMMARY_WIDE <- aamut_fitness_all_nr_NYSN_SUMMARY %>%
   pivot_wider(names_from = effect_class, values_from = c("count", "prop", "prop_nonneutral", "prop_95CI_lo", "prop_95CI_hi", "prop_nonneutral_95CI_lo", "prop_nonneutral_95CI_hi"))

# factor gene
aamut_fitness_all_nr_NYSN_SUMMARY_WIDE$gene <- factor(aamut_fitness_all_nr_NYSN_SUMMARY_WIDE$gene,
                                                      levels = c(genes_sorted, 'Genome'))
aamut_fitness_all_nr_NYSN_SUMMARY_WIDE <- arrange(aamut_fitness_all_nr_NYSN_SUMMARY_WIDE, gene)  # sorts by factor levels

# sort columns
names(aamut_fitness_all_nr_NYSN_SUMMARY_WIDE)
aamut_fitness_all_nr_NYSN_SUMMARY_WIDE <- dplyr::select(aamut_fitness_all_nr_NYSN_SUMMARY_WIDE,
                                                        c("gene", "total", "total_nonneutral", "count_lethal", "count_del", "count_neutral", "count_ben", "prop_lethal", "prop_95CI_lo_lethal", "prop_95CI_hi_lethal", "prop_del", "prop_95CI_lo_del", "prop_95CI_hi_del", "prop_neutral", "prop_95CI_lo_neutral", "prop_95CI_hi_neutral", "prop_ben", "prop_95CI_lo_ben", "prop_95CI_hi_ben", "prop_nonneutral_ben", "prop_nonneutral_95CI_lo_ben", "prop_nonneutral_95CI_hi_ben", "prop_nonneutral_lethal", "prop_nonneutral_95CI_lo_lethal", "prop_nonneutral_95CI_hi_lethal", "prop_nonneutral_del", "prop_nonneutral_95CI_lo_del", "prop_nonneutral_95CI_hi_del", "prop_nonneutral_neutral", "prop_nonneutral_95CI_lo_neutral", "prop_nonneutral_95CI_hi_neutral"))


# -----------------------------------------------------------------------------
# SLIDING WINDOWS - calculate these values in sliding windows
# Gu et al. used a sliding window size of 30 codons, step size of 1 codon, following our work in San et al.
(aamut_fitness_all_nr <- dplyr::select(aamut_fitness_all_nr, gene, aa_site, aa_mutation, everything()))
(aamut_fitness_all_nr <- arrange(aamut_fitness_all_nr, gene, aa_site, aa_mutation))

WINDOW_SIZE <- 30  # codons (amino acids)

# make sure our manually ordered gene list is exhaustive
all(genes_sorted %in% unique(aamut_fitness_all_nr$gene)) && all(unique(aamut_fitness_all_nr$gene) %in% genes_sorted)  # TRUE

# -----
# INITIALIZE
window_results <- tibble(gene = character(),
                         start = integer(),
                         end = integer(),
                         total = integer(), total_nonneutral = integer(),
                         ben_min_effect = numeric(), del_min_effect = numeric(), lethal_min_effect = numeric(), neutral_min_effect = numeric(),
                         ben_median_effect = numeric(), del_median_effect = numeric(), lethal_median_effect = numeric(), neutral_median_effect = numeric(),
                         ben_mean_effect = numeric(), del_mean_effect = numeric(), lethal_mean_effect = numeric(), neutral_mean_effect = numeric(),
                         ben_max_effect = numeric(), del_max_effect = numeric(), lethal_max_effect = numeric(), neutral_max_effect = numeric(),
                         ben_sd_effect = numeric(), del_sd_effect = numeric(), lethal_sd_effect = numeric(), neutral_sd_effect = numeric(),
                         ben_count = integer(), del_count = integer(), lethal_count = integer(), neutral_count = integer(),
                         ben_prop = numeric(), del_prop = numeric(), lethal_prop = numeric(), neutral_prop = numeric(),
                         
                         ben_prop_95CI_lo = numeric(), del_prop_95CI_lo = numeric(), lethal_prop_95CI_lo = numeric(), neutral_prop_95CI_lo = numeric(),
                         ben_prop_95CI_hi = numeric(), del_prop_95CI_hi = numeric(), lethal_prop_95CI_hi = numeric(), neutral_prop_95CI_hi = numeric(),
                         
                         ben_prop_nonneutral = numeric(), del_prop_nonneutral = numeric(), lethal_prop_nonneutral = numeric(), neutral_prop_nonneutral = numeric(),
                         
                         ben_prop_nonneutral_95CI_lo = numeric(), del_prop_nonneutral_95CI_lo = numeric(), lethal_prop_nonneutral_95CI_lo = numeric(), neutral_prop_nonneutral_95CI_lo = numeric(),
                         ben_prop_nonneutral_95CI_hi = numeric(), del_prop_nonneutral_95CI_hi = numeric(), lethal_prop_nonneutral_95CI_hi = numeric(), neutral_prop_nonneutral_95CI_hi = numeric())

# SLIDING WINDOW LOOP FOR FRACTIONS lethal/del/neutral/ben <== takes < 2 minutes
for (this_gene in genes_sorted) {
   # this_gene <- 'S'
   this_gene_data <- filter(aamut_fitness_all_nr, gene == this_gene)
   this_gene_MAX_SITE <- max(this_gene_data$aa_site)  # 1272
   
   # LOOP WINDOWS
   for (this_START in 1:(this_gene_MAX_SITE - WINDOW_SIZE + 1)) {
      # this_START <- 1
      this_END <- this_START + WINDOW_SIZE - 1
      
      # extract this window
      this_window_data <- filter(this_gene_data, aa_site >= this_START, aa_site <= this_END)
      
      # summarize window characteristics
      this_window_nsyn_results <- filter(this_window_data, mut_type == 'nsyn') %>%
         group_by(effect_class) %>%
         summarise(
            min_effect = min(delta_fitness),
            median_effect = median(delta_fitness),
            mean_effect = mean(delta_fitness),
            max_effect = max(delta_fitness),
            sd_effect = sd(delta_fitness),
            count = n()
         )
      
      # add missing classes
      for (this_effect in c('lethal', 'del', 'neutral', 'ben')) {
         if(nrow(filter(this_window_nsyn_results, effect_class == this_effect)) == 0) {
            
            # add
            this_window_nsyn_results <- rbind(this_window_nsyn_results,
                                              tibble(effect_class = this_effect, 
                                                     min_effect = NA, 
                                                     median_effect = NA,
                                                     mean_effect = NA,
                                                     max_effect = NA,
                                                     sd_effect = NA,
                                                     count = 0))
         }
      }
      
      # add totals
      this_TOTAL <- sum(this_window_nsyn_results$count)
      this_TOTAL_NONNEUTRAL <- sum(filter(this_window_nsyn_results, effect_class != 'neutral')$count)
      this_window_nsyn_results$total <- this_TOTAL
      this_window_nsyn_results$total_nonneutral <- this_TOTAL_NONNEUTRAL
      
      # props
      this_window_nsyn_results$prop <- this_window_nsyn_results$count / this_window_nsyn_results$total
      this_window_nsyn_results$prop_nonneutral <- this_window_nsyn_results$count / this_window_nsyn_results$total_nonneutral
      this_window_nsyn_results[this_window_nsyn_results$effect_class == 'neutral', ]$prop_nonneutral <- NA
      
      # add 95% CI
      this_window_nsyn_results$prop_95CI_lo <- as.numeric(NA)
      this_window_nsyn_results$prop_95CI_hi <- as.numeric(NA)
      this_window_nsyn_results$prop_nonneutral_95CI_lo <- as.numeric(NA)
      this_window_nsyn_results$prop_nonneutral_95CI_hi <- as.numeric(NA)
      
      for (i in 1:nrow(this_window_nsyn_results)) {
         # i <- 4
         prop_test_result <- prop.test(x = this_window_nsyn_results[i, ]$count, n = this_window_nsyn_results[i, ]$total, conf.level = 0.95)
         this_window_nsyn_results[i, ]$prop_95CI_lo <- prop_test_result$conf.int[1]
         this_window_nsyn_results[i, ]$prop_95CI_hi <- prop_test_result$conf.int[2]
         
         if(this_window_nsyn_results[i, ]$effect_class != 'neutral') {
            prop_test_nonneutral_result <- prop.test(x = this_window_nsyn_results[i, ]$count, n = this_window_nsyn_results[i, ]$total_nonneutral, conf.level = 0.95)
            this_window_nsyn_results[i, ]$prop_nonneutral_95CI_lo <- prop_test_nonneutral_result$conf.int[1]
            this_window_nsyn_results[i, ]$prop_nonneutral_95CI_hi <- prop_test_nonneutral_result$conf.int[2]
         }
      }
      
      # ---------------------------------------------------------------------
      # Append to results & increment
      
      this_row <- tibble(gene = as.character(this_gene), start = as.integer(this_START), end = as.integer(this_END),
                         total = as.integer(this_TOTAL), total_nonneutral = as.integer(this_TOTAL_NONNEUTRAL),
                         
                         ben_min_effect = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'ben',]$min_effect),
                         del_min_effect = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'del',]$min_effect),
                         lethal_min_effect = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'lethal',]$min_effect),
                         neutral_min_effect = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'neutral',]$min_effect),
                         
                         ben_median_effect = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'ben',]$median_effect),
                         del_median_effect = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'del',]$median_effect),
                         lethal_median_effect = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'lethal',]$median_effect),
                         neutral_median_effect = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'neutral',]$median_effect),
                         
                         ben_mean_effect = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'ben',]$mean_effect),
                         del_mean_effect = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'del',]$mean_effect),
                         lethal_mean_effect = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'lethal',]$mean_effect),
                         neutral_mean_effect = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'neutral',]$mean_effect),
                         
                         ben_max_effect = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'ben',]$max_effect),
                         del_max_effect = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'del',]$max_effect),
                         lethal_max_effect = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'lethal',]$max_effect),
                         neutral_max_effect = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'neutral',]$max_effect),
                         
                         ben_sd_effect = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'ben',]$sd_effect),
                         del_sd_effect = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'del',]$sd_effect),
                         lethal_sd_effect = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'lethal',]$sd_effect),
                         neutral_sd_effect = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'neutral',]$sd_effect),
                         
                         ben_count = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'ben',]$count),
                         del_count = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'del',]$count),
                         lethal_count = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'lethal',]$count),
                         neutral_count = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'neutral',]$count),
                         
                         ben_prop = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'ben',]$prop),
                         del_prop = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'del',]$prop),
                         lethal_prop = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'lethal',]$prop),
                         neutral_prop = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'neutral',]$prop),
                         
                         ben_prop_95CI_lo = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'ben',]$prop_95CI_lo),
                         del_prop_95CI_lo = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'del',]$prop_95CI_lo),
                         lethal_prop_95CI_lo = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'lethal',]$prop_95CI_lo),
                         neutral_prop_95CI_lo = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'neutral',]$prop_95CI_lo),
                         
                         ben_prop_95CI_hi = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'ben',]$prop_95CI_hi),
                         del_prop_95CI_hi = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'del',]$prop_95CI_hi),
                         lethal_prop_95CI_hi = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'lethal',]$prop_95CI_hi),
                         neutral_prop_95CI_hi = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'neutral',]$prop_95CI_hi),
                         
                         ben_prop_nonneutral = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'ben',]$prop_nonneutral),
                         del_prop_nonneutral = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'del',]$prop_nonneutral),
                         lethal_prop_nonneutral = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'lethal',]$prop_nonneutral),
                         neutral_prop_nonneutral = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'neutral',]$prop_nonneutral),
                         
                         ben_prop_nonneutral_95CI_lo = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'ben',]$prop_nonneutral_95CI_lo),
                         del_prop_nonneutral_95CI_lo = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'del',]$prop_nonneutral_95CI_lo),
                         lethal_prop_nonneutral_95CI_lo = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'lethal',]$prop_nonneutral_95CI_lo),
                         neutral_prop_nonneutral_95CI_lo = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'neutral',]$prop_nonneutral_95CI_lo),
                         
                         ben_prop_nonneutral_95CI_hi = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'ben',]$prop_nonneutral_95CI_hi),
                         del_prop_nonneutral_95CI_hi = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'del',]$prop_nonneutral_95CI_hi),
                         lethal_prop_nonneutral_95CI_hi = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'lethal',]$prop_nonneutral_95CI_hi),
                         neutral_prop_nonneutral_95CI_hi = as.numeric(this_window_nsyn_results[this_window_nsyn_results$effect_class == 'neutral',]$prop_nonneutral_95CI_hi))
      
      window_results <- rbind(window_results, this_row)
   }
}

# select columns for order
window_results <- dplyr::select(window_results,
                                c('gene', 'start', 'end', 'total', 'total_nonneutral', 'ben_count', 'del_count', 'lethal_count', 'neutral_count', 'ben_prop', 'ben_prop_95CI_lo', 'ben_prop_95CI_hi', 'del_prop', 'del_prop_95CI_lo', 'del_prop_95CI_hi', 'lethal_prop', 'lethal_prop_95CI_lo', 'lethal_prop_95CI_hi', 'neutral_prop', 'neutral_prop_95CI_lo', 'neutral_prop_95CI_hi', 'ben_min_effect', 'del_min_effect', 'lethal_min_effect', 'neutral_min_effect', 'ben_median_effect', 'del_median_effect', 'lethal_median_effect', 'neutral_median_effect', 'ben_mean_effect', 'del_mean_effect', 'lethal_mean_effect', 'neutral_mean_effect', 'ben_max_effect', 'del_max_effect', 'lethal_max_effect', 'neutral_max_effect', 'ben_sd_effect', 'del_sd_effect', 'lethal_sd_effect', 'neutral_sd_effect', 'ben_prop_nonneutral', 'del_prop_nonneutral', 'lethal_prop_nonneutral', 'neutral_prop_nonneutral', 'ben_prop_nonneutral_95CI_lo', 'del_prop_nonneutral_95CI_lo', 'lethal_prop_nonneutral_95CI_lo', 'neutral_prop_nonneutral_95CI_lo', 'ben_prop_nonneutral_95CI_hi', 'del_prop_nonneutral_95CI_hi', 'lethal_prop_nonneutral_95CI_hi', 'neutral_prop_nonneutral_95CI_hi'))

# SAVE or RELOAD
# write_tsv(window_results, "data/BloomFB_window_results.tsv")  # this already exists in repository for you
# window_results <- read_tsv("data/BloomFB_window_results.tsv")


# -----------------------------------------------------------------------------
# Scale the fitness effects to estimate selection coefficients

# selection coefficient
aamut_fitness_all_nr_NYSN <- filter(aamut_fitness_all_nr, mut_type == 'nsyn')  # 58,324 × 16
aamut_fitness_all_nr_NYSN$sel_coeff <- aamut_fitness_all_nr_NYSN$delta_fitness / max(abs(aamut_fitness_all_nr_NYSN$delta_fitness))
summary(aamut_fitness_all_nr_NYSN$sel_coeff)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -1.00000 -0.43140 -0.25304 -0.26371 -0.08557  0.86369 


# ------------------------------------------------------------------------------
# FIT WHOLE GENOME distributions

# ----------
# DEL: FIT GAMMA to del for ALL s<0
del_fit_gamma <- fitdistr(abs(filter(aamut_fitness_all_nr_NYSN, sel_coeff < 0)$sel_coeff), "gamma")
# shape         rate    
# 1.702711424   5.375853847 
# (0.009829622) (0.036030126)

# mean = shape * scale
mean(filter(aamut_fitness_all_nr_NYSN, sel_coeff < 0)$sel_coeff)  # -0.3167334
1.702711424 * (1 / 5.375853847)  # 0.3167332 MEAN QED

# check fit
hist(abs(filter(aamut_fitness_all_nr_NYSN, sel_coeff < 0)$sel_coeff), probability = TRUE, breaks = 30)
curve(dgamma(x, del_fit_gamma$estimate[1], del_fit_gamma$estimate[2]), col = "blue", add = TRUE, lwd = 2)

# ----------
# BEN: FIT EXP TO BEN
(lambda <- 1/mean(filter(aamut_fitness_all_nr_NYSN, sel_coeff > 0)$sel_coeff))
ben_fit_exponential <- fitdistr(filter(aamut_fitness_all_nr_NYSN, sel_coeff > 0)$sel_coeff, "exponential")
# rate   
# 11.5319164 
# ( 0.1317181)

hist(filter(aamut_fitness_all_nr_NYSN, sel_coeff > 0)$sel_coeff, probability = TRUE, breaks = 30)
curve(dexp(x, rate = lambda), col = "blue", add = TRUE, lwd = 2)

mean(filter(aamut_fitness_all_nr_NYSN, sel_coeff > 0)$sel_coeff)  # 0.08671586
1 / 11.5319164  # 0.08671586 QED

