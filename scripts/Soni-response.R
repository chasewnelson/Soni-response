# --------------------------------------------------------------------------------
# Analyse Soni et al. SLiM results
# Author: Chase W. Nelson, chase.nelson@nih.gov
# Cite: Nelson et al.; https://github.com/chasewnelson/Soni_response
# Description: wrangle and analyse results of Soni et al. SARS-CoV-2 intrahost evolution simulations
# UPDATED for Soni et al. revision
# --------------------------------------------------------------------------------

# Import libraries
library(tidyverse)
library(scales)
library(RColorBrewer)
library(boot)
library(zoo)

# set working directory (top-level)
setwd('Soni_response')  # change to wherever this repository is located on your machine


# --------------------------------------------------------------------------------
# Figure 1a - Recreate Soni et al. nsyn DFEs

# DFE 1 - Weakly deleterious background
# 10% f0, which is m1 (s = 0)
# 70% f1, which is m2 (s = -0.01 to -0.001)
# 10% f2, which is m3 (s = -0.1 to -0.01)
# 10% f3, which is m4 (s = -1 to -0.1)

# DFE 2 - Strongly deleterious background
# 10% f0, which is m1 (s = 0)
# 10% f1, which is m2 (s = -0.01 to -0.001)
# 10% f2, which is m3 (s = -0.1 to -0.01)
# 70% f3, which is m4 (s = -1 to -0.1)

# DFE 3 - Empirical with 1% beneficial - Soni et al. added in their v2
# 53.7% f0, which is m1 (s = 0) | .542/(.542+.112+.02+.326+.01)==0.5366337
# 11.1% f1, which is m2 (s = -0.01 to -0.001) | .112/(.542+.112+.02+.326+.01)==0.1108911
# 2.0% f2, which is m3 (s = -0.1 to -0.01) | .02/(.542+.112+.02+.326+.01)==0.01980198
# 32.3% f3, which is m4 (s = -1 to -0.1) | .326/(.542+.112+.02+.326+.01)==0.3227723
# 1.0% EXTRA fb, which is m0 (s = 0.05 to 0.13) | .01/(.542+.112+.02+.326+.01)==0.00990099

# DFE 4 - Empirical with 10% beneficial - Soni et al. added in their v2
# 44.5% f0, which is m1 (s = 0) | 0.445/(0.445+0.112+0.02+0.326+0.097)==0.445
# 11.2% f1, which is m2 (s = -0.01 to -0.001) | 0.112/(0.445+0.112+0.02+0.326+0.097)==0.112
# 2.0% f2, which is m3 (s = -0.1 to -0.01) | 0.02/(0.445+0.112+0.02+0.326+0.097)==0.02
# 32.6% f3, which is m4 (s = -1 to -0.1) | 0.326/(0.445+0.112+0.02+0.326+0.097)==0.326
# 9.7% EXTRA fb, which is m0 (s = 0.05 to 0.13) | 0.097/(0.445+0.112+0.02+0.326+0.097)==0.097

# DFE 5 - Bloom & Neher-based - we added in our updated response
# 0.09068 lethal, "m1", -1.0
# 0.42043 deleterious, "m2", gamma mean = -0.31673, shape = 1.70271
# 0.47396 neutral, "m3", 0.0
# 0.01493 beneficial, "m4", exponential mean = 0.08672

# Simulate 10k mutations using each of their DFEs, with Flynn beneficials at 1% and 10%, BloomDFE whole genome values
DFE_points_10000 <- tibble(DFE = c(rep('DFE1', 10000),
                                   rep('DFE2', 10000),
                                   rep('DFE3', 10000),
                                   rep('DFE4', 10000),
                                   rep('DFE5', 10000)),
                           s = c(c(runif(1000, -1, -0.1), runif(1000, -0.1, -0.01), runif(7000, -0.01, -0.001), rep(0, 1000)),
                                 c(runif(7000, -1, -0.1), runif(1000, -0.1, -0.01), runif(1000, -0.01, -0.001), rep(0, 1000)),
                                 c(runif(99, 0.05, 0.13), runif(3228, -1, -0.1), runif(198, -0.1, -0.01), runif(1109, -0.01, -0.001), rep(0, 5366)),
                                 c(runif(970, 0.05, 0.13), runif(3260, -1, -0.1), runif(200, -0.1, -0.01), runif(1120, -0.01, -0.001), rep(0, 4450)),
                                 c(rep(-1, 907), -rgamma(n = 4204, shape = 1.702711424, rate = 5.375853847), rep(0, 4740), rexp(149, rate = 1 / 0.08671586))
                           ))

# replace s < -1 with simply -1 (it's effectively lethal)
DFE_points_10000[DFE_points_10000$s < -1, ]$s <- -1

# PLOT
(Soni_DFE_VIOLINPLOT <- ggplot(filter(DFE_points_10000), aes(x = as.factor(1), y = s)) +
  
      # VIOLIN
      geom_violin(linewidth = 0.2, fill = 'lightgrey') +
      
      # s=0
      geom_hline(yintercept = 0, linetype = 'solid', color = 'black', size = .25) +
      
      facet_wrap(. ~ DFE, ncol = 1) +
      coord_flip() +
      theme_bw() +
      scale_x_discrete(expand = expansion(mult = c(0.5, 0.5))) +
      scale_y_continuous(limits = c(-1, 1), breaks = pretty_breaks(), expand = expansion(mult = c(0, 0))))

# SAVE SOURCE
# write_tsv(DFE_points_10000, 'Figure_DFEs_source.txt')


# -----------------------------------------------------------------------------
# IMPORT simulation results for analyses shown in Figures 1b, 1c, 2b

# Mutations
(DFE_mutations <- read_delim("data/DFE_mutations.txt",
                             col_names = c('file_tempID', 'ID', 'mut_type', 'pos', 's', 'h', 'subpop_ID', 'tick', 'freq'), 
                             delim = " "))
(DFE_beneficial_mutations <- read_delim("data/DFE_beneficial_mutations.txt", 
                                        col_names = c('file_tempID', 'ID', 'mut_type', 'pos', 's', 'h', 'subpop_ID', 'tick', 'freq'), 
                                        delim = " "))

# Added in Soni et al. updated response
(Flynn1pct_mutations <- read_delim("data/Flynn1pct_mutations.txt",
                             col_names = c('file_tempID', 'ID', 'mut_type', 'pos', 's', 'h', 'subpop_ID', 'tick', 'freq'), 
                             delim = " "))
(Flynn10pct_mutations <- read_delim("data/Flynn10pct_mutations.txt",
                                        col_names = c('file_tempID', 'ID', 'mut_type', 'pos', 's', 'h', 'subpop_ID', 'tick', 'freq'), 
                                        delim = " "))

# Added in response to reviewers - new estimated Bloom & Neher DFE
(BloomDFE_mutations <- read_delim("data/BloomDFE_mutations.txt",
                                   col_names = c('file_tempID', 'ID', 'mut_type', 'pos', 's', 'h', 'subpop_ID', 'tick', 'freq'), 
                                   delim = " "))

# --------------------------------------------------------------------------------
# MUTATION ANALYSIS

# Combine and wrangle
mutations <- rbind(DFE_mutations, DFE_beneficial_mutations, Flynn1pct_mutations, Flynn10pct_mutations, BloomDFE_mutations)

# extract metadata from first column
mutations$tempID <- as.integer(str_replace(mutations$file_tempID, ".*:", ""))
mutations$sample_size <- as.integer(str_replace(mutations$file_tempID, ".*_(\\d+)\\.out.*", "\\1"))
mutations$model_type <- str_replace(mutations$file_tempID, "^results\\/(\\w+)\\/.*", "\\1")
mutations$model_subtype <- str_replace(mutations$file_tempID, "^results.*\\/(\\w+)_rep.*", "\\1")
mutations$rep <- as.integer(str_replace(mutations$file_tempID, ".*_rep(\\d+)_.*", "\\1"))

# new model column
mutations$model <- as.character(NA)
mutations[mutations$model_type == 'DFE', ]$model <- paste0('Neu/', mutations[mutations$model_type == 'DFE', ]$model_subtype)
mutations[mutations$model_type == 'DFE_beneficial', ]$model <- paste0('Sel/', mutations[mutations$model_type == 'DFE_beneficial', ]$model_subtype)
mutations[mutations$model_type == 'Flynn1pct', ]$model <- paste0('Emp/', mutations[mutations$model_type == 'Flynn1pct', ]$model_subtype)
mutations[mutations$model_type == 'Flynn10pct', ]$model <- paste0('Emp/', mutations[mutations$model_type == 'Flynn10pct', ]$model_subtype)
mutations[mutations$model_type == 'BloomDFE', ]$model <- paste0('Emp/', mutations[mutations$model_type == 'BloomDFE', ]$model_subtype)
unique(mutations$model)  # "Neu/DFE1"       "Neu/DFE2"       "Sel/DFE1"       "Sel/DFE2"       "Emp/Flynn1pct"  "Emp/Flynn10pct" "Emp/BloomDFE"  

# new AF, MAF columns
mutations$AF <- mutations$freq / mutations$sample_size
mutations$MAF <- mutations$AF
mutations[mutations$AF > 0.5, ]$MAF <- 1 - mutations[mutations$AF > 0.5, ]$AF

# reorder
mutations <- dplyr::select(mutations, -file_tempID)
mutations <- dplyr::select(mutations, model, model_type, model_subtype, sample_size, rep, tempID, everything())

# verify
unique(mutations$model)  # "Neu/DFE1"       "Neu/DFE2"       "Sel/DFE1"       "Sel/DFE2"       "Emp/Flynn1pct"  "Emp/Flynn10pct" "Emp/BloomDFE"  
unique(mutations$model_type)  # "DFE"            "DFE_beneficial" "Flynn1pct"      "Flynn10pct"     "BloomDFE" 
unique(mutations$model_subtype)  # "DFE1"       "DFE2"       "Flynn1pct"  "Flynn10pct" "BloomDFE"  
unique(mutations$sample_size)  # 100 1000
length(unique(mutations$rep))  # 100 QED
unique(mutations$mut_type)  # "m1" "m2" "m3" "m4" "m0" 
# All good!

# -----
# Estimate fb proportion of beneficial mutations in DFE_beneficial_mutations
(mutations_DFE_ben_max_IDs <- filter(mutations, model_type == 'DFE_beneficial') %>%
   group_by(model_subtype, rep) %>%
   summarise(
      max_tempID = max(tempID),
      max_ID = max(ID)
   ))

1 / mean(mutations_DFE_ben_max_IDs$max_ID)  # 7.15365e-07 ==> 0.00007%
# mean == 1397888
# 7.153649e-07 the mean fraction of beneficial mutations here, i.e., 7.153649e-07 = 0.0000007 = 0.00007%

# ------
# add site
mutations$site <- mutations$pos + 1

# add site type - they model 2/3 NSYN, 1/3 SYN
mutations$site_type <- as.character(NA)
mutations[mutations$site %% 3 == 1, ]$site_type <- 'N'
mutations[mutations$site %% 3 == 2, ]$site_type <- 'N'
mutations[mutations$site %% 3 == 0, ]$site_type <- 'S'

# add codon number
mutations$codon <- as.integer(NA)
mutations$codon <- ceiling(mutations$site / 3)

# reorder
mutations <- dplyr::select(mutations, -file_tempID)  # redundant
mutations <- dplyr::select(mutations, model, model_type, model_subtype, sample_size, rep, tempID, ID, pos, site, everything())

# Filter to Gu specifications of MAF >= 0.025
mutations_filtered <- filter(mutations, MAF >= 0.025)

# multiallelic sites?
mutations_per_site <- mutations_filtered %>%
  group_by(model, sample_size, rep, site) %>%
  summarise(
   num_muts = n()
  )

nrow(mutations_filtered)  # 32774
nrow(filter(mutations_per_site, num_muts > 1))  # just 23 across all models/reps, can safely ignore multiallelic sites

# Mutations per replicate (analogous to per sample)
(mutations_per_rep <- mutations_filtered %>%
 group_by(model, sample_size, rep) %>%
 summarise(
  mut_count = n(),
  mean_MAF = mean(MAF),
  median_MAF = median(MAF),
  sd_MAF = sd(MAF),
  Q1_MAF = quantile(MAF, 0.25),
  Q3_MAF = quantile(MAF, 0.75)))

# ADD ROWS for reps that have NO (0) mutations left
for (this_model in unique(mutations_per_rep$model)) {
 # this_model <- 'Neu/DFE1'
 
 for (this_sample_size in unique (mutations_per_rep$sample_size)) {
  # this_sample_size <- 100
  
  for (this_rep in 1:100) {
   # this_rep <- 4
   
   if (nrow(filter(mutations_per_rep, model == this_model, sample_size == this_sample_size, rep == this_rep)) == 0) {
    mutations_per_rep <- rbind(mutations_per_rep,
                               tibble(model = this_model, sample_size = this_sample_size, rep = this_rep, mut_count = 0,
                                      mean_MAF = NA, median_MAF = NA, sd_MAF = NA, Q1_MAF = NA, Q3_MAF = NA))
   }
  }
 }
}


# --------------------------------------------------------------------------------
# NUCLEOTIDE DIVERISTY (PI) for each site

# pi numerator
mutations_filtered$num_pw_diffs <- (mutations_filtered$freq) * (mutations_filtered$sample_size - mutations_filtered$freq)

# pi denominator
mutations_filtered$num_pw_comps <- choose(mutations_filtered$sample_size, 2)

# pi
mutations_filtered$pi <- mutations_filtered$num_pw_diffs / mutations_filtered$num_pw_comps  # this is only at sites with variants


# --------------------------------------------------------------------------------
# SONI METHOD for piN/piS

# overall pi values for each rep by site type
(mutations_pi_byRep <- mutations_filtered %>%
 group_by(model, sample_size, rep, site_type) %>%
 summarise(
  mut_site_count = n(),
  pi_numerator = sum(pi)))

# ADD ROWS for reps that have NO (0) mutations; their pi is 0 - this takes < 60sec
for (this_model in unique(mutations_pi_byRep$model)) {
 # this_model <- 'Neu/DFE1'
 
 for (this_sample_size in unique (mutations_pi_byRep$sample_size)) {
  # this_sample_size <- 100
  
  for (this_rep in 1:100) {
   # this_rep <- 1
   
   for (this_site_type in unique(mutations_pi_byRep$site_type)) {
    # this_site_type <- 'S'
    
    if (nrow(filter(mutations_pi_byRep, model == this_model, sample_size == this_sample_size, rep == this_rep, site_type == this_site_type)) == 0) {
     mutations_pi_byRep <- rbind(mutations_pi_byRep,
                                 tibble(model = this_model, sample_size = this_sample_size, rep = this_rep, site_type = this_site_type, 
                                        mut_site_count = 0, pi_numerator = 0))
    }
   }
  }
 }
}

# SORT columns
mutations_pi_byRep <- arrange(mutations_pi_byRep, model, sample_size, rep, site_type)

# numbers of each site type
mutations_pi_byRep$site_count <- as.integer(NA)
mutations_pi_byRep[mutations_pi_byRep$site_type == 'N', ]$site_count <- 30000 * (2/3)  # according to Soni et al. model; ~3/4 in reality
mutations_pi_byRep[mutations_pi_byRep$site_type == 'S', ]$site_count <- 30000 * (1/3)  # according to Soni et al. model; ~1/4 in reality
mutations_pi_byRep$rep_pi <- mutations_pi_byRep$pi_numerator / mutations_pi_byRep$site_count

# take the means across those reps, i.e. reps (samples) as the unit
(mutations_pi_byRep_means <- mutations_pi_byRep %>%
      group_by(model, sample_size, site_type) %>%
      summarise(
         n = n(),
         mean_pi = mean(rep_pi),
         sd_pi = sd(rep_pi),
         SE_pi = sd_pi / sqrt(n),
         Q1_pi = quantile(rep_pi, 0.25),
         Q3_pi = quantile(rep_pi, 0.75)))


# --------------------------------------------------------------------------------
# GU METHOD for piN/piS

# INTIALIZE template
codon_based_pi <- tibble(model = c(rep("Neu/DFE1", 2000000), 
                                   rep("Neu/DFE2", 2000000), 
                                   rep("Sel/DFE1", 2000000), 
                                   rep("Sel/DFE2", 2000000),
                                   rep("Emp/Flynn1pct", 2000000),
                                   rep("Emp/Flynn10pct", 2000000),
                                   rep("Emp/BloomDFE", 2000000)
                                   ))

# add cols
codon_based_pi$sample_size <- rep(c(rep(100, 1000000), rep(1000, 1000000)), 7)
codon_based_pi$rep <- rep(1:100, 140000)

# sort columns
codon_based_pi <- arrange(codon_based_pi, model, sample_size, rep)

# add codon col
codon_based_pi$codon <- rep(1:10000, 1400)

# sort columns
codon_based_pi <- arrange(codon_based_pi, model, sample_size, rep, codon)

# JOIN PI SUMS
mutations_filtered_joiner <- dplyr::select(mutations_filtered, model, sample_size, rep, codon, site_type, pi)
mutations_filtered_joiner <- arrange(mutations_filtered_joiner, model, sample_size, rep, codon, site_type)

# join N_diffs
codon_based_pi <- left_join(codon_based_pi, 
                            dplyr::select(filter(mutations_filtered_joiner, site_type == 'N'), -site_type), 
                            by = c('model', 'sample_size', 'rep', 'codon'))
names(codon_based_pi)[names(codon_based_pi) == 'pi'] <- 'N_diffs'

# join S_diffs
codon_based_pi <- left_join(codon_based_pi, 
                            dplyr::select(filter(mutations_filtered_joiner, site_type == 'S'), -site_type), 
                            by = c('model', 'sample_size', 'rep', 'codon'))
names(codon_based_pi)[names(codon_based_pi) == 'pi'] <- 'S_diffs'

# set NA values to 0, they were absent
codon_based_pi[is.na(codon_based_pi$N_diffs), ]$N_diffs <- 0
codon_based_pi[is.na(codon_based_pi$S_diffs), ]$S_diffs <- 0

# initialize N and S sites, all the same in the model
codon_based_pi$N_sites <- 2
codon_based_pi$S_sites <- 1

# calculate means for each model/sample_size/rep
codon_based_pi_repMeans <- codon_based_pi %>%
 group_by(model, sample_size, rep) %>%
 summarise(
  N_diffs = sum(N_diffs),
  S_diffs = sum(S_diffs),
  N_sites = sum(N_sites),
  S_sites = sum(S_sites)
 )

# add piN, piS
codon_based_pi_repMeans$piN <- codon_based_pi_repMeans$N_diffs / codon_based_pi_repMeans$N_sites
codon_based_pi_repMeans$piS <- codon_based_pi_repMeans$S_diffs / codon_based_pi_repMeans$S_sites

# make long
codon_based_pi_repMeans_LONG <- pivot_longer(codon_based_pi_repMeans,
                                              cols = c('piN', 'piS'),
                                              names_to = "site_type", values_to = "pi")

# rep means
# Soni et al. 68% 'confidence interval' is a standard deviation
1 - (1-.95)/2  # 0.975
1 - (1-.68)/2  # 0.84

codon_based_pi_repMeans_LONG_SUMMARY <- codon_based_pi_repMeans_LONG %>%
  group_by(model, sample_size, site_type) %>%
  summarise(
   num_reps = n(),
   mean_pi = mean(pi),
   sd_pi = sd(pi),
   SE_pi = sd_pi / sqrt(num_reps),
   q95lo = quantile(pi, 0.025),
   Q1_pi = quantile(pi, 0.25),
   Q3_pi = quantile(pi, 0.75),
   q95hi = quantile(pi, 0.975),
   CI68_pi = SE_pi * qt(0.84, num_reps - 1),
   CI95_pi = SE_pi * qt(0.975, num_reps - 1))

# Wilcoxon Signed Rank tests using sample (replicate) as the unit
wilcox.test(filter(codon_based_pi_repMeans, model == 'Neu/DFE1', sample_size == 100)$piN, 
            filter(codon_based_pi_repMeans, model == 'Neu/DFE1', sample_size == 100)$piS, paired = TRUE)
# 1.567e-06

wilcox.test(filter(codon_based_pi_repMeans, model == 'Neu/DFE2', sample_size == 100)$piN, 
            filter(codon_based_pi_repMeans, model == 'Neu/DFE2', sample_size == 100)$piS, paired = TRUE)
# 3.26e-07

wilcox.test(filter(codon_based_pi_repMeans, model == 'Sel/DFE1', sample_size == 100)$piN, 
            filter(codon_based_pi_repMeans, model == 'Sel/DFE1', sample_size == 100)$piS, paired = TRUE)
# 8.298e-10

wilcox.test(filter(codon_based_pi_repMeans, model == 'Sel/DFE2', sample_size == 100)$piN, 
            filter(codon_based_pi_repMeans, model == 'Sel/DFE2', sample_size == 100)$piS, paired = TRUE)
# < 2.2e-16

wilcox.test(filter(codon_based_pi_repMeans, model == 'Emp/Flynn1pct', sample_size == 100)$piN, 
            filter(codon_based_pi_repMeans, model == 'Emp/Flynn1pct', sample_size == 100)$piS, paired = TRUE)
# 0.006706

wilcox.test(filter(codon_based_pi_repMeans, model == 'Emp/Flynn10pct', sample_size == 100)$piN, 
            filter(codon_based_pi_repMeans, model == 'Emp/Flynn10pct', sample_size == 100)$piS, paired = TRUE)
# < 2.2e-16

wilcox.test(filter(codon_based_pi_repMeans, model == 'Emp/BloomDFE', sample_size == 100)$piN, 
            filter(codon_based_pi_repMeans, model == 'Emp/BloomDFE', sample_size == 100)$piS, paired = TRUE)
# 0.01178

# POINT PLOT
pi_corr_factor = 1e5
(codon_based_pi_repMeans_LONG_SUMMARY_POINTPLOT <- ggplot(filter(codon_based_pi_repMeans_LONG_SUMMARY, sample_size == 100),  # ! model %in% c('Emp/Flynn1pct', 'Emp/Flynn10pct')
                                                          mapping = aes(x = as.factor(1), y = mean_pi * pi_corr_factor, color = site_type)) +
  geom_point(position = position_dodge(.2), size = 0.75) +
  
  # SD - which Soni et al. use
  # geom_errorbar(mapping = aes(ymin = mean_pi - sd_pi, ymax = mean_pi + sd_pi), linewidth = .35, width = 0, position = position_dodge(0.8)) +
  
  # 68% CI
  # geom_errorbar(mapping = aes(ymin = mean_pi - CI68_pi, ymax = mean_pi + CI68_pi), linewidth = .35, width = 0, position = position_dodge(0.8)) +
  
  # SE, the standard error - conventional
  geom_errorbar(mapping = aes(ymin = (mean_pi - SE_pi) * pi_corr_factor, ymax = (mean_pi + SE_pi) * pi_corr_factor), linewidth = .35, width = 0,
                position = position_dodge(.2)) +
  
  # IQR
  # geom_errorbar(mapping = aes(ymin = Q1_pi, ymax = Q3_pi), linewidth = .35, width = 0, position = position_dodge(0.8)) +
  
  # empirical 95% CI
  # geom_errorbar(mapping = aes(ymin = q95lo, ymax = q95hi), linewidth = .35, width = 0, position = position_dodge(0.8)) +
  
  # facet_wrap(. ~ factor(model,
  #                       levels = c("Neu/DFE1", "Sel/DFE1", "Neu/DFE2", "Sel/DFE2"),
  #                       labels = c("Weak del", "Weak del + ben", "Strong del", "Strong del + ben")),
  #            nrow = 2) +
      
      facet_wrap(. ~ model, nrow = 2) +
  
  ylab('pi x 10^5') +
  theme_bw() +
  scale_color_manual(values = c('#E74C3C', '#3498DB')) +
  scale_x_discrete(expand = expansion(mult = c(0.2, 0.2))) +
  scale_y_continuous(limits = c(0, NA), breaks = pretty_breaks(), expand = expansion(mult = c(0, 0.15))))

# POINT PLOT - not shown, ~identical to our method


# --------------------------------------------------------------------------------
# BOOTSTRAP PROCESS on codon_based_pi as Gu et al. did it

# calculate means for each model/sample_size/codon across all reps (equivalent of biological samples)
codon_based_pi_codonMeans <- codon_based_pi %>%
  group_by(model, sample_size, codon) %>%
  summarise(
   N_diffs = mean(N_diffs),
   S_diffs = mean(S_diffs),
   N_sites = mean(N_sites),
   S_sites = mean(S_sites)
  )

# MANUALLY ENSURE THIS IS THE CASE BEFOREHAND
codon_based_pi_codonMeans$num_defined_seqs <- codon_based_pi_codonMeans$sample_size


# --------------------------------------------------------------------------------
# *BASIC* BOOTSTRAP FUNCTION (dN - dS) for CODON UNIT, modified from https://github.com/chasewnelson/SNPGenie sliding windows script
dNdS_diff_boot_fun <- function(codon_results, numerator, denominator, num_replicates, num_cpus) {
 
 # Function for dN
 dN_function <- function(D, indices) {
  dN <- sum(D[indices, paste0(numerator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(numerator, "_sites")], na.rm = TRUE)
  return(dN)
 }
 
 # Function for dN
 dS_function <- function(D, indices) {
  dS <- sum(D[indices, paste0(denominator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(denominator, "_sites")], na.rm = TRUE)
  return(dS)
 }
 
 # Function for dN - dS
 dN_m_dS_function <- function(D, indices) {
  dN <- sum(D[indices, paste0(numerator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(numerator, "_sites")], na.rm = TRUE)
  dS <- sum(D[indices, paste0(denominator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(denominator, "_sites")], na.rm = TRUE)
  dN_m_dS <- dN - dS
  return(dN_m_dS)
 }
 
 # Function for dN/dS
 dN_over_dS_function <- function(D, indices) {
  dN <- sum(D[indices, paste0(numerator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(numerator, "_sites")], na.rm = TRUE)
  dS <- sum(D[indices, paste0(denominator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(denominator, "_sites")], na.rm = TRUE)
  dN_over_dS <- dN / dS
  return(dN_over_dS)
 }
 
 # CREATE FUNCTION FOR dN/dS TO CALCULATE ITS SE
 
 (dN <- sum(codon_results[ , paste0(numerator, "_diffs")], na.rm = TRUE) / sum(codon_results[ , paste0(numerator, "_sites")], na.rm = TRUE))
 (dS <- sum(codon_results[ , paste0(denominator, "_diffs")], na.rm = TRUE) / sum(codon_results[ , paste0(denominator, "_sites")], na.rm = TRUE))
 (dNdS <- dN / dS)
 
 # Run the BOOTSTRAPS
 # boot dN
 (boot_dN <- boot(data = codon_results, R = num_replicates, statistic = dN_function, parallel = 'multicore', ncpus = num_cpus))
 (dN <- boot_dN$t0)
 (boot_dN_SE <- sd(boot_dN$t))
 
 # boot dS
 (boot_dS <- boot(data = codon_results, R = num_replicates, statistic = dS_function, parallel = 'multicore', ncpus = num_cpus))
 (dS <- boot_dS$t0)
 (boot_dS_SE <- sd(boot_dS$t))
 
 # boot dN - dS
 (boot_dN_m_dS <- boot(data = codon_results, R = num_replicates, statistic = dN_m_dS_function, parallel = 'multicore', ncpus = num_cpus))
 (dN_m_dS <- boot_dN_m_dS$t0)
 (boot_dN_m_dS_SE <- sd(boot_dN_m_dS$t))
 (boot_dN_m_dS_Z <- dN_m_dS / boot_dN_m_dS_SE)
 (boot_dN_m_dS_P <- 2 * pnorm(-abs(boot_dN_m_dS_Z)))
 
 # boot dN/dS
 (boot_dN_over_dS <- boot(data = codon_results, R = num_replicates, statistic = dN_over_dS_function, parallel = 'multicore', ncpus = num_cpus))
 (dN_over_dS <- boot_dN_over_dS$t0)
 (boot_dN_over_dS_SE <- sd(boot_dN_over_dS$t))
 (boot_dN_over_dS_Z <- dN_over_dS / boot_dN_over_dS_SE)
 (boot_dN_over_dS_P <- 2 * pnorm(-abs(boot_dN_over_dS_Z)))
 
 ### NEW: ASL (acheived significance level)
 boot_dN_gt_dS_count <- sum(boot_dN_m_dS$t > 0) # 345
 boot_dN_eq_dS_count <- sum(boot_dN_m_dS$t == 0) # 0
 boot_dN_lt_dS_count <- sum(boot_dN_m_dS$t < 0) # 655
 ASL_dN_gt_dS_P <- boot_dN_lt_dS_count / (boot_dN_gt_dS_count + boot_dN_eq_dS_count + boot_dN_lt_dS_count)
 ASL_dN_lt_dS_P <- boot_dN_gt_dS_count / (boot_dN_gt_dS_count + boot_dN_eq_dS_count + boot_dN_lt_dS_count)
 
 return(paste(num_replicates, dN, dS, dNdS, dN_m_dS, boot_dN_SE, boot_dS_SE, boot_dN_over_dS_SE, boot_dN_over_dS_P, 
              boot_dN_m_dS_SE, boot_dN_m_dS_P, 
              boot_dN_gt_dS_count, boot_dN_eq_dS_count, boot_dN_lt_dS_count, ASL_dN_gt_dS_P, ASL_dN_lt_dS_P,
              sep = "\t"))
}


# --------------------------------------------------------------------------------
# ANALYSIS VARIABLES
MIN_DEFINED_CODONS <- 6
NBOOTSTRAPS <- 1000
NCPUS <- 10


# --------------------------------------------------------------------------------
# INITIALIZE DATA FRAME: modified for intrahost
codon_based_pi_results_bootstrap <- data.frame(
 model = character(),
 sample_size = integer(),
 num_bootstraps = integer(),
 min_defined_codons = integer(),
 num_codons = integer(),
 N_sites = numeric(),
 S_sites = numeric(),
 N_diffs = numeric(),
 S_diffs = numeric(),
 num_replicates = integer(),
 dN = numeric(),
 dS = numeric(),
 dNdS = numeric(),
 dN_m_dS = numeric(),
 boot_dN_SE = numeric(),
 boot_dS_SE = numeric(),
 boot_dN_over_dS_SE = numeric(),
 boot_dN_over_dS_P = numeric(),
 boot_dN_m_dS_SE = numeric(),
 P_value = numeric(),
 boot_dN_gt_dS_count = integer(), 
 boot_dN_eq_dS_count = integer(), 
 boot_dN_lt_dS_count = integer(), 
 ASL_dN_gt_dS_P = numeric(), 
 ASL_dN_lt_dS_P = numeric())

# LOOP EACH DATA SUBSET - takes < 60sec
for (this_model in sort(unique(codon_based_pi_codonMeans$model))) {
 # this_model <- 'Neu/DFE1'
 
 for (this_sample_size in sort(unique(codon_based_pi_codonMeans$sample_size))) {
  # this_sample_size <- 100
  
  # Filter by gene, frame, and minimum number of defined codons
  this_data <- filter(codon_based_pi_codonMeans, model == this_model, sample_size == this_sample_size, num_defined_seqs >= MIN_DEFINED_CODONS) 
  
  if(nrow(this_data) >= 2) {
   # LEADING SUMMARY COLUMNS:
   N_sites <- sum(this_data$N_sites, na.rm = T)
   S_sites <- sum(this_data$S_sites, na.rm = T)
   N_diffs <- sum(this_data$N_diffs, na.rm = T)
   S_diffs <- sum(this_data$S_diffs, na.rm = T)
   
   summary_data <- paste(nrow(this_data),
                         N_sites, S_sites, 
                         N_diffs, S_diffs, 
                         sep = "\t")
   
   # BOOTSTRAP THE ONE RATIO
   boot_dNdS <- dNdS_diff_boot_fun(this_data, 'N', 'S', NBOOTSTRAPS, NCPUS)
   
   # RECORD HEADER
   boot_vector_names <- c('num_codons', 'N_sites', 'S_sites', 'N_diffs', 'S_diffs',
                          'num_replicates', 
                          'dN', 'dS', 'dNdS', 'dN_m_dS', 'boot_dN_SE', 'boot_dS_SE', 'boot_dN_over_dS_SE', 'boot_dN_over_dS_P', 'boot_dN_m_dS_SE', 'P_value',
                          'boot_dN_gt_dS_count', 'boot_dN_eq_dS_count', 'boot_dN_lt_dS_count', 'ASL_dN_gt_dS_P', 'ASL_dN_lt_dS_P')
   
   boot_dNdS_vector <- unlist(str_split(string = paste(summary_data, boot_dNdS, sep = "\t"), pattern = "\t"))
   
   # Add names
   names(boot_dNdS_vector) <- boot_vector_names
   
   # Prepare additional rows
   codon_based_pi_results_bootstrap_ADDITION <- data.frame(model = this_model,
                                                           sample_size = this_sample_size,
                                                           num_bootstraps = NBOOTSTRAPS,
                                                           min_defined_codons = MIN_DEFINED_CODONS,
                                                           num_codons = as.integer(boot_dNdS_vector['num_codons']),
                                                           N_sites = as.numeric(boot_dNdS_vector['N_sites']),
                                                           S_sites = as.numeric(boot_dNdS_vector['S_sites']),
                                                           N_diffs = as.numeric(boot_dNdS_vector['N_diffs']),
                                                           S_diffs = as.numeric(boot_dNdS_vector['S_diffs']),
                                                           num_replicates = as.integer(boot_dNdS_vector['num_replicates']),
                                                           dN = as.numeric(boot_dNdS_vector['dN']),
                                                           dS = as.numeric(boot_dNdS_vector['dS']),
                                                           dNdS = as.numeric(boot_dNdS_vector['dNdS']),
                                                           dN_m_dS = as.numeric(boot_dNdS_vector['dN_m_dS']),
                                                           boot_dN_SE = as.numeric(boot_dNdS_vector['boot_dN_SE']),
                                                           boot_dS_SE = as.numeric(boot_dNdS_vector['boot_dS_SE']),
                                                           boot_dN_over_dS_SE = as.numeric(boot_dNdS_vector['boot_dN_over_dS_SE']),
                                                           boot_dN_over_dS_P = as.numeric(boot_dNdS_vector['boot_dN_over_dS_P']),
                                                           boot_dN_m_dS_SE = as.numeric(boot_dNdS_vector['boot_dN_m_dS_SE']),
                                                           P_value = as.numeric(boot_dNdS_vector['P_value']),
                                                           boot_dN_gt_dS_count = as.numeric(boot_dNdS_vector['boot_dN_gt_dS_count']),
                                                           boot_dN_eq_dS_count = as.numeric(boot_dNdS_vector['boot_dN_eq_dS_count']),
                                                           boot_dN_lt_dS_count = as.numeric(boot_dNdS_vector['boot_dN_lt_dS_count']),
                                                           ASL_dN_gt_dS_P = as.numeric(boot_dNdS_vector['ASL_dN_gt_dS_P']),
                                                           ASL_dN_lt_dS_P = as.numeric(boot_dNdS_vector['ASL_dN_lt_dS_P']))
   
   # Add the new rows to results
   codon_based_pi_results_bootstrap <- rbind(codon_based_pi_results_bootstrap, codon_based_pi_results_bootstrap_ADDITION)
  }
 }
}

# label columns
names(codon_based_pi_results_bootstrap) <- c('model', 'sample_size',
                                             'num_bootstraps', 'min_defined_codons', 'num_codons', 
                                             'N_sites', 'S_sites', 'N_diffs', 'S_diffs', 
                                             'num_replicates', 'dN', 'dS', 'dNdS', 'dN_m_dS', 
                                             'boot_dN_SE', 'boot_dS_SE', 'boot_dN_over_dS_SE', 'boot_dN_over_dS_P', 'boot_dN_m_dS_SE', 'P_value',
                                             'boot_dN_gt_dS_count', 'boot_dN_eq_dS_count', 'boot_dN_lt_dS_count', 'ASL_dN_gt_dS_P', 'ASL_dN_lt_dS_P')

# Manual 2-sided ASL P-value
codon_based_pi_results_bootstrap$P_ALS <- NA
codon_based_pi_results_bootstrap[! is.na(codon_based_pi_results_bootstrap$ASL_dN_gt_dS_P) & codon_based_pi_results_bootstrap$ASL_dN_gt_dS_P < codon_based_pi_results_bootstrap$ASL_dN_lt_dS_P, ]$P_ALS <- 
 2 * codon_based_pi_results_bootstrap[! is.na(codon_based_pi_results_bootstrap$ASL_dN_gt_dS_P) & codon_based_pi_results_bootstrap$ASL_dN_gt_dS_P < codon_based_pi_results_bootstrap$ASL_dN_lt_dS_P, ]$ASL_dN_gt_dS_P
codon_based_pi_results_bootstrap[! is.na(codon_based_pi_results_bootstrap$ASL_dN_gt_dS_P) & codon_based_pi_results_bootstrap$ASL_dN_gt_dS_P > codon_based_pi_results_bootstrap$ASL_dN_lt_dS_P, ]$P_ALS <- 
 2 * codon_based_pi_results_bootstrap[! is.na(codon_based_pi_results_bootstrap$ASL_dN_gt_dS_P) & codon_based_pi_results_bootstrap$ASL_dN_gt_dS_P > codon_based_pi_results_bootstrap$ASL_dN_lt_dS_P, ]$ASL_dN_lt_dS_P
codon_based_pi_results_bootstrap[! is.na(codon_based_pi_results_bootstrap$P_ALS) & codon_based_pi_results_bootstrap$P_ALS == 0, ]$P_ALS <- 1 / NBOOTSTRAPS

# FDR
codon_based_pi_results_bootstrap$Q_ASL_BH <- p.adjust(codon_based_pi_results_bootstrap$P_ALS, method = "BH")
codon_based_pi_results_bootstrap$Q_Z_BH <- p.adjust(codon_based_pi_results_bootstrap$P_value, method = "BH")

# PLOT again to see it matches the sample-based method
(codon_based_pi_results_bootstrap_LONG <- pivot_longer(codon_based_pi_results_bootstrap,
                                                       cols = c('dN', 'dS'),
                                                       names_to = "site_type",
                                                       values_to = "pi"))

codon_based_pi_results_bootstrap_LONG$pi_SE <- codon_based_pi_results_bootstrap_LONG$boot_dN_SE
codon_based_pi_results_bootstrap_LONG[codon_based_pi_results_bootstrap_LONG$site_type == 'dS', ]$pi_SE <- 
 codon_based_pi_results_bootstrap_LONG[codon_based_pi_results_bootstrap_LONG$site_type == 'dS', ]$boot_dS_SE

# SAVE SOURCE
# write_tsv(codon_based_pi_results_bootstrap_LONG, 'Figure2b_source.txt')

# RELOAD codon_based_pi_results_bootstrap_LONG from Figure2b_source.txt here
# codon_based_pi_results_bootstrap_LONG <- read_tsv("data/Figure2b_source.txt")

# PLOT
pi_corr_factor <- 1e5
(codon_based_pi_repMeans_LONG_SUMMARY_POINTPLOT2 <- ggplot(filter(codon_based_pi_results_bootstrap_LONG, sample_size == 100), 
                                                           mapping = aes(x = as.factor(1), y = pi * pi_corr_factor, color = site_type)) +
  geom_point(position = position_dodge(.2), size = 0.75) +
  
  # SD
  # geom_errorbar(mapping = aes(ymin = mean_pi - sd_pi, ymax = mean_pi + sd_pi), linewidth = .35, width = 0, position = position_dodge(0.8)) +
  
  # 68% CI
  # geom_errorbar(mapping = aes(ymin = mean_pi - CI68_pi, ymax = mean_pi + CI68_pi), linewidth = .35, width = 0, position = position_dodge(0.8)) +
  
  # SE, the standard error
  geom_errorbar(mapping = aes(ymin = (pi - pi_SE) * pi_corr_factor, ymax = (pi + pi_SE) * pi_corr_factor), linewidth = .35, width = 0, 
                position = position_dodge(.2)) +
  
  # IQR
  # geom_errorbar(mapping = aes(ymin = Q1_pi, ymax = Q3_pi), linewidth = .35, width = 0, position = position_dodge(0.8)) +
  
  # empirical 95% CI
  # geom_errorbar(mapping = aes(ymin = q95lo, ymax = q95hi), linewidth = .35, width = 0, position = position_dodge(0.8)) +
  
      facet_wrap(. ~ model, nrow = 2) +
  
  ylab('pi x 10^5') +
  theme_bw() +
  scale_color_manual(values = c('#E74C3C', '#3498DB')) +
  scale_x_discrete(expand = expansion(mult = c(0.2, 0.2))) +
  scale_y_continuous(limits = c(0, NA), breaks = pretty_breaks(), expand = expansion(mult = c(0, 0.15))))

# PLOT Figure 2 (revision)
pi_corr_factor2 <- 1e4
(codon_based_pi_repMeans_LONG_SUMMARY_POINTPLOT3 <- ggplot(filter(codon_based_pi_results_bootstrap_LONG, sample_size == 100,
                                                                  model %in% c('Emp/Flynn1pct', 'Emp/Flynn10pct', 'Emp/BloomDFE')), 
                                                           mapping = aes(x = as.factor(1), y = pi * pi_corr_factor2, color = site_type)) +
      geom_point(position = position_dodge(.2), size = 0.75) +
      
      # SD
      # geom_errorbar(mapping = aes(ymin = mean_pi - sd_pi, ymax = mean_pi + sd_pi), linewidth = .35, width = 0, position = position_dodge(0.8)) +
      
      # 68% CI
      # geom_errorbar(mapping = aes(ymin = mean_pi - CI68_pi, ymax = mean_pi + CI68_pi), linewidth = .35, width = 0, position = position_dodge(0.8)) +
      
      # SE, the standard error
      geom_errorbar(mapping = aes(ymin = (pi - pi_SE) * pi_corr_factor2, ymax = (pi + pi_SE) * pi_corr_factor2), linewidth = .35, width = 0, 
                    position = position_dodge(.2)) +
      
      # IQR
      # geom_errorbar(mapping = aes(ymin = Q1_pi, ymax = Q3_pi), linewidth = .35, width = 0, position = position_dodge(0.8)) +
      
      # empirical 95% CI
      # geom_errorbar(mapping = aes(ymin = q95lo, ymax = q95hi), linewidth = .35, width = 0, position = position_dodge(0.8)) +
      
      facet_wrap(. ~ model, ncol = 1) +
      
      ylab('pi x 10^4') +
      theme_bw() +
      scale_color_manual(values = c('#E74C3C', '#3498DB')) +
      scale_x_discrete(expand = expansion(mult = c(0.2, 0.2))) +
      scale_y_continuous(limits = c(0, NA), breaks = pretty_breaks(), expand = expansion(mult = c(0, 0.15))))

# view piN/piS values
distinct(dplyr::select(filter(codon_based_pi_results_bootstrap_LONG, sample_size == 100), model, dNdS, P_value))


# --------------------------------------------------------------------------------
# NEW in light of their updated response: Bloom et al. simulations, explore ben mutations and sliding windows

# count beneficials
mutations_filtered_countMutType <- filter(mutations_filtered) %>%
   group_by(model, sample_size, codon, mut_type) %>%
   summarize(
      mut_type_count = n()
   )

# initialize
mutations_filtered_countMutType$mut_class <- 'non-ben'

# mut classes for the Soni et al.
mutations_filtered_countMutType[mutations_filtered_countMutType$model != 'Emp/BloomDFE' & 
                                   mutations_filtered_countMutType$mut_type == 'm0',  ]$mut_class <- 'ben'

# mut classes for BloomDFE
mutations_filtered_countMutType[mutations_filtered_countMutType$model == 'Emp/BloomDFE' & 
                               mutations_filtered_countMutType$mut_type == 'm4',  ]$mut_class <- 'ben'

# COUNT BEN
mutations_filtered_countBen <- filter(mutations_filtered_countMutType, mut_class == 'ben') %>%
   group_by(model, sample_size, codon) %>%
   summarise(
      ben_count = n()
   )

# add rows for 0 beneficials
for (this_model in unique(mutations_filtered_countBen$model)) {
   # this_model <- 'Emp/BloomDFE'
   
   for (this_sample_size in unique(mutations_filtered_countBen$sample_size)) {
      # this_sample_size <- 100
      
      missing_codons <- setdiff(1:max(mutations_filtered$codon),
                                filter(mutations_filtered_countBen, model == this_model, sample_size == this_sample_size)$codon)
      missing_codons_count <- length(missing_codons)
      
      no_ben_rows <- tibble(model = rep(this_model, missing_codons_count),
                            sample_size = rep(this_sample_size, missing_codons_count),
                            codon = missing_codons,
                            ben_count = rep(0, missing_codons_count))
      
      mutations_filtered_countBen <- rbind(mutations_filtered_countBen, no_ben_rows)
   }
}

# sort by columns
mutations_filtered_countBen <- arrange(mutations_filtered_countBen, model, sample_size, codon)

# summarise
(mutations_filtered_countBen_SUMMARY <- mutations_filtered_countBen %>%
      group_by(model, sample_size) %>%
      summarize(
         total_codons = n(),
         ben_codons = sum(ben_count > 0),
         prop_codons_ben = ben_codons / total_codons,
         mean_ben = mean(ben_count)
      ))


# ------------------------------------------------------------------------------
# SLIDING WINDOWS

# CO-OP sliding window code from SNPGenie, https://github.com/chasewnelson/SNPGenie/blob/master/SNPGenie_sliding_windows.R

# set parameters
NUMERATOR <- "N"
DENOMINATOR <- "S"
WINDOW_SIZE <- 30
STEP_SIZE <- 1
NBOOTSTRAPS <- 1000
MIN_DEFINED_CODONS <- 10  # meaningless in this context
CORRECTION <- 'NONE'
NCPUS <- 10
# PREPEND_TO_OUTPUT <- ''

# ------------------------------------------------------------------------------
# BOOTSTRAP FUNCTION (dN - dS) for CODON UNIT, NO CORRECTION
dNdS_diff_boot_fun <- function(codon_results, numerator, denominator, num_replicates, num_cpus) {
   # # DUMMY
   #  codon_results <- codon_data_filtered
   #  numerator <- NUMERATOR
   #  denominator <- DENOMINATOR
   #  num_replicates <- NBOOTSTRAPS
   #  num_cpus <- NCPUS
   
   # Function for dN
   dN_function <- function(D, indices) {
      dN <- sum(D[indices, paste0(numerator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(numerator, "_sites")], na.rm = TRUE)
      return(dN)
   }
   
   # Function for dN
   dS_function <- function(D, indices) {
      dS <- sum(D[indices, paste0(denominator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(denominator, "_sites")], na.rm = TRUE)
      return(dS)
   }
   
   # Function for dN - dS
   dN_m_dS_function <- function(D, indices) {
      dN <- sum(D[indices, paste0(numerator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(numerator, "_sites")], na.rm = TRUE)
      dS <- sum(D[indices, paste0(denominator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(denominator, "_sites")], na.rm = TRUE)
      dN_m_dS <- dN - dS
      return(dN_m_dS)
   }
   
   # Function for dN/dS
   dN_over_dS_function <- function(D, indices) {
      dN <- sum(D[indices, paste0(numerator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(numerator, "_sites")], na.rm = TRUE)
      dS <- sum(D[indices, paste0(denominator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(denominator, "_sites")], na.rm = TRUE)
      dN_over_dS <- dN / dS
      return(dN_over_dS)
   }
   
   # CREATE FUNCTION FOR dN/dS TO CALCULATE ITS SE
   (dN <- sum(codon_results[ , paste0(numerator, "_diffs")], na.rm = TRUE) / sum(codon_results[ , paste0(numerator, "_sites")], na.rm = TRUE))
   (dS <- sum(codon_results[ , paste0(denominator, "_diffs")], na.rm = TRUE) / sum(codon_results[ , paste0(denominator, "_sites")], na.rm = TRUE))
   (dNdS <- dN / dS)
   
   # Run the BOOTSTRAPS
   # boot dN
   (boot_dN <- boot(data = codon_results, R = num_replicates, statistic = dN_function, parallel = 'multicore', ncpus = num_cpus))
   (dN <- boot_dN$t0)
   (boot_dN_SE <- sd(boot_dN$t))
   
   # boot dS
   (boot_dS <- boot(data = codon_results, R = num_replicates, statistic = dS_function, parallel = 'multicore', ncpus = num_cpus))
   (dS <- boot_dS$t0)
   (boot_dS_SE <- sd(boot_dS$t))
   
   # boot dN - dS
   (boot_dN_m_dS <- boot(data = codon_results, R = num_replicates, statistic = dN_m_dS_function, parallel = 'multicore', ncpus = num_cpus))
   (dN_m_dS <- boot_dN_m_dS$t0)
   (boot_dN_m_dS_SE <- sd(boot_dN_m_dS$t))
   (boot_dN_m_dS_Z <- dN_m_dS / boot_dN_m_dS_SE)
   (boot_dN_m_dS_P <- 2 * pnorm(-abs(boot_dN_m_dS_Z)))
   
   # boot dN/dS
   (boot_dN_over_dS <- boot(data = codon_results, R = num_replicates, statistic = dN_over_dS_function, parallel = 'multicore', ncpus = num_cpus))
   (dN_over_dS <- boot_dN_over_dS$t0)
   (boot_dN_over_dS_SE <- sd(boot_dN_over_dS$t))
   (boot_dN_over_dS_Z <- dN_over_dS / boot_dN_over_dS_SE)
   (boot_dN_over_dS_P <- 2 * pnorm(-abs(boot_dN_over_dS_Z)))
   
   ### NEW: ASL (acheived significance level)
   boot_dN_gt_dS_count <- sum(boot_dN_m_dS$t > 0)
   boot_dN_eq_dS_count <- sum(boot_dN_m_dS$t == 0)
   boot_dN_lt_dS_count <- sum(boot_dN_m_dS$t < 0)
   ASL_dN_gt_dS_P <- boot_dN_lt_dS_count / (boot_dN_gt_dS_count + boot_dN_eq_dS_count + boot_dN_lt_dS_count)
   ASL_dN_lt_dS_P <- boot_dN_gt_dS_count / (boot_dN_gt_dS_count + boot_dN_eq_dS_count + boot_dN_lt_dS_count)
   
   return(paste(num_replicates, dN, dS, dNdS, dN_m_dS, boot_dN_SE, boot_dS_SE, boot_dN_over_dS_SE, boot_dN_over_dS_P, 
                boot_dN_m_dS_SE, boot_dN_m_dS_P, 
                boot_dN_gt_dS_count, boot_dN_eq_dS_count, boot_dN_lt_dS_count, ASL_dN_gt_dS_P, ASL_dN_lt_dS_P,
                sep = "\t"))
}

# ADD COLUMNS TO DATA
codon_data <- codon_based_pi_codonMeans  # simply copy for consistent naming and no off-target ramifications
codon_data$codon_num <- codon_data$codon
RATIO_NAME <- paste0('pi', NUMERATOR, 'pi', DENOMINATOR)
codon_data$sw_ratio <- RATIO_NAME
codon_data$sw_start <- as.double(NA)
codon_data$sw_center <- as.double(NA)
codon_data$sw_end <- as.double(NA)
codon_data$sw_num_replicates <- as.integer(NA)
codon_data$sw_N_diffs <- as.double(NA)
codon_data$sw_S_diffs <- as.double(NA)
codon_data$sw_N_sites <- as.double(NA)
codon_data$sw_S_sites <- as.double(NA)
codon_data$sw_dN <- as.double(NA)
codon_data$sw_dS <- as.double(NA)
codon_data$sw_dNdS <- as.double(NA)
codon_data$sw_dN_m_dS <- as.double(NA)
codon_data$sw_boot_dN_SE <- as.double(NA)
codon_data$sw_boot_dS_SE <- as.double(NA)
codon_data$sw_boot_dN_over_dS_SE <- as.double(NA)
codon_data$sw_boot_dN_over_dS_P <- as.double(NA)
codon_data$sw_boot_dN_m_dS_SE <- as.double(NA)
codon_data$sw_boot_dN_m_dS_P <- as.double(NA)
codon_data$sw_boot_dN_gt_dS_count <- as.integer(NA)
codon_data$sw_boot_dN_eq_dS_count <- as.integer(NA)
codon_data$sw_boot_dN_lt_dS_count <- as.integer(NA)
codon_data$sw_ASL_dN_gt_dS_P <- as.double(NA)
codon_data$sw_ASL_dN_lt_dS_P <- as.double(NA)
codon_data$sw_ASL_dNdS_P <- as.double(NA)

# arrange
codon_data <- dplyr::select(codon_data, model, sample_size, codon, codon_num, everything())
codon_data <- dplyr::arrange(codon_data, model, sample_size, codon_num)

# PERFORM SLIDING WINDOW - just do 'Emp/BloomDFE' because we've run out of time <== takes ~2 hours per model/sample size
for(this_model in c('Emp/BloomDFE')) {  # in unique(codon_data$model)) {
   # this_model <- "Emp/Flynn1pct"
   
   for(this_sample_size in c(100)) {  # in unique(codon_data$sample_size)) {
      # this_sample_size <- 100
      
      for(i in seq(from = min(codon_data$codon_num), to = (max(codon_data$codon_num) - WINDOW_SIZE + 1), by = STEP_SIZE)) {  # always 1..10000
      # for(i in seq(from = 9121, to = (max(codon_data$codon_num) - WINDOW_SIZE + 1), by = STEP_SIZE)) {  # always 1..10000
         # i <- 1 # 111
         # cat(i, "")
         
         # Extract window; analyze
         window_codon_data <- filter(codon_data, model == this_model, sample_size == this_sample_size, codon_num >= i, codon_num <= (i + WINDOW_SIZE - 1))
         # view(window_codon_data)
         lowest_codon_num <- min(window_codon_data$codon_num)
         highest_codon_num <- max(window_codon_data$codon_num)
         
         if(nrow(window_codon_data) > 1) {
            
            # WARNING if missing codon numbers
            if(nrow(window_codon_data) != WINDOW_SIZE) {
               cat("### WARNING: codon", i, "is missing.")
               cat("### If this was unexpected or is common, there may be problems with your data.")
            }
            
            ### BOOTSTRAP EACH RATIO FOR THIS WINDOW. Columns:
            #[1] num_replicates
            #[2] dN
            #[3] dS
            #[4] dNdS
            #[5] dN_m_dS
            #[6] boot_dN_SE
            #[7] boot_dS_SE
            #[8] boot_dN_over_dS_SE
            #[9] boot_dN_over_dS_P
            #[10] boot_dN_m_dS_SE
            #[11] boot_dN_m_dS_P
            #[12] boot_dN_gt_dS_count
            #[13] boot_dN_eq_dS_count
            #[14] boot_dN_lt_dS_count
            #[15] ASL_dN_gt_dS_P
            #[16] ASL_dN_lt_dS_P
            
            # dN/dS bootstrap with appropriate correction
            boot_dNdS <- NA
            if(CORRECTION == "JC") {
               boot_dNdS <- dNdS_diff_boot_fun_JC(window_codon_data, NUMERATOR, DENOMINATOR, NBOOTSTRAPS, NCPUS)
            } else {
               boot_dNdS <- dNdS_diff_boot_fun(window_codon_data, NUMERATOR, DENOMINATOR, NBOOTSTRAPS, NCPUS)
            }
            boot_dNdS <- str_split(string = boot_dNdS, pattern = '\t')[[1]]
            
            # bootstrap results
            num_replicates <- as.integer(boot_dNdS[1])
            dN <- as.numeric(boot_dNdS[2])
            dS <- as.numeric(boot_dNdS[3])
            dNdS <- as.numeric(boot_dNdS[4])
            dN_m_dS <- as.numeric(boot_dNdS[5])
            boot_dN_SE <- as.numeric(boot_dNdS[6])
            boot_dS_SE <- as.numeric(boot_dNdS[7])
            boot_dN_over_dS_SE <- as.numeric(boot_dNdS[8])
            boot_dN_over_dS_P <- as.numeric(boot_dNdS[9])
            boot_dN_m_dS_SE <- as.numeric(boot_dNdS[10])
            boot_dN_m_dS_P <- as.numeric(boot_dNdS[11])
            boot_dN_gt_dS_count <- as.integer(boot_dNdS[12])
            boot_dN_eq_dS_count <- as.integer(boot_dNdS[13])
            boot_dN_lt_dS_count <- as.integer(boot_dNdS[14])
            ASL_dN_gt_dS_P <- as.numeric(boot_dNdS[15])
            ASL_dN_lt_dS_P <- as.numeric(boot_dNdS[16])
            ASL_dNdS_P <- 1
            
            if(! is.na(ASL_dN_gt_dS_P) && ! is.na(ASL_dN_lt_dS_P) && ASL_dN_gt_dS_P < ASL_dN_lt_dS_P) {
               ASL_dNdS_P <- 2 * ASL_dN_gt_dS_P
            } else if(! is.na(ASL_dN_gt_dS_P) && ! is.na(ASL_dN_lt_dS_P)) {
               ASL_dNdS_P <- 2 * ASL_dN_lt_dS_P
            }
            
            if(ASL_dNdS_P == 0) {
               ASL_dNdS_P <- 1 / NBOOTSTRAPS
            }
            
            # Add to table
            codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_start <- lowest_codon_num
            codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_center <- (lowest_codon_num + highest_codon_num) / 2
            codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_end <- highest_codon_num
            codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_num_replicates <- num_replicates
            codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_N_diffs <- sum(unname(unlist(window_codon_data[ , paste0(NUMERATOR, '_diffs')])), na.rm = TRUE)
            codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_S_diffs <- sum(unname(unlist(window_codon_data[ , paste0(DENOMINATOR, '_diffs')])), na.rm = TRUE)
            codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_N_sites <- sum(unname(unlist(window_codon_data[ , paste0(NUMERATOR, '_sites')])), na.rm = TRUE)
            codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_S_sites <- sum(unname(unlist(window_codon_data[ , paste0(DENOMINATOR, '_sites')])), na.rm = TRUE)
            codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_dN <- dN
            codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_dS <- dS
            codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_dNdS <- dNdS
            codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_dN_m_dS <- dN_m_dS
            codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_boot_dN_SE <- boot_dN_SE
            codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_boot_dS_SE <- boot_dS_SE
            codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_boot_dN_over_dS_SE <- boot_dN_over_dS_SE
            codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_boot_dN_over_dS_P <- boot_dN_over_dS_P
            codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_boot_dN_m_dS_SE <- boot_dN_m_dS_SE
            codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_boot_dN_m_dS_P <- boot_dN_m_dS_P
            codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_boot_dN_gt_dS_count <- boot_dN_gt_dS_count
            codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_boot_dN_eq_dS_count <- boot_dN_eq_dS_count
            codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_boot_dN_lt_dS_count <- boot_dN_lt_dS_count
            codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_ASL_dN_gt_dS_P <- ASL_dN_gt_dS_P
            codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_ASL_dN_lt_dS_P <- ASL_dN_lt_dS_P
            codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_ASL_dNdS_P <- ASL_dNdS_P
            
         } # else leave it as NA
      } # end last window
   }
}

# SAVE 
# write_tsv(filter(codon_data, model == "Emp/BloomDFE", sample_size == 100), "codon_data_sw_results_BloomDFE_100.tsv")

# SUBSET or RELOAD
# codon_data_BloomDFE_100 <- filter(codon_data, model == "Emp/BloomDFE", sample_size == 100)
codon_data_BloomDFE_100 <- read_tsv("data/codon_data_sw_results_BloomDFE_100.tsv")


# ------------------------------------------------------------------------------
# SELECTION CANDIDATES

# define positive selection candidates & conduct 'back of napkin' false-positive rate optimization
codon_data_BloomDFE_100$pos_sel_candidate <- codon_data_BloomDFE_100$sw_dN_m_dS > 0 & codon_data_BloomDFE_100$sw_boot_dN_m_dS_P < 0.00605  # from data exploration

# JOIN ben count
(mutations_filtered_countBen_BloomDFE_100 <- filter(mutations_filtered_countBen, model == "Emp/BloomDFE", sample_size == 100))

# sliding window of ben_count
mutations_filtered_countBen_BloomDFE_100 <- mutations_filtered_countBen_BloomDFE_100 %>%
   mutate(ben_sw30_sum = rollsum(ben_count, k = 30, align = 'left', fill = NA))

# ---
# JOIN
(codon_data_BloomDFE_100_JOINED <- left_join(codon_data_BloomDFE_100, mutations_filtered_countBen_BloomDFE_100, by = c('model', 'sample_size', 'codon')))

# ----------
# TRUTH is >= 1 ben

# summarize
codon_data_BloomDFE_100_JOINED %>% 
   group_by(ben_sw30_sum >= 1) %>%
   summarise(
      count = n()
   )
# `ben_sw30_sum >= 1` count
# 1 FALSE                1031
# 2 TRUE                 8940
# 3 NA                     29
# FP of random choice: 1031 / (1031 + 8940) == 0.1033999

filter(codon_data_BloomDFE_100_JOINED, sw_dN_m_dS > 0) %>% 
   group_by(pos_sel_candidate, ben_sw30_sum >= 1) %>%
   summarise(
      count = n()
   )
# pos_sel_candidate `ben_sw30_sum >= 1` count
# 1 FALSE             FALSE                 679
# 2 FALSE             TRUE                 3875
# 3 TRUE              TRUE                   80

# is that better than chance?
prop.test(x = c(8940, 80), n = c((1031 + 8940), 80))  # YES, p-value = 0.004358 QED


# --------------------------------------------------------------------------------
# SUMMARIZE SIMULATION METRICS
mutations_per_rep <- arrange(mutations_per_rep, model, sample_size, rep)  # these have ALREADY BEEN FILTERED for MAF >= 2.5%

# PLOT mutations per sample
mutations_per_rep$model <- factor(mutations_per_rep$model, levels = c("Neu/DFE1", "Sel/DFE1", "Neu/DFE2", "Sel/DFE2", "Emp/Flynn1pct", "Emp/Flynn10pct", "Emp/BloomDFE"))

# PLOT
(mutations_per_rep_BOXPLOT <- ggplot(filter(mutations_per_rep, sample_size == 100), aes(x = model, y = mut_count)) +
  
  # BOXPLOT
  geom_boxplot(outlier.size = 0.5, outlier.shape = 21, outlier.fill = NA, outlier.stroke = 0.2, linewidth = 0.25, fatten = 1) +
  stat_summary(fun = 'mean', geom = 'point', shape = 23, size = 1, fill = 'black') +
  
  # OBSERVED MEAN
  geom_hline(yintercept = 10.5, color = 'blue', linetype = 'dashed', linewidth = .2) +
  
  # OBSERVED MEDIAN
  geom_hline(yintercept = 5, color = 'blue', linetype = 'dashed', linewidth = .2) +
  
  theme_bw() +
  scale_y_continuous(limits = c(0, NA), breaks = pretty_breaks(), expand = expansion(mult = c(0, 0.1))))

# SAVE SOURCE
# write_tsv(mutations_per_rep, 'data/Figure1c_source.txt')

