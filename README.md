# Response to Soni et al.
Modified scripts, analysis code, and source data for response to Soni et al. 2023.

## Scripts

Two SLiM scripts (with and without a beneficial mutation introduced) were downloaded from [Soni et al.](https://github.com/vivaksoni/Gu_etal_2023_response/tree/main/scripts) (accessed 2023/06/16).

Updated versions that print additional output (`*_addOutput.slim`) were saved with the following lines of Eidos code added:

```eidos
// additional output
outputSample.genome1.output(filePath = getwd() + "/" + OUTPUTSTEM + "_100.out");
sampledIndividuals = sample(p1.individuals, 1000);
sampledIndividuals.genome1.output(filePath = getwd() + "/" + OUTPUTSTEM + "_1000.out");
```

Results in our response are based on the `*_100.out` data only, analogous to a uniform coverage of 100 effective sequencing reads, matching the method of Soni et al.

Each script was run with a weakly or strongly deleterious mutation background as follows:

```bash
# Weakly deleterious background with no beneficial mutations (Weak / -)
for i in $(seq 1 100); do slim -d GENOMESIZE=30000 -d MU=2.135e-6 -d INIT=1 -d K=1e5 -d REPRO=1 -d RUNTIME=168 -d R=5.5e-5 -d XI=0 -d BURSTN=100 -d "OUTPUTSTEM='results/DFE/DFE1_rep$i'" -d d_f0=0.1 -d d_f1=0.7 -d d_f2=0.1 -d d_f3=0.1 -d simID="$i" sc2_DFE_addOutput.slim > results/DFE/DFE1_rep${i}.log; done;

# Strongly deleterious background with no beneficial mutations (Strong / -)
for i in $(seq 1 100); do slim -d GENOMESIZE=30000 -d MU=2.135e-6 -d INIT=1 -d K=1e5 -d REPRO=1 -d RUNTIME=168 -d R=5.5e-5 -d XI=0 -d BURSTN=100 -d "OUTPUTSTEM='results/DFE/DFE2_rep$i'" -d d_f0=0.1 -d d_f1=0.1 -d d_f2=0.1 -d d_f3=0.7 -d simID="$i" sc2_DFE_addOutput.slim > results/DFE/DFE2_rep${i}.log; done;

# Weakly deleterious background with one beneficial mutation (Weak / +)
for i in $(seq 1 100); do slim -d GENOMESIZE=30000 -d MU=2.135e-6 -d INIT=1 -d K=1e5 -d REPRO=1 -d RUNTIME=168 -d R=5.5e-5 -d XI=0 -d BURSTN=100 -d "OUTPUTSTEM='results/DFE_beneficial/DFE1_rep$i'" -d d_f0=0.1 -d d_f1=0.7 -d d_f2=0.1 -d d_f3=0.1 -d simID="$i" sc2_DFE_beneficial_addOutput.slim > results/DFE_beneficial/DFE1_rep${i}.log; done;

# Strongly deleterious background with one beneficial mutation (Strong / +)
for i in $(seq 1 100); do slim -d GENOMESIZE=30000 -d MU=2.135e-6 -d INIT=1 -d K=1e5 -d REPRO=1 -d RUNTIME=168 -d R=5.5e-5 -d XI=0 -d BURSTN=100 -d "OUTPUTSTEM='results/DFE_beneficial/DFE2_rep$i'" -d d_f0=0.1 -d d_f1=0.1 -d d_f2=0.1 -d d_f3=0.7 -d simID="$i" sc2_DFE_beneficial_addOutput.slim > results/DFE_beneficial/DFE2_rep${i}.log; done;
```

Data corresponding to the `Mutation` and `Genome` blocks of SLiM output were extracted from the output files at the command line as follows:

```bash
# Mutation data - no beneficial mutation introduced (-)
grep -E '^[0-9]' results/DFE/*.out > DFE_mutations.txt

# Mutation data - one beneficial mutation introduced (+)
grep -E '^[0-9]' results/DFE_beneficial/*.out > DFE_beneficial_mutations.txt

# Genome data - no beneficial mutation introduced (-)
grep -E '^p' results/DFE/*.out | sed -E 's/ A/\tA/g' | sed -E 's/A /A\t/g' | sed -E 's/ /,/g' > DFE_genomes.txt

# Genome data - one beneficial mutation introduced (+)
grep -E '^p' results/DFE_beneficial/*.out | sed -E 's/ A/\tA/g' | sed -E 's/A /A\t/g' | sed -E 's/ /,/g' > DFE_beneficial_genomes.txt
```

These raw output files are available in the `/data/` directory.

All subsequent wrangling and analysis was performed manually in the R script `Soni_response.R`, R version 4.2.2 (2022-10-31), RStudio 2023.06.0+421. The following library packages were used: base, boot, datasets, dplyr, forcats, ggplot2, graphics, grDevices, lubridate, methods, purrr, RColorBrewer, readr, scales, stats, stringr, tibble, tidyr, tidyverse, utils.
