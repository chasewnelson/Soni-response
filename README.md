# Response to Soni et al.
Scripts, analysis code, input data, and intermediate files for response to Soni et al. 2024.

## <a name="citation"></a>Citation

When using this repository, please refer to and cite:

>Nelson CW, Poon LLM, Gu H. 2024. Reply to: Population genetic considerations regarding the interpretation of within-patient SARS-CoV-2 polymorphism data. *Nature Communications*, **15:** 3239. DOI: [10.1038/s41467-024-46262-3](https://doi.org/10.1038/s41467-024-46262-3)

and this page:

>https://github.com/chasewnelson/Soni-response

## <a name="description"></a>Description

DFE\* and Flynn\* SLiM scripts were downloaded from [Soni et al.](https://github.com/vivaksoni/Gu_etal_2023_response/tree/main/scripts).

Updated versions that print additional output (`*_addOutput.slim`) were saved with the following lines of Eidos code added:

```eidos
// additional output
outputSample.genome1.output(filePath = getwd() + "/" + OUTPUTSTEM + "_100.out");
sampledIndividuals = sample(p1.individuals, 1000);
sampledIndividuals.genome1.output(filePath = getwd() + "/" + OUTPUTSTEM + "_1000.out");
```

Results in our response are based on the `*_100.out` data only, analogous to a uniform coverage of 100 effective sequencing reads, matching the method of Soni et al.

Each DFE\* script (accessed 2023/06/16) was run with a weakly or strongly deleterious mutation background as follows:

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

Flynn\* (accessed 2023/09/26) and Bloom scripts were run as follows (note that Bloom DFE specifications are hard-coded; see script):

```bash
# Flynn 1% beneficial
for i in $(seq 1 100); do slim -d GENOMESIZE=30000 -d MU=2.135e-6 -d INIT=1 -d K=1e3 -d REPRO=1 -d RUNTIME=168 -d R=5.5e-5 -d XI=0 -d BURSTN=100 -d "OUTPUTSTEM='results/Flynn1pct/Flynn1pct_rep$i'" -d d_f0=0.542 -d d_f1=0.112 -d d_f2=0.02 -d d_f3=0.326 -d d_fb=0.01 -d simID="$i" sc2_Flynn_etal_DFE_addOutput.slim > results/Flynn1pct/Flynn1pct_rep${i}.log; done;

# Flynn 10% beneficial
for i in $(seq 1 100); do slim -d GENOMESIZE=30000 -d MU=2.135e-6 -d INIT=1 -d K=1e3 -d REPRO=1 -d RUNTIME=168 -d R=5.5e-5 -d XI=0 -d BURSTN=100 -d "OUTPUTSTEM='results/Flynn10pct/Flynn10pct_rep$i'" -d d_f0=0.445 -d d_f1=0.112 -d d_f2=0.02 -d d_f3=0.326 -d d_fb=0.097 -d simID="$i" sc2_Flynn_etal_DFE_addOutput.slim > results/Flynn10pct/Flynn10pct_rep${i}.log; done;

# Bloom & Neher DFE
for i in $(seq 1 100); do slim -d GENOMESIZE=30000 -d MU=2.135e-6 -d INIT=1 -d K=1e3 -d REPRO=1 -d RUNTIME=168 -d R=5.5e-5 -d "OUTPUTSTEM='results/BloomDFE/BloomDFE_rep$i'" -d simID="$i" sc2_BloomDFE.slim > results/BloomDFE/BloomDFE_rep${i}.log; done;
```

For each simulation above, data corresponding to the `Mutation` and `Genome` blocks of SLiM output were extracted from the output files at the command line as follows (replacing 'DFE' with the appropriate simulation type):

```bash
# Mutation data
grep -E '^[0-9]' results/DFE/*.out > DFE_mutations.txt

# Genome data
grep -E '^p' results/DFE/*.out | sed -E 's/ A/\tA/g' | sed -E 's/A /A\t/g' | sed -E 's/ /,/g' > DFE_genomes.txt
```

All raw output, random seeds used, and intermediate data files are available in the `/data/` directory of this repository. The file `aamut_fitness_all.csv` should be obtained from [Bloom & Neher](https://academic.oup.com/ve/article/9/2/vead055/7265011) at the study [GitHub link](https://raw.githubusercontent.com/jbloomlab/SARS2-mut-fitness/main/results/aa_fitness/aamut_fitness_all.csv) (public_2023-10-01 dataset; accessed 2023/10/05). Figure source data are provided as a Source Data file in the publication supplementary material.

All wrangling and analysis were performed manually in the R scripts `Soni-response.R` and `Bloom-SC2-DFE.R`, R version 4.2.2 (2022-10-31) and RStudio 2023.06.0+421, or R version 4.3.1 (2023-06-16) and RStudio 2023.06.1+524. R scripts are meant to be run interactively, line-by-line in RStudio. Figures were produced in R and annotated in PowerPoint. The following R libraries were used: base, boot, data.table, datasets, dplyr, forcats, ggplot2, graphics, grDevices, lubridate, MASS, methods, purrr, RColorBrewer, readr, scales, stats, stringr, tibble, tidyr, tidyverse, utils, zoo. 

## <a name="contact"></a>Contact
If you have questions about this repository, please click on the <a target="_blank" href="https://github.com/chasewnelson/Soni-response/issues">Issues</a> tab at the top of this page and begin a new thread, so that others might benefit from the discussion.
