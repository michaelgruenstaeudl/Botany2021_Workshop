### Botany 2021 - Workshop
Workshop W14 by Michael Gruenstaeudl during the conference "Botany 2021 Virtual"

# Assessing sequence coverage and inverted repeat annotations among complete plastid genomes

---

## Table of Contents
1. [Abstract](https://github.com/michaelgruenstaeudl/Botany2021_Workshop/blob/main/doc/abstract.md)
2. [Introduction](#introduction)
3. [Survey of IR annotations of archived plastid genomes](#survey-of-ir-annotations)
      1. Installation of airpg
      2. Application of airpg
3. [Visualization of sequencing depth and evenness of complete plastid genomes](#visualization-of-sequencing-depth-and-evenness)
      1. Installation of PACVr
      2. Mapping of sequence reads
      3. Application of PACVr

---

## Introduction
This workshop is intended to illustrate the application of computational methods that are helpful in assessing the quality of complete plastid genome sequences. Specifically, the workshop will illustrate the application of two software tools that (i) evaluate the annotations of inverted repeat annotations, and (ii) evaluate and visualize sequencing depth and evenness of the genome records, respectively. The workshop is intended for researchers with prior experience in the assembly and annotation of plastid genomes and assumes that participants have already assembled and annotated at least one plastid genome themselves. Participants will be guided through the application of both tools in a step-by-step process. All instructions presented in this workshop have been designed and customized for the UNIX command line (e.g., bash) and should be executed on a UNIX-compatible operating system (OS-X or Linux).

---

## Survey Of IR Annotations
The inverted repeats (IRs) are characteristic features of the great majority of land plant plastid genomes. High-quality plastid genome records should contain sequence annotations for these features. The software [airpg](https://pypi.org/project/airpg/) has been designed to automatically evaluate the presence of complete and correct IR sequence annotations in complete plastid genome records archived on GenBank.

Software needed: Python 3, [airpg](https://pypi.org/project/airpg/)

#### 1. Installation of airpg
```bash
$ pip install airpg
```
#### 2. Application of airpg
+ _Objective_: What proportion of all complete plastid genomes of all moss lineages (i.e., liverworts, hornworts, and mosses) submitted to NCBI Nucleotide since the beginning of 2000 does not have complete IR annotations?<br>
+ _Time needed_: ca. 8 min.

Identify the plastid genomes:
```bash
$ airpg_identify.py -q "complete genome[TITLE] AND \
(chloroplast[TITLE] OR plastid[TITLE]) AND \
2000/01/01:2021/05/31[PDAT] NOT partial[TITLE] \
AND (Marchantiophyta[ORGN] OR Bryophyta[ORGN] \
OR Anthocerotophyta[ORGN])" \
 -o airpg_SimpleExample_output1.tsv
 ```

Analyze their IR annotations:
 ```bash
 $ airpg_analyze.py -i airpg_SimpleExample_output1.tsv \
 -m john.smith@example.com -o airpg_SimpleExample_output2.tsv
  ```

Visualize the accumulation of plastid genomes of all moss lineages with and without complete IR annotations over time.
```bash
# Get number of genome records
$ NL=$(wc -l airpg_SimpleExample_output1.tsv | awk '{print $1}')
$ echo "$NL-1" | bc
# 76

# Get submission dates of oldest and newest genome record
$ awk -F'\t' '{print $6}' airpg_SimpleExample_output1.tsv | \
grep "^2" | sort -n | awk 'NR==1; END{print}'
# 2003-02-04
# 2021-04-24

# Adjust script airpg_SimpleExample_visualization.R manually and then run
$ Rscript /path_to_git_folder/extras/airpg_SimpleExample_visualization.R
 ```
![](https://github.com/michaelgruenstaeudl/Botany2021_Workshop/blob/main/extras/airpg_SimpleExample_visualization.png)

---

## Visualization Of Sequencing Depth And Evenness
The quality of complete plastid genome records often correlates with the records' sequence coverage. Specifically, high-quality plastid genome records often exhibit considerable sequencing depth and high sequencing evenness. The software [PACVr](https://cran.r-project.org/package=PACVr) has been designed to automatically evaluate and visualize sequencing depth and evenness of the genome records so that users can assess both coverage indicators for given plastid genome records.

Software: R, [PACVr](https://cran.r-project.org/package=PACVr)

#### Installation of PACVr
```bash
$ R
install.packages("PACVr")
```
#### Mapping of sequence reads
To measure sequencing depth and evenness of a plastid genome, the sequence reads that were originally used to assemble that genome must be mapped against it. The process of read mapping is explained [here](https://github.com/michaelgruenstaeudl/Botany2021_Workshop/blob/main/doc/mapping_reads_to_plastid_genome.md).

#### Application of PACVr
```R
if (!require("pacman")) {
  install.packages("pacman")
}
library(pacman)
pacman::p_load(RCircos,
               genbankr,
               devtools,
               tidyverse,
               vroom)
pacman::p_install_gh("nilsj9/PACVr")
library(PACVr)

#### FUNCTIONS ####
# ---------------------------------------------------------------------------- #
evennessScore <- function(coverage) {
  coverage_mean <- round(mean(coverage))
  D2 <- coverage[coverage <= coverage_mean]
  E <- 1 - (length(D2) - sum(D2) / coverage_mean) / length(coverage)
  return(E)
}

effSize_conv <- function(x) {
  return(2 * x / (sqrt(1 - x ^ 2)))
}


#### SAMPLE LIST ####
# ---------------------------------------------------------------------------- #
sample_file <- list.files(
  path = "./",
  recursive = TRUE,
  pattern = "^samples_list.csv",
  full.names = TRUE
)

if (length(file.exists(sample_file)) == 0) {
  acquire_data()
  data <- read.csv("samples_list.csv")
}

sample_list <- read.csv(sample_file[1], header = TRUE)

#### WORKING DIR ####
# ---------------------------------------------------------------------------- #
w_dir <- "/PATH/TO/WORKING/DIRECTORY/"

#### INPUT FILES ####
# ---------------------------------------------------------------------------- #
# coverage files
bed_files <- list.files(
  path = paste0(w_dir,"/samples"),
  full.names = TRUE,
  recursive = TRUE,
  pattern = ".bed"
)
# meta data files
metadata_files <- list.files(
  path = paste0(w_dir,"/samples"),
  full.names = TRUE,
  recursive = TRUE,
  pattern = "metadata.csv"
)
# annotation files
gb_files <- list.files(
  path = paste0(w_dir,"/samples"),
  full.names = TRUE,
  recursive = TRUE,
  pattern = "[0-9].gb"
)

### READ META DATA ###
# ---------------------------------------------------------------------------- #
# lineage information
sample_list$Lineage <-
  sapply(gb_files, function(x)
    gsub(".*Embryophyta ", "", paste(
      trimws(parseGenBank(x)$SOURCE$lineage), collapse = " "
    )))

meta_data <- purrr::map_df(metadata_files, ~ vroom::vroom(
  .x,
  delim = ",",
  col_types = c(
    Number_N = "i",
    Mismatches = "i",
    avgLength = "i",
    Model = "c",
    `Assembly Method` = "c",
    `Sequencing Technology` = "c"
  )
))
meta_data <- cbind(sample_list, meta_data)

### READ COVERAGE DATA ###
# ---------------------------------------------------------------------------- #
# coverage information
coverage_data <- purrr::map(bed_files, ~ vroom::vroom(
  .x,
  col_types = c(
    Chromosome = "c",
    chromStart = "i",
    chromEnd = "i",
    coverage = "i",
    lowCoverage = "c",
    gene = "c"
  )
))
# remove regions < 250
coverage_data <-
  lapply(coverage_data, function(x)
    x[-which(x$chromEnd - x$chromStart < 249), ])

# Assign names to list elements
names(coverage_data) <- unlist(unname(lapply(sample_list$Accession,
                                             function(x)
                                               c(
                                                 paste(x, "_genes", sep = ""),
                                                 paste(x, "_noncoding", sep = ""),
                                                 paste(x, "_regions", sep = "")
                                               ))))

# Split list by coverage type
coverage_regions  <-
  coverage_data[grep("regions", names(coverage_data))]
coverage_genes <- coverage_data[grep("genes", names(coverage_data))]
coverage_noncoding <-
  coverage_data[grep("noncoding", names(coverage_data))]

# Add coverage depth data
# plastid partition
regions <- lapply(coverage_regions, function(x)
  x %>%
    tidyr::separate_rows(Chromosome) %>%
    dplyr::group_by(Chromosome) %>%
    dplyr::summarise(
      lowCoverage = sum(lowCoverage == "*", na.rm = TRUE),
      .groups = "drop"
    ))
regions <- lapply(regions, function(x)
  x %>%
    tidyr::spread(., Chromosome, lowCoverage))
regions <- select(dplyr::bind_rows(regions), IRa, IRb, LSC, SSC)
meta_data <- dplyr::bind_cols(meta_data, regions)

# coding/non-coding partition
meta_data$coding <-
  unlist(unname(lapply(coverage_genes, function(x)
    sum(
      !is.na(x$lowCoverage == "*")
    ))))
meta_data$noncoding <-
  unlist(unname(lapply(coverage_noncoding, function(x)
    sum(!is.na(x$lowCoverage == "*")))))

meta_data <- rename(meta_data, Ns = Number_N)

### CALCULATE E-SCORE ###
# ---------------------------------------------------------------------------- #
# Calculate coverage Escore
meta_data$Escore <- unname(sapply(coverage_regions,
                                  function(x)
                                    evennessScore(x$coverage)))
```


---

#### 3. airpg - Complex example
Foo bar baz

```bash
# Get total number of families represented by these records
$ awk -F'\t' '{print $11}' airpg_ComplexExample_output1.tsv | \
awk -F';' '{ for(i=1;i<=NF;i++) print $i }' | grep "aceae" | \
sort -u | wc -l
# 308
```

#1>>airpg_SimpleExample_output1.log 2>&1
