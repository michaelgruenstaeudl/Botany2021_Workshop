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


#### Application of airpg - More tutorials

More complex tutorials regarding the application of airpg can be found [here](https://github.com/michaelgruenstaeudl/airpg) as well as - in a platform-independent way - [here](https://codeocean.com/capsule/6723913/tree/v1).


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

#### Application of PACVr - Sequencing depth only

```R
library(PACVr)

# Specify input files
gbkFile <- system.file("extdata", "NC_045072/NC_045072.gb", package="PACVr")
bamFile <- system.file("extdata", "NC_045072/NC_045072_PlastomeReadsOnly.sorted.bam",package="PACVr")

# Specify output file
outFile <- "NC_045072_AssemblyCoverage_viz.pdf"

# Run PACVr

PACVr.complete(gbk.file=gbkFile, bam.file=bamFile, windowSize=250,
    logScale=FALSE, threshold=0.5, syntenyLineType=3,relative=TRUE, 
    textSize=0.5, output=outFile)
```
![](https://github.com/michaelgruenstaeudl/Botany2021_Workshop/blob/main/extras/NC_045072_AssemblyCoverage_viz.png)


#### Application of PACVr - Sequencing depth and sequencing evenness

A more complex tutorial regarding the application of PACVr can be found [here](https://github.com/nilsj9/PlastidSequenceCoverage).

