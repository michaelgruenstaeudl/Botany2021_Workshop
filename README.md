### Botany 2021 - Workshop
Workshop W14 by Michael Gruenstaeudl during the conference "Botany 2021 Virtual"

# Assessing sequence coverage and inverted repeat annotations among complete plastid genomes

---

## Table of Contents
1. [Abstract](https://github.com/michaelgruenstaeudl/Botany2021_Workshop/blob/main/doc/abstract.md)
2. Introduction
3. [Survey of IR annotations of archived plastid genomes](#survey-of-ir-annotations)
      1. Installation of airpg
      2. Application of airpg
3. [Visualization of sequencing depth and evenness of complete plastid genomes](#visualization-of-sequencing-depth-and-evenness)
      1. Installation of PACVr
      2. Preparation of genomic data
      3. Application of PACVr

---

## Introduction
This workshop is intended to ...

The instructions of this workshop are designed to be executed on a Linux command line shell (e.g., bash).

---

## Survey Of IR Annotations
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
Software: R, [PACVr](https://cran.r-project.org/package=PACVr)

#### Installation of PACVr
```bash
$ R
install.packages("PACVr")
```
#### Preparation of genomic data
To measure sequencing depth and evenness of a plastid genome, the sequence reads that were originally used to assemble that genome must be mapped against it. The process of read mapping is explained [here](https://github.com/michaelgruenstaeudl/Botany2021_Workshop/blob/main/doc/mapping_reads_to_plastid_genome.md).

#### Application of PACVr
Foo bar baz


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
