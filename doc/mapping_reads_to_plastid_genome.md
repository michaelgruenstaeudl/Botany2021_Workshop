### Mapping sequence reads to a plastid genome
To measure sequencing depth and evenness of a plastid genome, the sequence reads that were originally used to assemble that genome must be mapped against it. This page explains the process of mapping sequence reads against a given plastid genome. To start, you need to know the GenBank accession number of the plastid genome ...

---

## Prerequisites

Software needed:
+ ncbi-entrez-direct (version number here)
+ sra-toolkit (>2.10.8)
+ samtools (>1.10)
+ BOWTIE2 (>2.3.4.1)
+ TRIMMOMATIC (version number here)

Other needed:
+ GenBank accession number

---




# DEFINITIONS
SAMPLE=$1
ACCESSION=$2
SRA=$3

# FOLDERS
mkdir -p $SAMPLELOC/$SAMPLE/ReadMapping

# download GenBank
${NCBIENT}esearch -db nucleotide -query $ACCESSION | ${NCBIENT}efetch -format gb > $SAMPLELOC/$SAMPLE/${ACCESSION}.gb

# download reference fasta
${NCBIENT}esearch -db nucleotide -query $ACCESSION | ${NCBIENT}efetch -format fasta > $SAMPLELOC/$SAMPLE/${ACCESSION}_ref.fasta

# download reads
cd $SAMPLELOC/$SAMPLE
${SRAT}prefetch --max-size 50000000 $SRA
cd /scratch/nilsj24/
${SRAT}fasterq-dump.2.10.8 --split-3 --skip-technical $SAMPLELOC/$SAMPLE/$SRA/$SRA.sra -O $SAMPLELOC/$SAMPLE/ReadMapping

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE $SAMPLELOC/$SAMPLE/ReadMapping/${SRA}.sra_1.fastq $SAMPLELOC/$SAMPLE/ReadMapping/${SRA}.sra_2.fastq  $SAMPLELOC/$SAMPLE/ReadMapping/${SRA}_1_PE.fastq $SAMPLELOC/$SAMPLE/ReadMapping/${SRA}_1_SE.fastq $SAMPLELOC/$SAMPLE/ReadMapping/${SRA}_2_PE.fastq $SAMPLELOC/$SAMPLE/ReadMapping/${SRA}_2_SE.fastq ILLUMINACLIP:/$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36

# conduct mapping of reads
REF=$SAMPLELOC/$SAMPLE/${ACCESSION}_ref.fasta
READSPER1=$SAMPLELOC/$SAMPLE/ReadMapping/${SRA}_1_PE.fastq
READSPER2=$SAMPLELOC/$SAMPLE/ReadMapping/${SRA}_2_PE.fastq

mkdir -p $SAMPLELOC/$SAMPLE/ReadMapping/db
${BOWTIE}bowtie2-build $REF $SAMPLELOC/$SAMPLE/ReadMapping/db/$ACCESSION
${BOWTIE}bowtie2 -x $SAMPLELOC/$SAMPLE/ReadMapping/db/$ACCESSION -1  $READSPER1 -2 $READSPER2 -S $SAMPLELOC/$SAMPLE/ReadMapping/${ACCESSION}_mapping.sam
${SAM}samtools view -Sb -F 0x04 $SAMPLELOC/$SAMPLE/ReadMapping/${ACCESSION}_mapping.sam > $SAMPLELOC/$SAMPLE/ReadMapping/${ACCESSION}_mapping_OneMoreLocations.bam
${SAM}samtools sort $SAMPLELOC/$SAMPLE/ReadMapping/${ACCESSION}_mapping_OneMoreLocations.bam > $SAMPLELOC/$SAMPLE/${ACCESSION}_mapping_OneMoreLocations.sorted.bam
${SAM}samtools index $SAMPLELOC/$SAMPLE/${ACCESSION}_mapping_OneMoreLocations.sorted.bam

rm -rf $SAMPLELOC/$SAMPLE/ReadMapping
rm -rf $SAMPLELOC/$SAMPLE/$SRA





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
#### Obtain/prepare genome data
Foo bar baz

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
