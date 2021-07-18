### Mapping sequence reads to a plastid genome
To measure sequencing depth and evenness of a plastid genome, the sequence reads that were originally used to assemble that genome must be mapped against its complete genome sequence. This page explains the process of mapping sequence reads against a given plastid genome. For the purpose of this tutorial, we will select a plastid genome record from GenBank as input.

---

## Prerequisites

Information needed:
+ GenBank accession number of plastid genome record (e.g., NC_026562)
+ Corresponding NCBI SRA number for raw reads (e.g., )

Software needed:
+ ncbi-entrez-direct (version number here)
+ sra-toolkit (>2.10.8)
+ samtools (>1.10)
+ bowtie2 (>2.3.4.1)
+ trimmomatic (version number here)

---

## Step-by-step workflow

##### Define target genome record and corresponding sequence reads
```bash
ACCESSION=NC_026562
SRA_NUMBER=$3
```

##### If software not installed globally, setting paths
```bash
PATH_ENTREZ=/PATH/TO/NCBI-ENTREZ-DIRECT/
PATH_SRATOOL=/PATH/TO/SRA-TOOLKIT/
PATH_TRIMMO=/PATH/TO/TRIMMOMATIC/
PATH_BOWTIE2=/PATH/TO/BOWTIE2/
PATH_SAMT=/PATH/TO/SAMTOOLS/
```

##### Download genome record as flatfile and FASTA file from NCBI GenBank
```bash
${PATH_ENTREZ}esearch -db nucleotide -query $ACCESSION | \
  ${PATH_ENTREZ}efetch -format gb > ${ACCESSION}.gb

${PATH_ENTREZ}esearch -db nucleotide -query $ACCESSION | \
  ${PATH_ENTREZ}efetch -format fasta > ${ACCESSION}.fasta
```

##### Download corresponding sequence reads from NCBI SRA
```bash
mkdir -p ReadMapping
${PATH_SRATOOL}prefetch --max-size 50000000 $SRA_NUMBER
${PATH_SRATOOL}fasterq-dump.2.10.8 --split-3 \
  --skip-technical $SRA_NUMBER/$SRA_NUMBER.sra -O ReadMapping
```

##### Trim raw sequence reads via Trimmomatic
```bash
java -jar ${PATH_TRIMMO}/trimmomatic-0.39.jar PE \
  ReadMapping/${SRA_NUMBER}.sra_1.fastq \
  ReadMapping/${SRA_NUMBER}.sra_2.fastq  \
  ReadMapping/${SRA_NUMBER}_1_PE.fastq \
  ReadMapping/${SRA_NUMBER}_1_SE.fastq \
  ReadMapping/${SRA_NUMBER}_2_PE.fastq \
  ReadMapping/${SRA_NUMBER}_2_SE.fastq \
  ILLUMINACLIP:/${PATH_TRIMMO}/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads \
  LEADING:3 TRAILING:3 MINLEN:36
```

##### Map reads via Bowtie2
```bash
mkdir -p ReadMapping/db

${PATH_BOWTIE2}bowtie2-build ${ACCESSION}.fasta ReadMapping/db/$ACCESSION

${PATH_BOWTIE2}bowtie2 -x ReadMapping/db/$ACCESSION \
  -1 ReadMapping/${SRA_NUMBER}_1_PE.fastq \
  -2 ReadMapping/${SRA_NUMBER}_2_PE.fastq \
  -S ReadMapping/${ACCESSION}_mapping.sam
```


##### Map reads via Bowtie2
```bash
${PATH_SAMT}samtools view -Sb -F 0x04 ReadMapping/${ACCESSION}_mapping.sam \
  > ReadMapping/${ACCESSION}_mapping_OneMoreLocations.bam

${PATH_SAMT}samtools sort ReadMapping/${ACCESSION}_mapping_OneMoreLocations.bam \
  > ${ACCESSION}_mapping_OneMoreLocations.sorted.bam

${PATH_SAMT}samtools index ${ACCESSION}_mapping_OneMoreLocations.sorted.bam
```
