# RNA-Seq Analysis Pipeline

This project provides a robust and flexible pipeline for processing RNA sequencing data—from raw reads to differential gene expression results. It is designed for ease of use, reproducibility, and customization in transcriptomic analyses.

---

## Pipeline Workflow

### 1. Data Retrieval
- **Download RNA-Seq datasets:** Access data from the [Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra).
- **Reference Files:** Obtain genome indices and gene annotation files (GTF) for downstream analysis.

### 2. Data Pre-processing
- **Format Conversion:** Convert SRA files to FASTQ format.
- **Quality Assessment:** Evaluate raw read quality using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and summarize reports with [MultiQC](https://multiqc.info/).
- **Trimming:** Remove adapters and low-quality bases using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic).

### 3. Alignment and Post-processing
- **Mapping:** Align reads to a reference genome using [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml).
- **File Processing:** Sort and index BAM files using [Samtools](http://www.htslib.org/).

### 4. Quantification and Analysis
- **Counting Reads:** Use [featureCounts](http://bioinf.wehi.edu.au/featureCounts/) to quantify gene-level read counts.
- **Differential Expression:** Perform analysis with [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) in [R](https://www.r-project.org/).
- **Visualization:** Generate plots such as volcano plots, MA plots, PCA, and heatmaps.

---

## Key Features

- **Modular Design:** Easily update or swap components of the pipeline.
- **Automation:** End-to-end processing from data download to analysis results.
- **Scalability:** Adaptable for both small-scale projects and large studies.
- **Reproducibility:** Ensured by version control and configuration files.

---

## System & Software Requirements

- **Operating System:** Linux (Ubuntu recommended)
- **Languages:** Bash and R
- **Essential Tools:**
  - [Git](https://git-scm.com)
  - [SRA Toolkit](https://github.com/ncbi/sra-tools)
  - [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
  - [MultiQC](https://multiqc.info/)
  - [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
  - [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml)
  - [Samtools](http://www.htslib.org/)
  - [featureCounts](http://bioinf.wehi.edu.au/featureCounts/)
  - [R](https://www.r-project.org/) with [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)

---

## Repository Structure

```plaintext
├── bin
│   └── DE_seq.R
├── Data
│   ├── data.sh
│   └── Readme.pdf
├── Docs
│   └── toolsandsteps.pdf
├── Environment
│   ├── environment.yml
│   └── Readme.pdf
├── Fastqc Analysis
│   ├── controlled1_trimmed_fastqc.pdf
│   ├── controlled_trimmed_fastqc.pdf
│   ├── infected1_trimmed_fastqc.pdf
│   ├── infected_trimmed_fastqc.pdf
│   ├── SRR11412217_pass_fastqc.pdf
│   ├── SRR11412218_pass_fastqc.pdf
│   ├── SRR11412227_pass_fastqc.pdf
│   └── SRR11412228_pass_fastqc.pdf
├── Pipeline
│   └── rna_seq.sh
└── Results
    ├── DE_results.csv
    ├── Ma plot.png
    ├── pca plot.png
    └── volcano plot.png
```
---
# 1. Clone the Repository:
git clone https://github.com/yourusername/RNASeq-Analysis-Pipeline.git
cd RNASeq-Analysis-Pipeline

# 2. Install Dependencies:
# Ensure all required tools and R packages are installed.
# Refer to the System & Software Requirements section above.

# 3. Run the Pipeline:

# Data Retrieval:
# Execute the 'data.sh' script in the Data directory to download RNA-Seq datasets,
# reference genome, and annotation files.
bash Data/data.sh

# Processing:
# Run the 'rna_seq.sh' script in the Pipeline directory to perform quality control,
# alignment, and quantification.
bash Pipeline/rna_seq.sh

# Analysis:
# Use the R script 'DE_seq.R' in the bin directory (or follow the in-repository documentation)
# to perform differential expression analysis and generate visual reports.
Rscript bin/DE_seq.R
