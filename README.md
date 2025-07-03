# Collembola Genome Assembly and Annotation Pipeline

This repository contains a reproducible pipeline for reference genome assembly and functional annotation of a Collembola species using pooled sequencing (pool-seq) data.

## 🧬 Project Overview

- **Species**: Entomobryo sp. (Collembola)
- **Input**: Paired-end Illumina FASTQ files
- **Tools**: Conda, FastQC, fastp, MEGAHIT, SPAdes, QUAST, BUSCO, RepeatModeler, RepeatMasker, BRAKER2, InterProScan, eggNOG-mapper
- **Platform**: bwUniCluster 3.0 (HPC environment)
- **Workflow Management**: Bash + SLURM

## 📁 Directory Structure

```
entomobryo_project/
├── 00_tools/             # Installed tools and environments
├── 01_raw_data/          # Raw FASTQ files
├── 02_qc/                # Quality control reports
├── 03_trimmed_reads/     # Trimmed reads
├── 04_kmer_analysis/     # Jellyfish + GenomeScope
├── 05_assembly/          # Assembly results (MEGAHIT, SPAdes)
├── 06_polishing/         # Polishing (optional)
├── 07_quality_check/     # QUAST + BUSCO
├── 08_repeats/           # Repeat masking
├── 09_annotation/        # BRAKER2 gene models
├── 10_functional/        # InterProScan + eggNOG-mapper
├── 11_logs/              # All logs
├── genome_pipeline.sh    # Full pipeline script
└── submit_pipeline.sh    # SLURM submission wrapper
```

## 🚀 How to Run (On HPC)

1. Clone this repo to your HPC workspace.
2. Upload your raw data to `01_raw_data/`.
3. Edit email and resource settings in `submit_pipeline.sh`.
4. Make scripts executable:
    ```bash
    chmod +x genome_pipeline.sh
    chmod +x submit_pipeline.sh
    ```
5. Submit to SLURM:
    ```bash
    sbatch submit_pipeline.sh
    ```

## 📦 Outputs

- Full genome assemblies per sample
- Gene predictions (GFF/GTF)
- Functional annotations (InterPro, eggNOG)
- Summary stats and BUSCO completeness

## 👤 Author

Mohammad Jamil Shuvo – jamilshuvo94@gmail.com
