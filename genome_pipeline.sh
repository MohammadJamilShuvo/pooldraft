#!/bin/bash

# ===============================
# Collembola Reference Genome Pipeline
# ===============================

# === Activate conda environment ===
source /pfs/work9/workspace/scratch/fr_ms2252-collembola/00_tools/miniconda3/etc/profile.d/conda.sh
conda activate genome_env

# === Variables ===
WORKDIR="/pfs/work9/workspace/scratch/fr_ms2252-collembola/entomobryo_project"
THREADS=16
KMER_SIZE=21
MIN_CONTIG_LEN=1000
LINEAGE="arthropoda_odb10"

# === Step 0: Setup Environment (SKIPPED) ===
# cd "$TOOLDIR"
# wget ...
# bash miniconda.sh ...
# conda create ...
# conda install ...

# === Step 1: QC ===
fastqc "$WORKDIR/01_raw_data/"*.fastq.gz -o "$WORKDIR/02_qc"
cd "$WORKDIR/02_qc"
multiqc . -o .

# === Step 2: Trimming ===
cd "$WORKDIR/03_trimmed_reads"
for R1 in "$WORKDIR/01_raw_data/"*_R1_001.fastq.gz; do
  SAMPLE=$(basename "$R1" _R1_001.fastq.gz)
  R2="$WORKDIR/01_raw_data/${SAMPLE}_R2_001.fastq.gz"
  if [[ -f "$R2" ]]; then
    fastp -i "$R1" -I "$R2" -o "${SAMPLE}_R1_trimmed.fastq.gz" -O "${SAMPLE}_R2_trimmed.fastq.gz" \
          -h "${SAMPLE}_fastp.html" -j "${SAMPLE}_fastp.json"
  fi
done

# === Step 3: K-mer analysis ===
cd "$WORKDIR/04_kmer_analysis"
zcat "$WORKDIR"/03_trimmed_reads/*_R1_trimmed.fastq.gz "$WORKDIR"/03_trimmed_reads/*_R2_trimmed.fastq.gz > all_trimmed_reads.fastq
jellyfish count -C -m $KMER_SIZE -s 100M -t $THREADS all_trimmed_reads.fastq -o reads.jf
jellyfish histo -t $THREADS reads.jf > kmer_histogram.txt

# === Step 4: Assembly (MEGAHIT & SPAdes) ===
cd "$WORKDIR/05_assembly"

R1_A="$WORKDIR/03_trimmed_reads/2317_S1_L002_R1_trimmed.fastq.gz"
R2_A="$WORKDIR/03_trimmed_reads/2317_S1_L002_R2_trimmed.fastq.gz"
megahit -1 "$R1_A" -2 "$R2_A" -o megahit_plotA --presets meta-sensitive --min-contig-len $MIN_CONTIG_LEN -t $THREADS
spades.py -1 "$R1_A" -2 "$R2_A" -o spades_plotA --careful -t $THREADS -m 64

R1_B="$WORKDIR/03_trimmed_reads/2324_S2_L002_R1_trimmed.fastq.gz"
R2_B="$WORKDIR/03_trimmed_reads/2324_S2_L002_R2_trimmed.fastq.gz"
megahit -1 "$R1_B" -2 "$R2_B" -o megahit_plotB --presets meta-sensitive --min-contig-len $MIN_CONTIG_LEN -t $THREADS
spades.py -1 "$R1_B" -2 "$R2_B" -o spades_plotB --careful -t $THREADS -m 64

# === Step 5: QUAST + BUSCO ===
cd "$WORKDIR/07_quality_check"
quast.py "$WORKDIR/05_assembly/megahit_plotA/final.contigs.fa" "$WORKDIR/05_assembly/megahit_plotB/final.contigs.fa" \
         "$WORKDIR/05_assembly/spades_plotA/scaffolds.fasta" "$WORKDIR/05_assembly/spades_plotB/scaffolds.fasta" \
         -o quast_results -t $THREADS --labels megahit_A,megahit_B,spades_A,spades_B

for A in "$WORKDIR"/05_assembly/*/*.{fa,fasta}; do
  [[ -f "$A" ]] && busco -i "$A" -o "busco_$(basename $(dirname "$A"))" -l $LINEAGE -m genome -c $THREADS
done

# === Step 6: Repeat Masking ===
cd "$WORKDIR/08_repeats"
for A in "$WORKDIR"/05_assembly/*/*.{fa,fasta}; do
  SAMPLE=$(basename $(dirname "$A"))
  cp "$A" "${SAMPLE}.fa"
  BuildDatabase -name "${SAMPLE}_rmdb" "${SAMPLE}.fa"
  RepeatModeler -database "${SAMPLE}_rmdb" -pa $THREADS -LTRStruct -dir "./${SAMPLE}_modeler"
  RepeatMasker -pa $THREADS -lib "${SAMPLE}_modeler/consensi.fa.classified" -dir "./${SAMPLE}_masked" -xsmall "${SAMPLE}.fa"
done

# === Step 7: BRAKER2 Gene Prediction ===
cd "$WORKDIR/09_annotation"
for MASKED in "$WORKDIR"/08_repeats/*.fa.masked; do
  SAMPLE=$(basename "$MASKED" .fa.masked)
  braker.pl --genome="$MASKED" --species="$SAMPLE" --softmasking --cores $THREADS \
            --workingdir "$WORKDIR/09_annotation/braker_${SAMPLE}" --skipAllTraining
done

# === Step 8: Functional Annotation ===
cd "$WORKDIR/10_functional"
for SAMPLE in megahit_plotA megahit_plotB spades_plotA spades_plotB; do
  GTF="$WORKDIR/09_annotation/braker_${SAMPLE}/braker.gtf"
  GENOME="$WORKDIR/08_repeats/${SAMPLE}.fa.masked"
  gffread "$GTF" -g "$GENOME" -x "${SAMPLE}_cds.fasta" -y "${SAMPLE}_prot.fasta"
  interproscan.sh -i "${SAMPLE}_prot.fasta" -f tsv -o "${SAMPLE}_interpro.tsv" -goterms -pa -dp -cpu $THREADS
  emapper.py -i "${SAMPLE}_prot.fasta" --itype proteins -o "${SAMPLE}_eggnog" --cpu $THREADS --output_dir "$WORKDIR/10_functional"
done

# === Step 9: Summary & Archive ===
cd "$WORKDIR"
for SAMPLE in megahit_plotA megahit_plotB spades_plotA spades_plotB; do
  gffread "$WORKDIR/09_annotation/braker_${SAMPLE}/braker.gtf" -T -o "$WORKDIR/09_annotation/${SAMPLE}_braker.gff"
done

tar -czf collembola_results.tar.gz 02_qc 03_trimmed_reads 04_kmer_analysis 05_assembly 07_quality_check 08_repeats 09_annotation 10_functional
tar -czf logs.tar.gz 11_logs
