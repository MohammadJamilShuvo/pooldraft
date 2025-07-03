#!/bin/bash
#SBATCH --job-name=collembola_pipeline
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=3-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@uni-freiburg.de
#SBATCH --output=11_logs/%x_%j.out
#SBATCH --error=11_logs/%x_%j.err

cd /pfs/work9/workspace/scratch/fr_ms2252-collembola/entomobryo_project

LOGDIR="11_logs"
mkdir -p "$LOGDIR"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOGFILE="$LOGDIR/full_pipeline_$TIMESTAMP.log"

bash ./run_colenv.sh 2>&1 | tee "$LOGFILE"
