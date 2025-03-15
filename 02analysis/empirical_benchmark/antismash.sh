#!/bin/bash
#SBATCH --job-name=antismash_array
#SBATCH --output=antismash_%A_%a.out
#SBATCH --error=antismash_%A_%a.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --array=0-5
#SBATCH --time=24:00:00
#SBATCH --partition=all

# Load Conda (Make sure your system supports this)
source /miniconda3/etc/profile.d/conda.sh  # Adjust path if needed
conda activate antismash

FILES=(
    "/results/assembly-annotation-eval/OAK003.A0101.raw-raw.proteins.carpedeam.configSafe0/OAK003.A0101.raw-raw.proteins.carpedeam.configSafe0.fna"
    "/results/assembly-annotation-eval/OAK003.A0101.raw-raw.proteins.carpedeam.configUnsafe0/OAK003.A0101.raw-raw.proteins.carpedeam.configUnsafe0.fna"
    "/results/assembly-annotation-eval/OAK003.A0101.raw-raw.proteins.megahit.config10/OAK003.A0101.raw-raw.proteins.megahit.config10.fna"
    "/results/assembly-annotation-eval/OAK003.A0101.raw-raw.proteins.penguin.config10/OAK003.A0101.raw-raw.proteins.penguin.config10.fna"
    "/results/assembly-annotation-eval/OAK003.A0101.raw-raw.proteins.spades.config10/OAK003.A0101.raw-raw.proteins.spades.config10.fna"
)

# Get the filename based on SLURM_ARRAY_TASK_ID
INPUT_FILE=${FILES[$SLURM_ARRAY_TASK_ID]}

# Extract the basename for output naming
BASENAME=$(basename "$INPUT_FILE" .fna)

# Define output directory
OUTPUT_DIR="/antismash/processAll/$BASENAME"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Run antiSMASH
antismash -c 16 \
    --databases /results/antismashdb \
    --genefinding-tool prodigal-m \
    --output-dir "$OUTPUT_DIR" \
    "$INPUT_FILE"

