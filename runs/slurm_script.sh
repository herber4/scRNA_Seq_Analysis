
#!/bin/bash
#SBATCH --job-name=cellranger
#SBATCH -n 16
#SBATCH --partition=compute
#SBATCH --time=100:00:00
#SBATCH --mem=128G
#SBATCH --mail-type=all
#SBATCH --output=/austin/pbmc/runs/cellranger.%j.out
#SBATCH --error=/pbmc/runs/cellranger.%j.err
#SBATCH --mail-user=herber4@clemson.edu
# Load change dir

cd /austin/pbmc/runs

cellranger multi --localcores 16 --localmem 128 --disable-ui --id pbmc_slurm_test --csv ../multiconfig.csv
