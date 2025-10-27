#!/bin/bash
#SBATCH --job-name=bTB_wildlife
#SBATCH --output=logs/bTB_wildlife_%A_%a.out
#SBATCH --error=logs/bTB_wildlife_%A_%a.err
#SBATCH --partition=...
#SBATCH --qos=...
#SBATCH --array=5
#SBATCH --nodes=1
#SBATCH --ntasks=1                # one R job per array element
#SBATCH --cpus-per-task=16        # use 16 cores inside R (doParallel)
#SBATCH --mem=256G                # safely under 16 GB/core × 16 cores
#SBATCH --account=...
#SBATCH --time=24:00:00

echo "Running on host: $(hostname)"
echo "Job started at: $(date)"

module load anaconda
module load gcc/10.3.0
module load gsl/2.7
conda activate btbwl

mkdir -p logs data
cd "$SLURM_SUBMIT_DIR"

echo "Compiling C++ wildlife models..."
g++ -O3 -std=gnu++11 -Wall \
  -I/curc/sw/install/gsl/2.7/gcc/10.3.0/include \
  -L/curc/sw/install/gsl/2.7/gcc/10.3.0/lib \
  bTB_wildlifeModel_CTMC.cpp -lgsl -lgslcblas -lm -o wl_model_CTMC

g++ -O3 -std=gnu++11 -Wall \
  -I/curc/sw/install/gsl/2.7/gcc/10.3.0/include \
  -L/curc/sw/install/gsl/2.7/gcc/10.3.0/lib \
  bTB_wildlifeModel_DTMC.cpp -lgsl -lgslcblas -lm -o wl_model_DTMC

chmod 755 wl_model_CTMC wl_model_DTMC
ls -l wl_model_CTMC wl_model_DTMC
echo "Compilation complete."

# If your cluster has a fast, exec-enabled scratch dir, run from there:
if [[ -n "${SLURM_TMPDIR:-}" ]]; then
  cp wl_model_CTMC wl_model_DTMC "$SLURM_TMPDIR"/
  export WL_MODEL_CTMC="$SLURM_TMPDIR/wl_model_CTMC"
  export WL_MODEL_DTMC="$SLURM_TMPDIR/wl_model_DTMC"
else
  export WL_MODEL_CTMC="$SLURM_SUBMIT_DIR/wl_model_CTMC"
  export WL_MODEL_DTMC="$SLURM_SUBMIT_DIR/wl_model_DTMC"
fi

echo "Launching R job for section $SLURM_ARRAY_TASK_ID..."
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

Rscript --vanilla hpc.R "$SLURM_ARRAY_TASK_ID"

echo "Job finished at: $(date)"
