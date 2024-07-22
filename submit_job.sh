#!/bin/bash
#SBATCH --job-name=YQ-interpolator
#SBATCH --account=qgp
#SBATCH --partition=qgp
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --nodelist=ccc[0224]
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output="/projects/jnorhos/jordi/feyisola/YQ_InterpolatorAndExtrapolator/output/interpolation.out"
#SBATCH --error="/projects/jnorhos/jordi/feyisola/YQ_InterpolatorAndExtrapolator/output/interpolation.err"

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
echo "OpenMP will use ${SLURM_CPUS_PER_TASK} nodes"

./YQ_InterpolatorAndExtrapolator

