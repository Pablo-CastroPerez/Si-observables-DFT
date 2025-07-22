#!/bin/bash
#
# Parallel script produced by bolt
#        Resource: ARCHER2 (HPE Cray EX (128-core per node))
#    Batch system: Slurm
#
#
#SBATCH --job-name=QE_lecture
#SBATCH --nodes=1
#SBATCH --tasks-per-node=36
#SBATCH --exclusive
#SBATCH --time=00:15:0

#SBATCH --account=account_code_here
#SBATCH --partition=standard
#SBATCH --qos=short #standard 


#module load epcc-job-env
export OMP_NUM_THREADS=1

module load quantum_espresso

# Run the parallel program
srun -n 36 pw.x -npool 6 -pd .true. -in si.scf.in > si.scf.out

