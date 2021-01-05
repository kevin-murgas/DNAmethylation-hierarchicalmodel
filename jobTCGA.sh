#!/bin/bash
#SBATCH --mem=512
#SBATCH --array=1-87
#SBATCH --error=jParErr-%j.err
#SBATCH --partition=common

export TMPDIR=/dscrhome/kam109/temp
export PATH=/opt/apps/rhel7/R-3.4.4/bin:$PATH
echo "working directory = "$SLURM_SUBMIT_DIR
echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURM_ARRAY_TASK_COUNT="$SLURM_ARRAY_TASK_COUNT
echo "Slurm Array Task ID="${SLURM_ARRAY_TASK_ID}
Rscript StanCParallel.R ${SLURM_ARRAY_TASK_ID} ~/ResultsTCGA/StanCParResults_${SLURM_ARRAY_TASK_ID}.Rdata > ~/ParOut/taskOut${SLURM_ARRAY_TASK_ID}.txt

