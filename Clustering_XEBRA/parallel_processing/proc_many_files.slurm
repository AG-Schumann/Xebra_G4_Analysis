#!/bin/bash
##SBATCH --output=/home/atp/ybiondi/slurmlog
##SBATCH --error=/home/atp/ybiondi/slurmlog
#SBATCH --nodes=2

echo "slurm host is"
echo $SLURMD_NODENAME
echo "slurm queue is"
echo $SLURM_JOB_PARTITION
echo "slurm working directory absolute is"
echo $SLURM_SUBMIT_DIR
echo "slurm batch id"
echo $SLURM_JOB_ID
echo "slurm job name is"
echo $SLURM_JOB_NAME

source $HOME/.bashrc
/bin/echo "my username is " $USER

mkdir -p /scratch/$USER/$SLURM_JOB_ID
cd /scratch/$USER/$SLURM_JOB_ID

echo "$binary_file -d "$work_dir" -i "$inputfile" "
$binary_file -d "$work_dir" -i "$inputfile"

cp -p /scratch/$USER/$SLURM_JOB_ID/"$inputfile" /disk/bulk_atp/ybiondi/darwin/
cd ../
rm -rf /scratch/$USER/$SLURM_JOB_ID

wait

#sleep 300
