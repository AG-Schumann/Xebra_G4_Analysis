#!/bin/bash

#SBATCH --job-name=XeBRA_MC
#SBATCH --output=job-prompt.txt

#SBATCH --ntasks=1
#SBATCH --array=1-10

#source ~/.bashrc
cd
#loadG4
user=$(whoami)
mkdir /scratch/$user -p
tmp=/scratch/$user/"$2_$3_noscint_${SLURM_ARRAY_TASK_ID}.root"
dst=$1/"$2_$3_noscint_${SLURM_ARRAY_TASK_ID}.root"
source /opt/geant/v10.3.3/bin/geant4.sh && source /opt/geant/v10.3.3/share/Geant4-10.3.3/geant4make/geant4make.sh && export G4WORKDIR=.
~/notebooks/Xebra_G4/bin/Linux-g++/xebra_G4 -p ~/notebooks/Xebra_G4/macros/preinit.mac -f $1/"src_Pointsources_DP.mac" -n $3 -o $tmp
mv $tmp $dst

