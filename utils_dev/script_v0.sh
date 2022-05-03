#!/bin/bash

# Job name
#SBATCH --job-name=test_1

# Number pf nodes
#SBATCH --nodes=1

# Number of processor per node
#SBATCH --ntasks-per-node=16

# Nomber of RAM per node
#SBATCH --mem=200Go

# Type of machines requested (type1 or 2)
#SBATCH --partition=type_2

# Name of output file
#SBATCH --output=test_1.out

# Calculation times
#SBATCH --time=00-02:00:00

# loading modules
module load userspace/tr17.10
modue load python/3.6.3




########### script ##########
date
echo -e "\nStart calcul cubs\n"

path_folder_pID=(/mnt/beegfs/pturk/PID_couple)
path_folder_data_split=(/mnt/beegfs/pturk/PfamSplit_50)
path_new_folder=(/mnt/beegfs/pturk/test_v1)

cd /trinity/home/pturk
python ./scriptCubs.py $path_folder_pID $path_folder_data_split $path_new_folder
date

echo -e "\nDone\n"
