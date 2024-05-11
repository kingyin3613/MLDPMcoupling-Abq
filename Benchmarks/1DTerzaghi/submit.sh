#!/bin/bash
#SBATCH --account=p12345
#SBATCH --partition=normal
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=8G
#SBATCH --job-name=abaqus_coupled
#SBATCH --output=abaqus_%j.out
#SBATCH --error=abaqus_%j.err

# Change directory to your working directory
cd /projects/p12345/your/abaqus/working/directory
# Load Abaqus module (replace "abaqus/2020" with "abaqus" if use abaqus 6.13)
module purge
module load abaqus/2020

#------------------------------------------------------------------------------------
# Run Abaqus jobs (always create pipes first, then run VUEL, finally run UEL)
rm *.lck
rm LDPM2FLM.pipe
rm FLM2LDPM.pipe
mkfifo LDPM2FLM.pipe
mkfifo FLM2LDPM.pipe

abaqus interactive analysis job=prism_100mmx500mm input=prism_100mmx500mm.inp USER=VUEL_MLDPM.for double=both ask_delete=OFF &
PID1=$!
echo $PID1

abaqus interactive analysis job=prism_100mmx500mm-edgeEle input=prism_100mmx500mm-edgeEle.inp USER=UEL_FLM.for ask_delete=OFF
PID2=$!
echo $PID2

wait $PID1
wait $PID2

#------------------------------------------------------------------------------------
# Always delete the pipe files once a two-way coupling analysis is finished.


#------------------------------------------------------------------------------------
# In case of any job is aborted, check if there exists a .lck file, if one exists, run the corresponding Abaqus termination command as follows, then delete the .lck and .pipe files
# abaqus job=prism_100mmx500mm terminate
# abaqus job=prism_100mmx500mm-edgeEle terminate
