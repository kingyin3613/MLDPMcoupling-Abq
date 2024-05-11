In order to run these simulations, you must use Abaqus on UNIX systems, and the Abaqus should be successfully linked to Intel Fortran compiler. We suggest copying the relevant scripts (FORTRAN user subroutine codes `UEL.for`, `VUEL.for`, as well as the geometry and input files) to the working directory.

Here is a simple guide for running two-way coupling on Northwestern Quest:

Option 1: using interactive commands

Login to FastX or PuTTY or other SSH clients with your NetID account
Use the command: 

```bash
srun --account=p12345 --partition=short -N 1 -n 4 --mem=8G --time=04:00:00 --pty bash –l
```

to run an interactive bash session on a single compute node with 4 cores, and access to 8 GB of RAM for up to 4 hours, debited to the p12345 account.

In bash command line session, run the following commands:

```bash
cd /your/abaqus/working/directory on Quest
module load abaqus/2020
mkfifo LDPM2FLM.pipe
mkfifo FLM2LDPM.pipe
abaqus job=LDPMjobname input=LDPMjobname.inp USER=VUEL_MLDPM.for double=both ask_delete=OFF
abaqus job=FLMjobname input=FLMjobname.inp USER=UEL_FLM.for ask_delete=OFF
```

Option 2: using a batch file to submit a coupling job

Edit and upload your batch file (e.g. `submit.sh` attached) on Quest

```bash
cd /your/batchfile/location on Quest
```

Login to FastX or PuTTY or other SSH clients with your NetID account
Use the command: 

```bash
sbatch submit.sh
```