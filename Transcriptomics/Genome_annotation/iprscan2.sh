#!/bin/bash
#################################################################################################################
#SBATCH --job-name=sbatchTemplate ## Name of your job
#SBATCH --ntasks=64 ## number of cpu's to allocate for a job
#SBATCH --ntasks-per-node=64 ## number of cpu's to allocate per each node
#SBATCH --nodes=1 ## number of nodes to allocate for a job
#SBATCH --mem=256G ## memory to allocate for your job in MB
#SBATCH --time=1-00:00:00 ## time to allocate for your job in format: DD-HH:MM:SS
#SBATCH --error=%J.errors ## stderr file name(The %J will print job ID number)
#SBATCH --output=%J.output ## stdout file name(The %J will print job ID number)
#SBATCH --mail-type=NONE ## Send your job status via e-mail: Valid type values are NONE, BEGIN, END, FAIL, REQUEUE, ALL
########### Job information #############
echo "================================"
echo "Start at `date`"
echo "Job id is $SLURM_JOBID"
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NTASKS processors."
echo "================================"
#########################################

######## Load required modules ##########
#. /etc/profile.d/modules.sh # Required line for modules environment to work
#module load openmpi/1.8.4 python/2.7 # Load modules that are required by your program
#conda init
#conda activate multiqc
#source /lustre1/home/mass/eskalon/miniconda/bin/activate agat
source /lustre1/home/mass/eskalon/miniconda/bin/activate funannotate
#########################################

### Below you can enter your program job command ###

#agat_convert_sp_gff2gtf.pl \
#-gff Stylophora_pistillata_gca002571385v1.Stylophora_pistillata_v1.60.gff3 \
#-o Stylophora_pistillata.ensembl.gtf


# input
FAA=protein.faa

# outputs
OUTDIR=/lustre1/home/mass/eskalon/Oculina/genomes/annot/interproscan_ocupat_v1.1
TMPDIR=/lustre1/home/mass/eskalon/Oculina/genomes/tmp/interproscan_spis_v1.1

mkdir -p "$OUTDIR" "$TMPDIR"

/lustre1/home/mass/eskalon/interpro/interproscan/interproscan.sh \
  -i "$FAA" \
  -appl Pfam,PANTHER,NCBIFAM,SUPERFAMILY,CDD \
  -d "$OUTDIR" \
  -T "$TMPDIR" \
  -f TSV,XML \
  -cpu 64 \
  -goterms \
  -pa \
  --iprlookup
