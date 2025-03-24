# STARTING DIRECTORY: mirza@graham.computecanada.ca:/home/mirza/scratch



# 1. CREATE FOLDER STRUCTURE AND EXPLORE DATA FILES

# Make and access folder called Genomics_Assignment_3 
mkdir Genomics_Assignment_3
cd Genomics_Assignment_3

# Copy the Assignment3_W25_files from laptop to this Genomics_Assignment_3 folder using Scp - once this is done, cd into folder
cd Assignment3_W25_files

# Within this directory, make appropriate sub-directories to store corresponding files
mkdir analysis
cd analysis
mkdir FastQC_out # This will store fastqc files downstream
cd ..

mkdir data

# Go into the 'data' directory - this is where all my data files will be stored
cd data

# Take a peak at the shortreads files using less and ls
ls -lh 
less 136x2_R1.fastq
less 136x2_R2.fastq 
less 136x2_longreads.fastq



# 2. FASTQC & TRIMMING

# Conduct fastqc on our data files so that we can gauge quality score metrics - also, start storing commands and output into analysis_log.txt
module load fastqc
fastqc 136x2_R1.fastq 136x2_R2.fastq -o ../analysis/FastQC_out >
../analysis/analysis_log.txt

# Trim fastq files
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \
136x2_R1.fastq 136x2_R2.fastq \
136x2_R1-pe.fastq /dev/null \
136x2_R2-pe.fastq /dev/null \
ILLUMINACLIP:/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 \
SLIDINGWINDOW:4:15 MINLEN:50 |& tee -a ../analysis/analysis_log.txt

# conduct fastqc on trimmed files in order to re-gauge quality score metrics
fastqc 136x2_R1-pe.fastq 136x2_R2-pe.fastq -o ../analysis/FastQC_out >> ../analysis/analysis_log.txt

# Delete the zipped files within the FastQC_out folder
cd ..
cd analysis
cd FastQC_out
rm ./*.zip

# Take a look at the HTML files now for pre trimmed and post trimmed


# 3. SPADES ASSEMBLY

# Make a new folder in analysis which will hold genome assembly info
cd ..
mkdir spades_genome_assembly

cd ..
cd data

nano run_spades.sh # Create a .sh file that will run the spades genome assembly

    #!/bin/bash
    #SBATCH --output=spades_output.log
    #SBATCH --time=4:00:00
    #SBATCH --cpus-per-task=16
    #SBATCH --mem=16G
    #SBATCH --account=def-skremer
    #SBATCH --mail-type=END,FAIL
    #SBATCH --mail-user=mahmad15@uoguelph.ca

    module load spades

    spades.py -1 ./136x2_R1_pe.fastq \
            -2 ./136x2_R2_pe.fastq \
            -o ../analysis/spades_genome_assembly >> ../analysis/analysis_log.txt



# 4. QUAST GENOME COMPARISON

#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --account=def-skremer

# Load the following modules for Quast
module load StdEnv/2020
module load gcc/9.3.0
module load python
module load quast

# Use Quast to retrieve and compare the genome assemblies between SPADES, FLYE and UNICYCLER assemblies
quast -o ../analysis/quast_output/ ../analysis/spades_genome_assembly/contigs.fasta \
./Flye_136x2/30-contigger/contigs.fasta \
./Unicycler_136x2/Unicycler_assembly.fasta >> ../analysis/analysis_log.txt



# 5. ABRICATE GENOME COMPARISON

#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --account=def-skremer

# Load the following modules for Abricate
module load StdEnv/2020
module load gcc/9.3.0
module load python
module load abricate

# Use abricate on all of the genome assemblies

# Create naming conventions
Spades_assembly="/home/mirza/scratch/Genomics_Assignment_3/Assignment3_W25_files/analysis/spades_genome_assembly/scaffolds.fasta"
Flye_assembly="/home/mirza/scratch/Genomics_Assignment_3/Assignment3_W25_files/data/Flye_136x2/Flye_assembly.fasta"
Unicycler_assembly="/home/mirza/scratch/Genomics_Assignment_3/Assignment3_W25_files/data/Unicycler_136x2/Unicycler_assembly.fasta"

# Make output directory to hold analysis results
mkdir -p ../analysis/abricate_analysis_outputs
output_dir=../analysis/abricate_analysis_outputs

# Make subdirectories for each assembly
mkdir -p $output_dir/spades_abricate_output
mkdir -p $output_dir/flye_abricate_output
mkdir -p $output_dir/unicycler_abricate_output


# Run Abricate for SPADES
abricate --db plasmidfinder $Spades_assembly > $output_dir/spades_abricate_output/Spades_plasmidfinder_output.txt 2>> ../analysis/analysis_log.txt
abricate --db resfinder $Spades_assembly > $output_dir/spades_abricate_output/Spades_resfinder_output.txt 2>> ../analysis/analysis_log.txt
abricate --db vfdb $Spades_assembly > $output_dir/spades_abricate_output/Spades_vfdb_output.txt 2>> ../analysis/analysis_log.txt

# Run Abricate for FLYE
abricate --db plasmidfinder $Flye_assembly > $output_dir/flye_abricate_output/Flye_plasmidfinder_output.txt 2>> ../analysis/analysis_log.txt
abricate --db resfinder $Flye_assembly > $output_dir/flye_abricate_output/Flye_resfinder_output.txt 2>> ../analysis/analysis_log.txt
abricate --db vfdb $Flye_assembly > $output_dir/flye_abricate_output/Flye_vfdb_output.txt 2>> ../analysis/analysis_log.txt

# Run Abricate for UNICYCLER
abricate --db plasmidfinder $Unicycler_assembly > $output_dir/unicycler_abricate_output/Unicycler_plasmidfinder_output.txt 2>> ../analysis/analysis_log.txt
abricate --db resfinder $Unicycler_assembly > $output_dir/unicycler_abricate_output/Unicycler_resfinder_output.txt 2>> ../analysis/analysis_log.txt
abricate --db vfdb $Unicycler_assembly > $output_dir/unicycler_abricate_output/Unicycler_vfdb_output.txt 2>> ../analysis/analysis_log.txt