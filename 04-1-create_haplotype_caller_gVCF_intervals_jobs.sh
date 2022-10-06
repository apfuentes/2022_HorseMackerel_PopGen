#!/bin/bash

# Receive the file path from user input in the Terminal window (first position after the script call).
INTERVAL_FILES_LIST="$1"

while read -r line; do

echo -e "#!/bin/bash
" > job_HOM.$line.sh
echo -e "#SBATCH -A snic2019-3-66" >> job_HOM.$line.sh
echo -e "#SBATCH -M snowy" >> job_HOM.$line.sh
echo -e "#SBATCH -p core -n 6" >> job_HOM.$line.sh
echo -e "#SBATCH -t 10-00:00:00" >> job_HOM.$line.sh
echo -e "#SBATCH -J $line.HOM_gVCF" >> job_HOM.$line.sh
echo -e "#SBATCH -e $line.HOM_gVCF_%J_%A_%a.err" >> job_HOM.$line.sh
echo -e "#SBATCH -o $line.HOM_gVCF_%J_%A_%a.out" >> job_HOM.$line.sh
echo -e "#SBATCH --mail-type=all" >> job_HOM.$line.sh
echo -e "#SBATCH --mail-user=angela.fuentespardo@gmail.com" >> job_HOM.$line.sh
echo -e "#SBATCH -a 1-12
" >> job_HOM.$line.sh

echo -e "# Load required software." >> job_HOM.$line.sh
echo -e "module load bioinfo-tools" >> job_HOM.$line.sh
echo -e "module load GATK/3.8-0
" >> job_HOM.$line.sh

echo -e "# Directory where the gVCF files will be stored." >> job_HOM.$line.sh
echo -e "cd /proj/uppstore2017191/nobackup/private/UserDirectories/angela/HorseMackerel/data/04-gVCF-files
" >> job_HOM.$line.sh

echo -e "# Set path to required files and directories as environment variables." >> job_HOM.$line.sh
echo -e "echo This is array job: \$SLURM_ARRAY_TASK_ID
" >> job_HOM.$line.sh

echo -e "BAM_FILE='./bam_files_HC_gVCF.txt'" >> job_HOM.$line.sh
echo -e "ID_FILE='./sample_IDs_HC_gVCF.txt'
" >> job_HOM.$line.sh

echo -e "BAM_target=\$(sed -n \"\$SLURM_ARRAY_TASK_ID\"p \$BAM_FILE)" >> job_HOM.$line.sh
echo -e "ID_target=\$(sed -n \"\$SLURM_ARRAY_TASK_ID\"p \$ID_FILE)
" >> job_HOM.$line.sh

echo -e "echo -e This is the BAM target: \${BAM_target}" >> job_HOM.$line.sh
echo -e "echo -e This is the ID target: \${ID_target}
" >> job_HOM.$line.sh

echo -e "# Create sample subdirectory if it does not exist already." >> job_HOM.$line.sh
echo -e "if [ -d \${ID_target} ]; then echo \"\${ID_target} dir exists\"; else mkdir \${ID_target}; fi" >> job_HOM.$line.sh
echo -e "# Create directory for temporal files." >> job_HOM.$line.sh
echo -e "if [ -d \$SNIC_TMP/tmp ]; then echo \"tmp/ exists in SNIC_TMP\"; else mkdir \$SNIC_TMP/tmp; fi
" >> job_HOM.$line.sh

echo -e "# Run GATK HaplotypeCaller to generate a single gVCF file for each pool sample and scaffold interval." >> job_HOM.$line.sh
echo -e "java -Djava.io.tmpdir=\$SNIC_TMP/tmp -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx32G -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar \\" >> job_HOM.$line.sh
echo -e "-T HaplotypeCaller \\" >> job_HOM.$line.sh
echo -e "-R /proj/uppstore2017191/nobackup/private/UserDirectories/angela/HorseMackerel/data/00-genome/v.1.0/fTraTra1.PB.asm1.purge1.scaff1.fa \\" >> job_HOM.$line.sh
echo -e "-I \${BAM_target} \\" >> job_HOM.$line.sh
echo -e "--emitRefConfidence GVCF \\" >> job_HOM.$line.sh
echo -e "-mbq 20 \\" >> job_HOM.$line.sh
echo -e "-minPruning 5 \\" >> job_HOM.$line.sh
echo -e "-L $line.list \\" >> job_HOM.$line.sh
echo -e "-o ./\${ID_target}/\${ID_target}.$line.g.vcf
" >> job_HOM.$line.sh

sbatch job_HOM.$line.sh

done <"$INTERVAL_FILES_LIST"