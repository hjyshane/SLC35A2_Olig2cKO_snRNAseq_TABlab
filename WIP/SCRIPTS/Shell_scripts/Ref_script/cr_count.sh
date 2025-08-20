#$ -j y
#$ -cwd
#$ -S /bin/bash
#$ -q all.q
#$ -pe smp 16
#$ -V

sample=24-2687_JacksonLabs_Pg_Spleen
fastq_path=/igm/efs-globus/globus/Chris_Lauber/241016_Lauber_GSL-AB-4127
/igm/apps/10X_chromium/cellranger-8.0.0/bin/cellranger count --id=$sample \
                                                                 --transcriptome=/igm/apps/10X_chromium/refdata-gex-GRCm39-2024-A \
                                                                 --fastqs=${fastq_path} \
                                                                 --sample=$sample \
                                                                 --create-bam=true \
								 --localcores=16 \
                                                                 --localmem=64 \
