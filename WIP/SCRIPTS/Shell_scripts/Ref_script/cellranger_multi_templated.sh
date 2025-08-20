
#$ -j y
#$ -cwd
#$ -S /bin/bash
#$ -N FF006_multi
#$ -q all.q
#$ -pe smp 16
#$ -V

/igm/apps/10X_chromium/cellranger-8.0.0/bin/cellranger multi --id=FF006 \
                                                                 --csv=FF006_multi_template_multiplex.csv \
                                                                 --localcores=16 \
                                                                 --localmem=64 \
								 