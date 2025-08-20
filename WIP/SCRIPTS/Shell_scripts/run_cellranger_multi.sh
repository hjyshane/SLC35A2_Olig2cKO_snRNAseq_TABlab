#$ -q all.q
#$ -j y
#$ -cwd 
#$ -S /bin/bash
#$ -terse
#$ -pe smp 16
#$ -N cellranger_mkfastq
#$ -m beas                   # Email at the beginning and end of the job
#$ -M june.yoon@nationwidechildrens.org

module load cellranger_9.0.0

cellranger multi \
  --id=FF064 \
  --csv=/igm/projects/250718_GSL-EM-4483/multi_config_FF064.csv \
  --localcores=16 \
  --localmem=64