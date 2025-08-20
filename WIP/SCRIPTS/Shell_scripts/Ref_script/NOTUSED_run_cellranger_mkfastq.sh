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

cellranger mkfastq \
--samplesheet=/igm/runs/IGM_Seq05/250721_A01461_0814_AH32TKDRX7/SampleSheet.csv \
--run=/igm/runs/IGM_Seq05/250721_A01461_0814_AH32TKDRX7 \
--use-bases-mask=Y28,I8,Y91 >& bclconv.out
