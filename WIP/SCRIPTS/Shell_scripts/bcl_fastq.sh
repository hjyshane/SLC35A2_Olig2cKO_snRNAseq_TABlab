#$ -q all.q
#$ -j y
#$ -cwd 
#$ -S /bin/bash
#$ -terse
#$ -pe smp 16
#$ -N bcl_convert
#$ -m beas                   # Email at the beginning and end of the job
#$ -M june.yoon@nationwidechildrens.org

/igm/apps/bclconvert/bcl-convert-4.1.7/bcl-convert --bcl-input-directory /igm/runs/IGM_Seq05/250721_A01461_0814_AH32TKDRX7 \
--output-directory /igm/projects/250718_GSL-EM-4483/FASTQ_out \
--sample-sheet /igm/projects/250718_GSL-EM-4483/SampleSheet.csv