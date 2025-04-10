source /swst/apps/anaconda3/2021.05_gcc-8.5.0/etc/profile.d/conda.sh

conda activate synSg

directory="/storage/home/yur97/yur97/synSg/Yiyun_S1"

# Loop through the files in the directory
for file in "$directory"/out*_1.fq; do
    bowtie -n 1 -p 8 -q --norc CABE $file > $file.CABE.align1.txt
    # bowtie -n 1 -p 8 -q --norc CCBE $file > $file.CCBE.align1.txt
done