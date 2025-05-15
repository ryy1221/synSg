source /swst/apps/anaconda3/2021.05_gcc-8.5.0/etc/profile.d/conda.sh

conda activate synSg

directory="/storage/group/epo2/default/yur97/github/synSg/data/sequencing/Library"

# # Loop through the files in the directory
# for file in "$directory"/out*.fq; do
#   if [[ "$file" == *"KA1"* || "$file" == *"KA2"* ]]; then
#     echo $file
#     bowtie -n 1 -p 8 -q --norc CABE $file > $file.align1.txt
#   elif [[ "$file" == *"KC1"* || "$file" == *"KC2"* ]]; then
#     echo $file
#     bowtie -n 1 -p 8 -q --norc CCBE $file > $file.align1.txt
#   fi
# done

# Loop through the files in the directory that start with out.D0_JC2 and end with .fq
for file in "$directory"/out.D0_JC2*.fq; do
    echo $file
    bowtie -n 1 -p 8 -q --norc CCBE "$file" > "$file.align1.txt"
done
