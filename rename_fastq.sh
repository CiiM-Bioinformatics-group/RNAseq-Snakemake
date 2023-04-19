#for file in /vol/projects/CIIM/SIcohort/RNAseq/raw/RA*fastq.gz; do
  #mv "$file" "$(echo $file | sed 's/_S[0-9]*//g')"
#  mv "$file" "$(echo $file | sed -E 's/(RA[0-9]+-[0-9]+)(_S[0-9]+)(_L[0-9]+)?(_R[0-9]+)(_001)/\1\4/')"
#done

for file in /vol/projects/CIIM/SIcohort/RNAseq/raw/RC*fastq.gz; do
  mv "$file" "$(echo $file | sed -E 's/RC([0-9]+).+(R[0-9]+).+.fastq.gz/RA\1-0_\2.fastq.gz/')"
done

# from
# RA0001-5_S94_R1.fastq.gz
# to 
# RA0001-5_R1.fastq.gz
# and
# RA0062-2_L002_R1_001.fastq.gz
# to
# RA0062-2_R1.fastq.gz
