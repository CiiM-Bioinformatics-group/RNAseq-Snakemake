
for file in RA*_R*.fastq.gz; do
  sampleid=$(echo "$file" | sed 's/^\(RA[0-9]*\)-[0-9]*_.*/\1/')
  mkdir -p "$sampleid"
  mv "$file" "$sampleid"
done