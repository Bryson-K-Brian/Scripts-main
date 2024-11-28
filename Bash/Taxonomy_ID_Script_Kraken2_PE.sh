for f in *R1_001.fastq.gz;
do

mkdir -p taxid

#Define Read1 and Read2
RD1=${f%R1_001.fastq.gz}R1_001.fastq.gz
RD2=${f%R1_001.fastq.gz}R2_001.fastq.gz

echo ""
echo ""
echo "Read1 is $RD1" 
echo "Read2 is $RD2"
echo ""
echo ""

#Taxonomic identification of reads
kraken2 --db viral_db --threads 20 --output ${f%_S*}.tsv --report --minimum-base-quality 20 --paired --use-names $RD1 $RD2 

mv *.tsv taxid

done