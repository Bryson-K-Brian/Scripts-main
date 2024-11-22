mkdir -p metagenomic_assembly

cd trimmed_reads

for f in *R1_001_trimmed.fastq.gz;
do

#Define Read1 and Read2
RD1=${f%R1_001_trimmed.fastq.gz}R1_001_trimmed.fastq.gz
RD2=${f%R1_001_trimmed.fastq.gz}R2_001_trimmed.fastq.gz

echo ""
echo ""
echo "Read1 is $RD1" 
echo "Read2 is $RD2"
echo ""
echo ""

#Meatgenomic Assembly with Megahit
megahit -1 $RD1 -2 $RD2 -o ../metagenomic_assembly/${f%_S*}

echo ""
echo ""
echo "Sample ${f%_S*} Complete"
echo ""
echo ""

done

cd ..

