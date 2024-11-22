mkdir -p trimmed_reads


for f in raw_reads/*R1_001.fastq.gz;
do

#Define Read1 and Read2
RD1=${f%R1_001.fastq.gz}R1_001.fastq.gz
RD2=${f%R1_001.fastq.gz}R2_001.fastq.gz

echo ""
echo ""
echo "Read1 is $RD1" 
echo "Read2 is $RD2"
echo ""
echo ""

trm1=${f%R1_001.fastq.gz}R1_001_trimmed.fastq.gz
trm2=${f%R1_001.fastq.gz}R2_001_trimmed.fastq.gz

#Quality Control using BBDUK
bbduk.sh in=$RD1 in2=$RD2 out=$trm1 out2=$trm2 trimq=20 qtrim=rl minlength=50 minbasequality=0

echo ""
echo ""
echo "Sample ${f%_S*} Complete"
echo ""
echo ""

done

mv raw_reads/*trimmed* trimmed_reads
