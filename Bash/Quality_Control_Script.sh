mkdir -p trimmed_reads


for f in raw_reads/*01.fastq.gz;
do

#Define Read1 and Read2
RD1=${f%01.fastq.gz}01.fastq.gz
RD2=${f%01.fastq.gz}02.fastq.gz

echo ""
echo ""
echo "Read1 is $RD1" 
echo "Read2 is $RD2"
echo ""
echo ""

trm1=${f%01.fastq.gz}01_trimmed.fastq.gz
trm2=${f%01.fastq.gz}02_trimmed.fastq.gz

#Quality Control using BBDUK
bbduk.sh in=$RD1 in2=$RD2 out=$trm1 out2=$trm2 trimq=20 qtrim=rl minlength=50 minbasequality=0

echo ""
echo ""
echo "Sample ${f%_S*} Complete"
echo ""
echo ""

done

mv raw_reads/*trimmed* trimmed_reads
