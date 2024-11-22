mkdir -p readmapping

for f in *_R1_001.fastq.gz;
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

#Assembly with Spades
bwa mem -t 6 P_mirabilis_chromosome_reference.fna $RD1 $RD1 -o ${f%_T}.sam
samtools view -b ${f%_T}.sam > ${f%_T}.bam
samtools view -bF 12 ${f%_T}.bam >${f%_T}_12F.bam
samtools sort ${f%_T}_12F.bam > ${f%_T}_12F_sorted.bam
samtools depth ${f%_T}_12F_sorted.bam > ${f%_T}_12F_sorted.bam.txt
python Draw_SequencingDepth.py ${f%_T}_12F_sorted.bam.txt Depth


echo ""
echo ""
echo "Sample ${f%_T} Complete"
echo ""
echo ""

mkdir -p readmapping/${f%_T}
mv ${f%_T}* readmapping/${f%_T}

done
