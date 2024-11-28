mkdir -p taxid

for f in *.fa;
do


echo ""
echo "contigs are ${f%.fa}" 
echo ""

#Taxonomic identification of reads
kraken2 --db viral_db --threads 20 --output ${f%.fa}.tsv --report --use-names $f 

mv ${f%.fa}.tsv taxid

done