mkdir abricate

for f in *.fasta;
do

abricate --csv --db card $f > abricate/${f%.fasta}_card.csv

echo "sample ${f%.fasta} card done"

done

for f in *.fasta;
do

abricate --csv --db plasmidfinder $f > abricate/${f%.fasta}_plasmidfinder.csv

echo "sample ${f%.fasta} plasmidfinder done"

done

for f in *.fasta;
do

abricate --csv --db vfdb $f > abricate/${f%.fasta}_vfdb.csv

echo "sample ${f%.fasta} vfdb done"

done

for f in *.fasta;
do

abricate --csv --db resfinder $f > abricate/${f%.fasta}_resfinder.csv

echo "sample ${f%.fasta} vfdb done"

done