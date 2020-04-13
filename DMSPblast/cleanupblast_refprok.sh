#!/bin/sh
FILE=${1?Error: no name given}
NAME=${2?Error: no name given}

# add column for percent query length
awk '{print $6/$2}' $FILE > q
cut -f16 $FILE > seq
cut -f1-15 $FILE > alnresults
paste alnresults q seq > tmp
rm alnresults q seq

# filter to keep anything with evalue <1e-30 and query length >=70%
awk '$13 <=1e-30 && $16 >=0.7' tmp > tmp_filter
rm tmp

# transcribe blastnames
cut -f3 tmp_filter |awk -F\| '{print $4}' | while read headername
do
	grep $headername refprokids.txt >> temp_names.txt || printf $headername'\n' >> temp_names.txt
done
sed -i -e 's/>//g' temp_names.txt

paste tmp_filter temp_names.txt > tmp
echo "qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore score pqlen sseq mmetspid" | tr " " "\t" > header
cat header tmp > "$NAME"prok_uniqsp.txt

rm temp_names* tmp tmp_filter header