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
cut -f3 tmp_filter |awk -F- '{print $1}' | while read headername
do
	grep $headername mmetspids.txt >> temp_names.txt || printf $headername'\n' >> temp_names.txt
done
sed -i -e 's/>//g' temp_names.txt

paste tmp_filter temp_names.txt > tmp
echo "qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore score pqlen sseq mmetspid" | tr " " "\t" > header
cat header tmp > "$NAME"_hits.txt
rm temp_names.txt* tmp tmp_filter header

# keep only the best hit from each transcriptomes
# find the unique transcriptomes in blast results, calculate the length of results and length of unique names, if hits is longer than unique names proceed
cut -f18 "$NAME"_hits.txt |sort |uniq > unique_names.txt
a=`(wc -l < "$NAME"_hits.txt | awk '{print $1}')`
a="$(($a-1))"
b=`(wc -l < unique_names.txt | awk '{print $1}')`

if [ $a -gt $b ]
then
# loop to find the results for each transcriptome
	awk '{print $1}' unique_names.txt | while read headername
	do
	grep $headername "$NAME"_hits.txt > tmp
# first determine if any evalues ==0, if so collect these results
	awk '{if($13 == 0.0)print $0}' tmp > tmp0
# if there were evalues==0 then...
	length0=`(wc -l < tmp0 |awk '{print $1}')`
	if [ $length0 -gt 0 ]
	then
# sort based on the percent query coverage and select the result with the highest value
		search=`(cut -f16 tmp0 | sort -n | tail -n1)`
	else
# if no 0 evalues, select the result with the smallest evalue
		search=`(cut -f13 tmp | sort -g | head -n1)`
	fi
# find your numeric value of interest
	grep "$search" tmp > final
# some of these are still identical, if so...
	lengthf=`(wc -l < final |awk '{print $1}')`
	if [ $lengthf -gt 0 ]
	then
# keep the first result
		head -n1 final >> "$NAME"_uniq.txt
	else
# if not then just keep the single result
		cat final >> "$NAME"_uniq.txt
	fi
	done
fi
rm final tmp tmp0 unique_names.txt


# Finally just keep the best result for each strain (some transcriptomes for the same strain in different conditions
cut -f18 "$NAME"_uniq.txt | awk -F\| '{print $2}' |sort |uniq |grep "^[A-Z]" > uniq.txt

awk '{print $1}' uniq.txt |while read headername
do
# find all of the best hits for this strain
	grep $headername "$NAME"_uniq.txt > tmp
# first determine if any evalues ==0, if so collect these results
	awk '{if($13 == 0.0)print $0}' tmp > tmp0
# if there were evalues==0 then...
	length0=`(wc -l < tmp0 |awk '{print $1}')`
	if [ $length0 -gt 0 ]
	then
# sort based on the percent query coverage and select the result with the highest value
		search=`(cut -f16 tmp0 | sort -n | tail -n1)`
	else
# if no 0 evalues, select the result with the smallest evalue
		search=`(cut -f13 tmp | sort -g | head -n1)`
	fi
# find your numeric value of interest
	grep "$search" tmp > final
# some of these are still identical, if so...
	lengthf=`(wc -l < final |awk '{print $1}')`
	if [ $lengthf -gt 1 ]
	then
# keep the first result
		head -n1 final >> "$NAME"_uniqsp.txt
	else
# if not then just keep the single result
		cat final >> "$NAME"_uniqsp.txt
	fi
done
rm uniq.txt tmp tmp0 final

