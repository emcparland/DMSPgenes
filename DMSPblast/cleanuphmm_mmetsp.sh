#!/bin/sh
FILE=${1?Error: no name given}
NAME=${2?Error: no name given}

# filter to keep anything with evalue <1e-30
awk '$5 <=1e-30' $FILE > tmp_filter

# transcribe blastnames
cut -f1 -d" " tmp_filter |awk -F\| '{print $1}' | while read headername
do
	grep $headername mmetspids.txt >> temp_names.txt || printf $headername'\n' >> temp_names.txt
done
sed -i -e 's/>//g' temp_names.txt

paste tmp_filter temp_names.txt > tmp
echo "target - query - E-value score bias E-value score bias exp reg clu ov env dom rep inc targetdescrip mmetspid" > header
cat header tmp > "$NAME"_hmm.txt
rm temp_names.txt* tmp tmp_filter header

# keep only the best hit from each transcriptomes
# find the unique transcriptomes in blast results, calculate the length of results and length of unique names, if hits is longer than unique names proceed
awk '{print $27}' "$NAME"_hmm.txt |sort |uniq > unique_names.txt
a=`(wc -l < "$NAME"_hmm.txt | awk '{print $1}')`
a="$(($a-1))"
b=`(wc -l < unique_names.txt | awk '{print $1}')`

if [ $a -gt $b ]
then
# loop to find the results for each transcriptome
	awk '{print $1}' unique_names.txt | while read headername
	do
	grep "$headername" "$NAME"_hmm.txt | head -n1 >> "$NAME"_hmmuniq.txt
	done
fi
rm unique_names.txt


# Finally just keep the best result for each strain (some transcriptomes for the same strain in different conditions
awk '{print $27}' "$NAME"_hmmuniq.txt | awk -F\| '{print $2}' |sort |uniq |grep "^[A-Z]" > uniq.txt

awk '{print $1}' uniq.txt |while read headername
do
# find all of the best hits for this strain
grep "$headername" "$NAME"_hmmuniq.txt > tmp
# sort and select the result with the smallest evalue
search=`(awk '{print $5}' tmp | sort -g | head -n1)`
# find your numeric value of interest
grep "$search" tmp > final
# some of these are still identical, if so...
lengthf=`(wc -l < final |awk '{print $1}')`
if [ $lengthf -gt 1 ]
then
# keep the first result
	head -n1 final >> "$NAME"_hmmuniqsp.txt
else
# if not then just keep the single result
	cat final >> "$NAME"_hmmuniqsp.txt
fi
done
rm uniq.txt tmp final

