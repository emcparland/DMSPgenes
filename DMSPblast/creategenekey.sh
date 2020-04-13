#!/bin/sh
# Create a key for genes
awk -F\| '{print $2}' mmetspids.txt | sort |uniq > mmetsp_uniqsp.txt
awk '{print $0}' mmetsp_uniqsp.txt | while read headername
do
	cut -f2-18 DSYB_uniqsp.txt | grep $headername   > tmp
	length=`(wc -l < tmp | awk '{print $1}')`
	if [ $length -gt 0 ]
	then
		echo 1 >> DSYB
	else
		echo 0 >> DSYB
	fi
	
	cut -f2-18 tpmt1_uniqsp.txt | grep $headername > tmp
	length=`(wc -l < tmp | awk '{print $1}')`
	if [ $length -gt 0 ]
	then
		echo 1 >> TpMT1
	else
		echo 0 >> TpMT1
	fi
	
	cut -f2-18 tpmt2_uniqsp.txt | grep $headername > tmp
	length=`(wc -l < tmp | awk '{print $1}')`
	if [ $length -gt 0 ]
	then
		echo 1 >> TpMT2
	else
		echo 0 >> TpMT2
	fi
done
paste DSYB TpMT2 TpMT1 > key
rm tmp DSYB TpMT2 TpMT1
paste key mmetsp_uniqsp.txt > tmp
echo "DSYB	TpMT2	TpMT1	MMETSPstrain" | tr " " "\t" > header
cat header tmp > mmetspkey.txt
rm tmp header

# adding a space holder for the prokaryotes with dsyb
awk -F"\t" '{print $18}' dsybprok_uniqsp.txt | awk -F" " '{print $2,$3}' > prok_sp.txt
awk '{print $0}' prok_sp.txt | while read headername
do
	echo 1 >> DSYB
	echo 0 >> TpMT2
	echo 0 >> TpMT1
done
paste DSYB TpMT2 TpMT1 prok_sp.txt > key
cat mmetspkey.txt key > genekey_final.txt
rm DSYB TpMT2 TpMT1 key mmetspkey.txt key