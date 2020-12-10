#!/bin/sh
# Create a key for genes
awk -F\| '{print $2}' mmetspids.txt | sort |uniq > mmetsp_uniqsp.txt
printf Fragilariopsis-cylindrus-CCMP1102\\nChrysochromulina_tobin_CCMP291\\nPhaeodactylum-tricornutum-CCAP1055\\nThalassiosira-pseudonana-CCMP1335\\n >> mmetsp_uniqsp.txt

awk '{print $0}' mmetsp_uniqsp.txt | while read headername
do
	grep $headername final_DSYB.txt > tmp
	length=`(wc -l < tmp | awk '{print $1}')`
	if [ $length -gt 0 ]
	then
		echo 1 >> DSYB
	else
		echo 0 >> DSYB
	fi
	
	grep $headername final_TpMT2.txt > tmp
	length=`(wc -l < tmp | awk '{print $1}')`
	if [ $length -gt 0 ]
	then
		echo 1 >> TpMT2
	else
		echo 0 >> TpMT2
	fi
done
paste DSYB TpMT2 > key
rm tmp DSYB TpMT2
paste key mmetsp_uniqsp.txt > tmp
echo "DSYB TpMT2 MMETSPstrain" | tr " " "\t" > header
cat header tmp > mmetspkey.txt
rm tmp header

# adding a space holder for the prokaryotes with dsyb
awk -F"\t" '{print $18}' ../DMSPblast/dsybprok_uniqsp.txt | awk -F" " '{print $2,$3}' > prok_sp.txt
awk '{print $0}' prok_sp.txt | while read headername
do
	echo 1 >> DSYB
	echo 0 >> TpMT2
done
paste DSYB TpMT2 prok_sp.txt > key
cat mmetspkey.txt key > genekey_final.txt
rm DSYB TpMT2 key mmetspkey.txt
