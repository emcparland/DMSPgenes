#!/bin/sh
NAME=${1?Error: no name given}
touch tmp.txt

cut -f4 $NAME | grep -v "MMETSPstrain" |awk '{print $1}' |while read headername

do
line1="grep '$headername' mmetsp_all_lineages.csv |head -n1 |cut -d "," -f6-13 >> tmp.txt"

eval $line1

done

echo "genus,family,order,class,norank,norank,kingdom,generic" | cat - tmp.txt > tmp_lineages.txt
rm tmp.txt