
Install: [muscle](https://www.drive5.com/muscle/) with conda and [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/)

# Step 1: Align all sequences with muscle

```cat DSYBall_hmmuniqsp.fa TpMT2_hmmuniqsp.fa > all_final.fa```

Transcribe ORF names to MMETSP species
```awk '{print $1}' mmetspids.txt | while read headername```
```do```
```search=`(echo $headername | awk -F "\|" '{print $1}')` ```
```sed -i -e "s/$search.*/$headername/g" all_final.fa```
```done```

Make alignment with all sequences:
```muscle -in all_final.fa -out all_align.fa -seqtype protein```

Make alignment with only sequences with a phenotype and genotype:

```awk '{print $1}' names_phenogeno.txt | while read headername```
```do```
```search=`(echo $headername | awk -F "|" '{print $2}')` ```
```grep -A1 $search all_final.fa >> all_pheno.fa```
```done```

```muscle -in all_pheno.fa -out all_pheno_align.fa -seqtype protein```

# Step 2: Trim with Geneious or Jalview and refine alignment
```muscle -in all_align.fa -out all_align_final.fa -seqtype protein -refine```
```muscle -in all_pheno_align.fa -out all_pheno_align_final.fa -seqtype protein -refine```

Convert to phylip file for raxml
```perl Fasta2Phylip.pl all_align_final.fa all_align_final.phy```
```perl Fasta2Phylip.pl all_pheno_align_final.fa all_pheno_align_final.phy```

# Step 3: Build tree with RAXML

## First construct the ML tree, using the AUTO function which automiatcally determines the best protein model with respect to the likelihood on a fixed, reasonable tree
Note: I tried to include outgroup here but I don't think it liked the name and still trying to figure out why. Either way, outgroup effects drawing not actual tree

```raxmlHPC -p 12345 -s $phy_file -n A1 -m PROTGAMMAAUTO -# 1```

## Now repeat with the selected model, but also find the best-scoring tree by generating 20 trees from distinct starting trees
```raxmlHPC -p 12345 -s $phy_file -n A2 -m PROTGAMMALG -# 20```

## Next get support values, I used the frequenced based stopping criterion instead of setting chosen bootstraps
```raxmlHPC -p 12345 -s $phy_file -n A3 -m PROTGAMMALG -b 12345 -# autoFC```

## Finally draw bootstrap values onto the maximum likelihood tree to create bipartitions (support values assigned to nodes)
```raxmlHPC -p 12345 -n A4 -m PROTGAMMALG -f b -t RAxML_bestTree.A2 -z RAxML_bootstrap.A3```

## Can also build a consensus tree with bootstraps (strict (STRICT), majority (MR) or extended majority (MRE) rules)
```raxmlHPC -p 12345 -n A5 -m PROTGAMMALG -J MR -z RAxML_bootstrap.A3```
