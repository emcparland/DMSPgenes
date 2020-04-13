
Install: [muscle](https://www.drive5.com/muscle/) with conda and [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/)

# Step 1: Align sequences with muscle
```muscle -in DSYB_final.fa -out DSYB_align.fa -seqtype protein```
```muscle -in TpMT2_final.fa -out TpMT2_align.fa -seqtype protein```

# Step 2: Trim with Geneious or Jalview and refine alignment
```muscle -in DSYB_align_trim.fa -out DSYB_align_final.fa -seqtype protein -refine```
```muscle -in TpMT2_align_trim.fa -out TpMT2_align_final.fa -seqtype protein -refine```

Convert to phylip file for raxml
```perl Fasta2Phylip.pl DSYB_align_final.fa DSYB_align_final.phy```
```perl Fasta2Phylip.pl TpMT2_align_final.fa TpMT2_align_final.phy```

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
