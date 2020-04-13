# Survey of DMSP synthesis genes

## Download local blast if you don't have it already with conda

```conda create --name localblasting```
```conda install -c bioconda blast```
```conda activate localblasting```

# tblastn DMSP synthesis genes to survey eukaryotes in MMETSP transcriptomes
The MMETSP is the most diverse database of marine protists. Here, I'm using the most up to date MMETSP transcriptomes, which were reanalyzed by Johnson et al. and includes transcriptomes (n=678) from primarily cultured unique strains (n=395) (some transcriptomes are from the same organism in different culturing conditions). 

## First create custom blast db with MMETSP transcriptomes
Download all of the transcriptomes (I used the [zenodo-get](https://zenodo.org/record/1261813#.XhkWCC-ZPq0) python script which will recursively download all of the transcriptomes)
```zenodo_get.py -r 1212585```

Check that you have 678 transcriptomes
```ls MMETSP* | wc -l```

Combine into a single file for blast
```cat *.fasta > mmetsp.fa```

Make your custom nt MMETSP db
```makeblastdb -dbtype nucl -in mmetsp.fa -out mmetsp_nt.db -title mmetsp_nt```

## tblastn sequences of interest
I used an indivdual sequence for TpMT2 (NCBI accession #) and TpMT1 (NCBI accession #) from [(Kageyama et al. 2018)](https://www.sciencedirect.com/science/article/abs/pii/S0003986118300080?via%3Dihub) and all of the DSYB sequences reported by [Curson et al. 2018](https://www.nature.com/articles/s41564-018-0119-5#Sec27) in the supplement

I was using the poseidon cluster at WHOI so request time with srun first

```srun -p scavenger --time=00:30:00 --ntasks-per-node=6 --mem=8gb --pty bash```

```tblastn -num_threads 4 -query TpMT2.fa -db ~/mmetsp/mmetsp_nt.db -out tpmt2_blastraw.txt -outfmt "6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore score sseq"```

```tblastn -num_threads 4 -query TpMT1.fa -db ~/mmetsp/mmetsp_nt.db -out tpmt1_blastraw.txt -outfmt "6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore score sseq"```

```tblastn -num_threads 6 -query allDSYB.fa -db ~/mmetsp/mmetsp_nt.db -out DSYB_blastraw.txt -outfmt "6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore score sseq"```

## Clean up tblastn results
*I kept results if they passed a criterion of e-value <= 1e-30 and query coverage >=70%*

Run shell script, input raw tblastn results and tell it name of gene you want to label files with:
```./cleanupblast_mmetsp.sh tpmt2_blastraw.txt tpmt2```

Note that this script: 
- Calculates the percent query coverage (divides column 6 by column 2, or alignment length/query sequence length). 
- Filters based on criterion. 
- Transcribes the names in blast results (which are just mmetsp ids and contig numbers) to MMETSP id and species name. 
- Some MMETSP species had multiple hits for the same protein (particularly DSYB because I searched with more than one query sequence). If the length of the output is > the length of the unique subject names then this script also filters to produce just one hit for each MMETSP id (the 'best' hit for each MMETSP id is chosen based on highest e-value, then highest query coverage). This produces a txt file named with your gene of interest (e.g. DSYB_uniq.txt). 
- Finally if a species has multiple MMETSP id's (meaning it was grown in different conditions and the resulting transcriptomes are represented in database) then filters to keep the top hit for each species (again chosen based on highest e-value, then highest query coverage). This produces an additional txt file (e.g. DSYB_uniqsp.txt)*.
- *This is information that could be of importance so make sure to edit a script like this to reflect your scientific questions.

# For the prokaryotes

```srun -p scavenger --time=00:30:00 --ntasks-per-node=6 --mem=8gb --pty bash```

```tblastn -query dsyb.fa -remote -db ref_prok_rep_genomes -max_target_seqs 500 -out dsyb_refprokgen_blastraw.txt -outfmt "6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore score sseq"```

```tblastn -query TpMT2.fa -remote -db ref_prok_rep_genomes -max_target_seqs 500 -out tpmt2_refprokgen_blastraw.txt -outfmt "6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore score sseq"```

```tblastn -query TpMT1.fa -remote -db ref_prok_rep_genomes -max_target_seqs 500 -out tpmt1_refprokgen_blastraw.txt -outfmt "6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore score sseq"```

## Need to get readable names for hits (output gives NZ id's)
Get the names of your significant hits
```awk '$13 <=1e-30 && $16 >=0.7' dsyb_refprokgen_blastraw.txt | cut -f3 | awk -F\| '{print $4}' > prok_search.txt
awk '$13 <=1e-30 && $16 >=0.7' tpmt2_refprokgen_blastraw.txt | cut -f3 | awk -F\| '{print $4}' >> prok_search.txt
awk '$13 <=1e-30 && $16 >=0.7' tpmt1_refprokgen_blastraw.txt | cut -f3 | awk -F\| '{print $4}' >> prok_search.txt
sort prok_search.txt | uniq > Prok_NZ_search.txt
rm prok_search.txt```

Search in nucleotide database on NCBI with [batch entrez](https://www.ncbi.nlm.nih.gov/sites/batchentrez)
```grep "^>" sequencetmp.fasta | awk -F" " '{print $1,$2,$3}' > refprokids.txt
rm Prok_NZ_search.txt```

## Clean up results
*I again kept results if they passed a criterion of e-value <= 1e-30 and query coverage >=70%*

Run shell script, input raw tblastn results and tell it name of gene you want to label files with:
```./cleanupblast_refprok.sh dsyb_refprokgen_blastraw.txt dsyb```

Note that this script: 
- Calculates the percent query coverage (divides column 6 by column 2, or alignment length/query sequence length). 
- Filters based on criterion. 
- Transcribes the NZ ids in blast results to NZ id + species name.
*TpMT2 had no significant prokaryote hits*

# Add manual searches
Here I also added a few manual searches that weren't represented in MMETSP (see Methods in paper for details).
```cat dsybprok_uniqsp.txt extra_proks_sig.txt
cat DSYB_uniqsp.txt extra_DSYB_sig.txt
cat tpmt2_uniqsp.txt extra_TpMT2_sig.txt```

# Cat proks and euks together and make fasta files
```cat DSYB_uniqsp.txt dsybprok_uniqsp.txt > DSYB_final.txt```
(There are no prok sequences for TpMT2)
```cp tpmt2_uniqsp.txt TpMT2_final.txt```

This command will grab your headers and sequences, remove any "-" in the MMETSP sequences, and print a fasta file with >
```cat TpMT2_final.txt  | awk -F "\t" '{print $18,$17}' | awk '{gsub(/-/,"",$2)}1'|awk 'BEGIN {OFS="\n"}{print ">"$1,$2}' > TpMT2_final.fa```

# Create a key for all of the mmetsp species
Binary columns where 1=gene present in strain

```./creategenekey.sh```

# Make Figure 1 and pheno-genotype correlation
I used this output to make Figure 1 and calculate the correlation of the phenotype and genotype in R
```source(Fig1code.R)```

# Recursive blast search
As we know the MMETSP transcriptomes were not axenic, we wanted to check if any of the sequences could be prokaryote contamination. To do this, I recursively searched each of the significant TpMT2 and DSYB eukaryotic sequences against the nr NCBI database, kept the top 100 hits for each sequence and determined if they were primarily prokaryote or eukaryote.

```srun -p scavenger --time=02:00:00 --ntasks-per-node=6 --mem=8gb --pty bash```

```tblastn -db nr -query TpMT2_final.fa -remote -out TpMT2_reblast.txt -max_hsps 1 -max_target_seqs 100 -outfmt "6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore score"```

```tblastn -db nr -query DSYB_uniqsp.fa -remote -out DSYB_reblast.txt -max_hsps 1 -max_target_seqs 100 -outfmt "6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore score"```

Get the id's and search on [batch entrez](https://www.ncbi.nlm.nih.gov/sites/batchentrez)
```cat *_reblast_all.txt | cut -f3 | awk -F \| '{print $4}' |sort | uniq >tmp.txt```
Download summary file and pull out names
```awk 'NR%4==1' tmp_summary.txt | awk -F" " '{print $2,$3}' >tmp_taxon.txt```
Then search for [taxon id](https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi)
```awk -F\| '{print $4}' tax_report.txt >taxid_search.txt```
Retrieve taxonomy from batch entrez again and download summary to get taxonomy result, then parse to make easily searchable
```awk 'NR%2==1' taxonomy_result.txt | awk -F" " '{print $2,$3}' > tmp1```
```awk 'NR%2==0' taxonomy_result.txt |awk -F"," '{print $2}' > tmp2```
```paste tmp1 tmp2 > search1.txt```

Now run a for loop to assign kingdom to blast hit:
Get list of accession numbers and species names
```awk 'NR%4==3' tmp_summary.txt | awk -F" " '{print $1}' >tmp_accesion.txt```
```paste tmp_accesion.txt tmp_taxon.txt >search2.txt```
Clean up
```rm tmp*```





