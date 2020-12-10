# Survey of DMSP synthesis genes

# prokaryote searches in NCBI ref seq genomes

Create a local conda environment with blast to access the NCBI db remotely.

## tblastn ref seq genomes at ncbi

```srun -p scavenger --time=00:30:00 --ntasks-per-node=6 --mem=8gb --pty bash```

```tblastn -query dsyb.fa -remote -db ref_prok_rep_genomes -max_target_seqs 500 -out dsyb_refprokgen_blastraw.txt -outfmt "6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore score sseq"```

```tblastn -query TpMT2.fa -remote -db ref_prok_rep_genomes -max_target_seqs 500 -out tpmt2_refprokgen_blastraw.txt -outfmt "6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore score sseq"```

```tblastn -query TpMT1.fa -remote -db ref_prok_rep_genomes -max_target_seqs 500 -out tpmt1_refprokgen_blastraw.txt -outfmt "6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore score sseq"```

### Need to get readable names for hits (output gives NZ id's)
Get the names of your significant hits
```awk '$13 <=1e-30 && $16 >=0.7' dsyb_refprokgen_blastraw.txt | cut -f3 | awk -F\| '{print $4}' > prok_search.txt```

```awk '$13 <=1e-30 && $16 >=0.7' tpmt2_refprokgen_blastraw.txt | cut -f3 | awk -F\| '{print $4}' >> prok_search.txt```

```awk '$13 <=1e-30 && $16 >=0.7' tpmt1_refprokgen_blastraw.txt | cut -f3 | awk -F\| '{print $4}' >> prok_search.txt```

```sort prok_search.txt | uniq > Prok_NZ_search.txt```

```rm prok_search.txt```

Search in nucleotide database on NCBI with [batch entrez](https://www.ncbi.nlm.nih.gov/sites/batchentrez)
```grep "^>" sequencetmp.fasta | awk -F" " '{print $1,$2,$3}' > refprokids.txt```

```rm Prok_NZ_search.txt```

### Clean up results
*For prokaryotes, I kept results if they passed a criterion of e-value <= 1e-30 and query coverage >=70%*

Run shell script, input raw tblastn results and tell it name of gene you want to label files with:
```./cleanupblast_refprok.sh dsyb_refprokgen_blastraw.txt dsyb```

this script: 
- Calculates the percent query coverage (divides column 6 by column 2, or alignment length/query sequence length). 
- Filters based on criterion. 
- Transcribes the NZ ids in blast results to NZ id + species name.

*TpMT2 had no significant prokaryote hits*


# Survey eukaryotes in MMETSP transcriptomes
The MMETSP is the most diverse database of marine protists. Here, I'm using the most up to date MMETSP transcriptomes, which were reanalyzed by Johnson et al. and includes transcriptomes (n=678) from primarily cultured unique strains (n=395) (some transcriptomes are from the same organism in different culturing conditions). 

### Download local blast if you don't have it already with conda

```conda create --name localblasting```
```conda install -c bioconda blast```
```conda activate localblasting```

### Create custom blast db with MMETSP transcriptomes
Download all of the transcriptomes (I used the [zenodo-get](https://zenodo.org/record/1261813#.XhkWCC-ZPq0) python script which will recursively download all of the transcriptomes)
```zenodo_get.py -r 1212585```

Check that you have 678 transcriptomes
```ls MMETSP* | wc -l```

Combine into a single file for blast
```cat *.fasta > mmetsp.fa```

Make your custom nt MMETSP db
```makeblastdb -dbtype nucl -in mmetsp.fa -out mmetsp_nt.db -title mmetsp_nt```

## Part 1: tblastn sequences of interest
I used an indivdual sequence for TpMT2 (NCBI accession #) and TpMT1 (NCBI accession #) from [(Kageyama et al. 2018)](https://www.sciencedirect.com/science/article/abs/pii/S0003986118300080?via%3Dihub) and all of the DSYB sequences reported by [Curson et al. 2018](https://www.nature.com/articles/s41564-018-0119-5#Sec27) in the supplement

```tblastn -num_threads 4 -query TpMT2.fa -db ~/mmetsp/mmetsp_nt.db -out tpmt2_blastraw.txt -outfmt "6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore score sseq"```

```tblastn -num_threads 4 -query TpMT1.fa -db ~/mmetsp/mmetsp_nt.db -out tpmt1_blastraw.txt -outfmt "6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore score sseq"```

```tblastn -num_threads 6 -query allDSYB.fa -db ~/mmetsp/mmetsp_nt.db -out DSYB_blastraw.txt -outfmt "6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore score sseq"```

### Clean up tblastn results
*For eukaryotes, I used thresholds of e-value <=2e-15 and query coverage >=75% to get the seed sequences for an hmm profile*

Run shell script, input raw tblastn results and tell it name of gene you want to label files with:
```./cleanupblast_mmetsp.sh tpmt2_blastraw.txt tpmt2```
```./cleanupblast_mmetsp.sh DSYB_blastraw.txt DSYB```

This script will: 
- Calculate the percent query coverage (divides column 6 by column 2, or alignment length/query sequence length). 
- Filter based on criterion. 
- Transcribe the names in blast results (which are just mmetsp ids and contig numbers) to MMETSP id and species name. 
- Some MMETSP species had multiple hits for the same protein (particularly DSYB because I searched with more than one query sequence). If the length of the output is > the length of the unique subject names then this script also filters to produce just one hit for each MMETSP id (the 'best' hit for each MMETSP id is chosen based on highest e-value, then highest query coverage). This produces a txt file named with your gene of interest (e.g. DSYB_uniq.txt). 
- Finally, if a species has multiple MMETSP id's (meaning it was grown in different conditions and the resulting transcriptomes are represented in database) then filters to keep the top hit for each species (again chosen based on highest e-value, then highest query coverage). This produces an additional txt file (e.g. DSYB_uniqsp.txt)*.

## Part 2: Create hmm profiles and perform hmm searches of MMETSP proteins
### Create fasta files
```cat TpMT2_final.txt  | awk -F "\t" '{print $18,$17}' | awk '{gsub(/-/,"",$2)}1'|awk 'BEGIN {OFS="\n"}{print ">"$1,$2}' > TpMT2_final.fa```

```cat DSYB_final.txt | awk -F "\t" '{print $18,$17}' | awk '{gsub(/-/,"",$2)}1'|awk 'BEGIN {OFS="\n"}{print ">"$1,$2}' > DSYB_final.fa```

### Create alignment with muscle and then manually clean up with jalview
```muscle -in TpMT2_final.fa -out TpMT2_align.fa -seqtype protein```

```muscle -in DSYB_final.fa -out DSYB_align.fa -seqtype protein```

### Build hmm profiles
```hmmbuild TpMT2.hmm TpMT2_align.fa```

```hmmbuild DSYB.hmm DSYB_align.fa```

### Create MMETSP protein database for hmmsearch
Download the MMETSP proteins
```wget https://zenodo.org/record/3247846/files/mmetsp_dib_trinity2.2.0_pep_zenodo.tar.gz```

Compile a fasta file of all the proteins. I add the MMETSP identifier at the beginning of each transcript so I can identify them after they are compiled:

```ls MMETSP* | while read headername```
```do```
	```name=$(awk '{print FILENAME}' $headername |head -n1 | awk -F'.' '{print $1}')```
	```sed "s/>/\>$name|/g" $name.trinity_out_*.pep >> mmetsp_pep.fa```
```done```

### Search proteins
```hmmsearch --tblout TpMT2.out TpMT2.hmm mmetsp_pep.fa```

```hmmsearch --tblout DSYB.out DSYB.hmm mmetsp_pep.fa```

### removing some of the extra table parts and add back headers
```head -n13421 TpMT2.out | tail -n13418 > tmp```
```echo "target - query - E-value score bias E-value score bias exp reg clu ov env dom rep inc targetdescrip"> header```
```cat header tmp > TpMT2_rawhmm.txt```
```rm header tmp```

### filter for e <=1e-30 and keep only the unique species (similar to cleanup script above)
```./cleanuphmm_mmetsp.sh TpMT2_rawhmm.txt TpMT2```

```./cleanuphmm_mmetsp.sh DSYB_rawhmm.txt DSYB```

### Include the manual searches for species not represented in MMETSP.
```cat extra_proks_sig.txt >> dsybprok_uniqsp.txt```

```cat extra_TpMT2_sig.txt >> TpMT2_hmmuniqsp.txt```

```cat extra_DSYB_sig.txt >> DSYB_hmmuniqsp.txt```

### Cat proks and euks together and make fasta files
```cat DSYB_hmmuniqsp.txt dsybprok_uniqsp.txt > DSYBall_hmmuniqsp.txt```

### create corresponding fasta file
First create a temporary searchable file (the original pep file has multiple lines of sequences)

```tr -d '\n' < mmetsp_pep.fa > tmp2```
```tr ">" "\n>" <tmp2 > tmp3```
```sed -e 's/)/)\n/g' tmp3 > pep_search```

```awk '{print $1}' TpMT2_hmmuniqsp.txt | while read headername```
```do```
```grep -A2 $headername pep_search | sed -e '2,2d' >> TpMT2_hmmuniqsp.fa```
Add a carot
```sed -i -e 's/MMETSP/>MMETSP/g' TpMT2_hmmuniqsp.fa```
Repeat for DSYB

Make sure to include the manual searches. This command will grab your headers and sequences, remove any "-" in the MMETSP sequences, and print a fasta file with >

```cat extra_TpMT2_sig.txt  | awk -F "\t" '{print $18,$17}' | awk '{gsub(/-/,"",$2)}1'|awk 'BEGIN {OFS="\n"}{print ">"$1,$2}' >> TpMT2_hmmuniqsp.fa```

# Recursive blast search
We wanted to check if any of the MMETSP sequences could be prokaryote contamination. To do this, I recursively searched each of the significant TpMT2 and DSYB eukaryotic sequences against the nr NCBI database, kept the top 100 hits for each sequence and determined if they were primarily prokaryote or eukaryote.

```gene=TpMT2``` (repeat with gene=DSYB)

```grep "^>" "$gene"_hmmuniqsp.fa > headers.txt```

```awk '{print $0}' headers.txt | while read headername```
```do```
```grep -A1 "$headername" "$gene"_hmmuniqsp.fa > tmp.fa```
```blastp -db nr -query tmp.fa -remote -out result.txt -max_hsps 1 -max_target_seqs 100 -outfmt "6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore score"```
```cat result.txt >> "$gene"_reblast.txt```

### Assign taxonomy 
I used the R package[Taxonomzir](https://cran.r-project.org/web/packages/taxonomizr/index.html). This requires some previous setup and long downloads.

```library(taxonomizr)```

```results <- read.delim("TpMT2_reblast.txt", header = F,stringAsFactors=F)```

```accessid <- as.vector(results$accession)```

```taxaId <- accessionToTaxa(accessid, 'accessionTaxa.sql')```

```taxa <- getTaxonomy(taxaId,'accessionTaxa.sql')```

```df <- cbind(accessid, taxaId,taxa)```
```write.csv(df,file = "TpMT2_reblasttaxa.csv")```

A few id's are missing. NCBI's [batch entrez](https://www.ncbi.nlm.nih.gov/sites/batchentrez) is helpful for checking those manually. Finally, remove any results deemed contamination by the analysis.

# Create a key for all of the mmetsp species
Binary columns where 1=gene present in strain, use to make Figure 1

```./creategenekey.sh```

# Make Figure 1, Supp Table 1 and perform and pheno-genotype correlation in R:
```source(Fig1code.R)```