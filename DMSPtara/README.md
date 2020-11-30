
Install: hmmer

## Make a hmm profile for DSYB and TpMT2
I used the alignment files created for the phylogenetic trees to build the hmm profiles:

```hmmbuild DSYB.hmm DSYB_align_final.fa```

```hmmbuild TpMT2.hmm TpMT2_align_final.fa```

I then used the [Ocean Gene Atlas](http://tara-oceans.mio.osupytheas.fr/ocean-gene-atlas/) to search for an hmm profile search of the MATOU (eukaryotic Unigenes built from metatranscriptomes), and the two R scripts provided to filter the results and make maps (maps were made modifying Meren lab script found here: https://github.com/merenlab/world-map-r).

The Ocean Gene Atlas: exploring the biogeography of plankton genes online. E. Villar, T. Vannier, C. Vernette, M. Lescot, M. Cuenca, A. Alexandre, P. Bachelerie, T. Rosnet, E. Pelletier, S. Sunagawa, P. Hingamp. (2018) Nucleic Acids Research. 
doi: 10.1093/nar/gky376
