#!/bin/bash

### For DB access:
# shopt -s expand_aliases
# source ~/.bash_aliases

ofile=statistic_TTA-vs-FS_genes.tsv

[ -s $ofile ] && mv $ofile ${ofile}~

header='--header'

for gbff_file in `ls -1 ALL_genomes/*.gbff`
do
	echo Run: $gbff_file

	./search_TTA_genes.pl $gbff_file --auto --threads=15 -output=stdout $header >> $ofile
	header=''
done
