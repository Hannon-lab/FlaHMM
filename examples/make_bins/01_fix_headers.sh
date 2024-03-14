#!/bin/bash
# EDTA will cut chromosomes longer than a fixed length. Here we restore the chromosome names.

# Author: Susanne Bornelöv
# Last edited: 2023-02-21 (documentation and comments)

mkdir -p out
mkdir -p out/01_fix_header

genomes=`ls -d species/*/*`

for genome in $genomes
do
   echo GENOME: $genome

   species=`echo $genome | cut -d'/' -f2`
   build=`echo $genome | cut -d'/' -f3`

   if [ ! -e "$genome/annotation/EDTA/genome.fa.mod.EDTA.anno/genome.fa.mod.out" ]; then
      echo Input file missing for $genome
      continue;
   fi

   mkdir -p out/01_fix_header
   cat $genome/annotation/EDTA/genome.fa.mod.EDTA.anno/genome.fa.mod.out > out/01_fix_header/$species.$build.out

	if [ -e "$genome/annotation/EDTA/genome.fa" ]; then
	   # Ugly hack, but replaces headers that EDTA shortened with original contig name
   	cat $genome/annotation/EDTA/genome.fa | grep "^>" | sed -e 's/>//' | cut -d' ' -f1 > headers1.txt
	   cat $genome/annotation/EDTA/genome.fa.mod | grep "^>" | sed -e 's/>//' | cut -d' ' -f1 > headers2.txt
   	paste headers1.txt headers2.txt > headers.txt
	   perl scripts/fixHeaders.pl out/01_fix_header/$species.$build.out headers.txt
	else
		echo Skipping header correction step due to $genome/annotation/EDTA/genome.fa missing...
		echo This is usually not a problem unless the genome assembly contains unusually long contig names.
	fi
done

rm -f headers.txt headers{1,2}.txt
