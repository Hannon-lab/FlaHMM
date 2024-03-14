#!/bin/bash
# Convert tab-separated files to bed, but keep only LTR transposons

# Author: Susanne Bornelöv
# Last edited: 2023-02-21 (documentation and comments)

mkdir -p out/03_bed/{plus,minus}
mkdir -p out/03_bed/{plus,minus}/merged

genomes=`ls -d species/*/*`

for genome in $genomes
do
   echo GENOME: $genome

   species=`echo $genome | cut -d'/' -f2`
   build=`echo $genome | cut -d'/' -f3`

	# Remove "Unknown" repeats that may be repeat genes (e.g., Hist*)
	scripts/toBed6+10.pl out/02_tab/plus/$species.$build.tab | LC_COLLATE=C sort -k1,1 -k2,2n | awk -v OFS="\t" '$13=="Gypsy"' > out/03_bed/plus/$species.$build.bed
	scripts/toBed6+10.pl out/02_tab/minus/$species.$build.tab | LC_COLLATE=C sort -k1,1 -k2,2n | awk -v OFS="\t" '$13=="Gypsy"' > out/03_bed/minus/$species.$build.bed

#	scripts/toBed6+10.pl out/02_tab/plus/$species.$build.tab | LC_COLLATE=C sort -k1,1 -k2,2n | awk -v OFS="\t" '$12=="LTR"' > out/03_bed/plus/$species.$build.bed
#	scripts/toBed6+10.pl out/02_tab/minus/$species.$build.tab | LC_COLLATE=C sort -k1,1 -k2,2n | awk -v OFS="\t" '$12=="LTR"' > out/03_bed/minus/$species.$build.bed

	bedtools merge -i out/03_bed/plus/$species.$build.bed > out/03_bed/plus/merged/$species.$build.bed
	bedtools merge -i out/03_bed/minus/$species.$build.bed > out/03_bed/minus/merged/$species.$build.bed
done
