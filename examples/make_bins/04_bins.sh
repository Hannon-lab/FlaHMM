#!/bin/bash
# Prepare 5kb genomic bins and calculate (LTR) coverage within each bin

# Author: Susanne Bornelöv
# Last edited: 2023-02-21 (documentation and comments)
# Last edited: 2024-07-10 (removed if clause to ensure bins are always recalculated)

mkdir -p out/04_bins/{bins,overlap}

genomes=`ls -d species/*/*`

for genome in $genomes
do
   echo GENOME: $genome

   species=`echo $genome | cut -d'/' -f2`
   build=`echo $genome | cut -d'/' -f3`

	# Create windows
	for size in 2500 5000 10000
	do
		if [[ "$size" == 10000 ]]; then
		   bedtools makewindows -w $size -s 5000 -g $genome/genome.fa.fai | awk -v OFS="\t" '$3-$2>0' | gzip -c > out/04_bins/bins/$species.$build.$size.bed.gz
		else
		   bedtools makewindows -w $size -g $genome/genome.fa.fai | awk -v OFS="\t" '$3-$2>0' | gzip -c > out/04_bins/bins/$species.$build.$size.bed.gz
		fi

		grp=$(echo "scale=1; $size/1000" | bc)
		grp=`echo $grp | sed -e 's/\\.0//'`

		mkdir -p out/04_bins/overlap/bins_${grp}k/{plus,minus}
		bedtools coverage -a <(zcat out/04_bins/bins/$species.$build.$size.bed.gz) -b out/03_bed/plus/merged/$species.$build.bed > out/04_bins/overlap/bins_${grp}k/plus/$species.$build.bed
		bedtools coverage -a <(zcat out/04_bins/bins/$species.$build.$size.bed.gz) -b out/03_bed/minus/merged/$species.$build.bed > out/04_bins/overlap/bins_${grp}k/minus/$species.$build.bed
	done

done
