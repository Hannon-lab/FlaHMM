#!/bin/bash
# Convert the .out file to tab-separated format

# Author: Susanne Bornelöv
# Last edited: 2023-02-21 (documentation and comments)

mkdir -p out/02_tab/{all,plus,minus}

genomes=`ls -d species/*/*`

for genome in $genomes
do
   echo GENOME: $genome

   species=`echo $genome | cut -d'/' -f2`
   build=`echo $genome | cut -d'/' -f3`

   hgLoadOut -tabFile=out/02_tab/all/$species.$build.tab -nosplit test out/01_fix_header/$species.$build.out

	cat out/02_tab/all/$species.$build.tab | awk -v OFS="\t" '$10=="+"' > out/02_tab/plus/$species.$build.tab
	cat out/02_tab/all/$species.$build.tab | awk -v OFS="\t" '$10=="-"' > out/02_tab/minus/$species.$build.tab
done
