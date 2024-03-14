# run EDTA

These instructions are based on https://github.com/susbo/Drosophila_unistrand_clusters/tree/main/Transposon_libraries

The examples below assumes that you are in `species/Dmel/dm6/annotation/EDTA` (replace Dmel/dm6 for your species/assembly) and that `species/Dmel/dm6/genome.fa` is available.

Since we are only going to use LTR/Gypsy annotations, you may run EDTA end-to-end

```
EDTA.pl --genome ../../genome.fa --sensitive 1 --anno 1 --evaluate 1 -t 10 --overwrite 1 --force 1
```

or only run the LTR module by using:
```
EDTA_raw.pl --genome ../../genome.fa --type ltr -t 10 --overwrite 0
EDTA.pl --genome ../../genome.fa --sensitive 1 --anno 1 --evaluate 1 -t 10 --overwrite 0 --force 1
```

Our version of EDTA crashed whenever a Penelope elements was detected. The avoid this, one of the scripts was slightly modified. This is described in `EDTA_fix`.
