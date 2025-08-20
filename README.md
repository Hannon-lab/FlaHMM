# FlaHMM 
FlaHMM: unistrand *flamenco*-like piRNA cluster prediction in *Drosophila* species

### Installation

FlaHMM can be cloned from our GitHub repository using the following command
```
git clone https://github.com/Hannon-lab/FlaHMM.git
```

To install all Python dependencies, we recommend creating a conda environment and install them using pip (or pip3, if applicable)

```
conda create -n FlaHMM python=3.8
conda activate FlaHMM

pip install pandas
pip install scikit-learn
pip install "hmmlearn==0.2.7"
pip install matplotlib
pip install tqdm
pip install seaborn
pip install openpyxl
```
FlaHMM has primarily been tested using hmmlearn v0.2.7. If needed, FlaHMM can also run with hmmlearn 0.2.8 or newer, however, please note that MultinomialHMM has been renamed to CategoricalHMM.

### Downloading transposon annotations

FlaHMM requires information about Gypsy transposon content across the target genome. The easiest option will be to download our pre-computed transposon content files using the following command.

```
cd FlaHMM
wget https://content.cruk.cam.ac.uk/ghlab/Susanne/FlaHMM/bins.tar.gz
tar xf bins.tar.gz
```

This includes pre-processed files for 193 genome assemblies from 119 species. All full list of all pre-processed genome assemblies and species is available [here](data/precomputed_species_list.txt). Alternatively, you can prepare custom transposon content annotations using the [instructions below](#optional-preparing-custom-transposon-annotations).

### Running FlaHMM

To run FlaHMM, go to the FlaHMM directory and run FlaHMM.py with the following parameters:
- Drosophila species genome assembly to make predictions on (e.g., Dmel.dm6, Dfic.GCF_018152265 or Dper.d101g; please find a full list [here](data/precomputed_species_list.txt)).
- Bin size in kb (2.5, 5 or 10)
- LTR/Gypsy threshold (0.025, 0.05, 0.075, 0.1 or 0.2)

For example:

```
python FlaHMM.py --species Dmel.dm6 --bins 5 --threshold 0.025
```
Results are saved in results folder and can be visualised with the Results_Visualisation.ipynb Jupyter notebook found in the [supplementary repository](https://github.com/Hannon-lab/FlaHMM-supplement/tree/main/08_HTML_visual).

### Output of FlaHMM

Running FlaHMM with the command above would give the following results

```
$ python FlaHMM.py  --species Dmel.dm6 --bins 5 --threshold .025
=========================================================================================
FlaHMM v1.0.0

A tool for flam-like cluster predictions based on LTR/Gypsy transposon content.

If you encounter any issues or have any suggestions, please visit our GitHub issues page:
https://github.com/Hannon-lab/FlaHMM/issues
=========================================================================================
Formatting test data....
100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 29918/29918 [00:02<00:00, 11260.50it/s]
Making flaHMM predictions for Dmel.dm6...

Done. Predictions saved in results folder as results/results_Dmel.dm6_Bin_5k_threshold_0.025.txt.
```

Here, the output file `results/results_Dmel.dm6_Bin_5k_threshold_0.025.txt` is indicated, and has the following structure

```
$ head results/results_Dmel.dm6_Bin_5k_threshold_0.025.txt 
        chr     bin_start       bin_end strand  species_test    region_binary   emission_0.025  pred_0.025      proba_NoCluster_0.025   proba_Cluster_0.025     proba_Centromere_0.025
0       chr2L   0       5000    plus    Dmel.dm6        0       0       0       0.9999719472651957      7.436615274089606e-06   2.0616119612784525e-05
1       chr2L   5000    10000   plus    Dmel.dm6        0       0       0       0.9999847199204476      2.347948548920461e-06   1.2932131042633266e-05
2       chr2L   10000   15000   plus    Dmel.dm6        0       0       0       0.9999911049428629      7.586888007424614e-07   8.13636828857453e-06
3       chr2L   15000   20000   plus    Dmel.dm6        0       0       0       0.9999945944408082      2.6234102246346134e-07  5.143218177404513e-06
4       chr2L   20000   25000   plus    Dmel.dm6        0       0       0       0.9999966175533153      1.0732452310906451e-07  3.2751220558534618e-06
5       chr2L   25000   30000   plus    Dmel.dm6        0       0       0       0.9999978318904452      5.8910512990356805e-08  2.109198985609179e-06
6       chr2L   30000   35000   plus    Dmel.dm6        0       0       0       0.9999985746911687      4.378999174827752e-08   1.3815187994768812e-06
7       chr2L   35000   40000   plus    Dmel.dm6        0       0       0       0.9999990335761717      3.9067539081585073e-08  9.273563859768871e-07
8       chr2L   40000   45000   plus    Dmel.dm6        0       0       0       0.9999993185044835      3.759258416502015e-08   6.439028927488428e-07
```

Specifically, the output data contains the following columns:
1. chr - Chromosome
2. bin_start - Start coordinate
3. bin_end - End coordinate
4. strand - Strand
5. species_test - Name of species (and assembly)
6. region_binary - Known annotation if available (0 = none, 1 = flam(like), 2 = centromere)
7. emission_0.025 - Emission (0 = TEs below threshold, 1 = TEs on one strand, 2 = TEs on both strands)
8. pred_0.025 - Predicted region (0 = none, 1 = flam(like), 2 = centromere)
9. proba_NoCluster_0.025 - Predicted probability for none (0)
10. proba_Cluster_0.025 - Predicted probability for flam(like) (1)
11. proba_Centromere_0.025 - Predicted probability for centromere (2)

### Optional: Preparing custom transposon annotations

To run FlaHMM, you need to calculate transposon content for the genome of interest using RepeatMasker with a suitable transposon library (FASTA). If no transposon library is available for your species of interest, we recommend running EDTA. Alternatively, you could use an existing library (e.g., RepBase) from a closely related species. Detailed instructions are found below:
* [examples/run_EDTA](examples/run_EDTA) - Scripts/instructions to generate EDTA transposon libraries (optional)
* [examples/make_bins](examples/make_bins) - Scripts to calculate LTR/Gypsy content genome-wide per 2.5 kb, 5 kb, or 10 kb bins from RepeatMasker output

If you use your own BED-like input files, put them into the `bins` directory using the following directory structure
```
$ tree bins
bins
├── bins_10k
│   ├── minus
│   │   └── Dmel.dm6.bed
│   ├── plus
│   │   └── Dmel.dm6.bed
└── bins_5k
    ├── minus
    │   └── Dmel.dm6.bed
    └── plus
        └── Dmel.dm6.bed
```

and ensure that the BED-like files have [the following format](examples/make_bins#output-file-format).

### Citation

If you find FlaHMM useful for your research, please cite:

FlaHMM: unistrand *flamenco*-like piRNA cluster prediction in *Drosophila* species using hidden Markov models<br>
Maria-Anna Trapotsi, Jasper van Lopik, Gregory J Hannon, Benjamin Czech Nicholson, Susanne Bornelöv<br>
*NAR Genomics and Bioinformatics* 6:3:lqae119; doi: [https://doi.org/10.1093/nargab/lqae119](https://doi.org/10.1093/nargab/lqae119)
