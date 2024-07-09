# FlaHMM 
FlaHMM: unistrand flamenco-like piRNA cluster prediction in Drosophila

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

FlaHMM requires Gypsy transposon annotations across the target genomes. The easiest option will be to download our pre-computed annotations using the following command.

```
cd FlaHMM
wget https://content.cruk.cam.ac.uk/ghlab/Susanne/FlaHMM/bins.tar.gz
tar xf bins.tar.gz
```

This files includes 193 genome assemblies for 119 species. All full list of all pre-processed genome assemblies is available [here](data/precomputed_species_list.txt). Alternatively, you can prepare custom transposon annotations using the [instructions below](#optional-preparing-custom-transposon-annotations).

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

### Optional: Preparing custom transposon annotations

To run the FlaHMM tool you need to calculate transposon content for the genome of interest using RepeatMasker. If no transposon library is present for your species of interest, you may wish to run EDTA to construct transposon libraries; alternatively, you could use RepBase annotations from a closely related species. Detailed instructions are found below:
* examples/run_EDTA - Scripts/instructions to generate EDTA transposon libraries (optional)
* examples/make_bins - Scripts to calculate LTR/Gypsy content genome-wide per 2.5 kb, 5 kb, or 10 kb bins from RepeatMasker output

If you would like to use your own bins, put them into the bins directory and ensure they follow the following structure:
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

### Citation

If you find FlaHMM use, please consider citing it:

FlaHMM: unistrand flamenco-like piRNA cluster prediction in Drosophila species using hidden Markov models

Maria-Anna Trapotsi, Jasper van Lopik, Gregory J Hannon, Benjamin Czech Nicholson, Susanne Bornelöv

*bioRxiv* 2024.05.14.592433; doi: [https://doi.org/10.1101/2024.05.14.592433](https://doi.org/10.1101/2024.05.14.592433)
