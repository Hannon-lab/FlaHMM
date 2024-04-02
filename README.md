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

Alternatively, you can prepare custom transposon annotations using the instructions below.

### Running FlaHMM

To run FlaHMM, go to the FlaHMM directory and run FlaHMM.py with the following parameters:
- Drosophila genome assembly to make predictions on
- Bin size in kb (2.5, 5 or 10)
- LTR/Gypsy threshold (0.025, 0.05, 0.075, 0.1, or 0.2)

For example:

```
python FlaHMM.py --species Dmel.dm6 --bins 5 --threshold 0.025
```
Results are saved in results folder and can be visualised with the Results_Visualisation.ipynb Jupyter notebook found in the [supplementary repository](https://github.com/Hannon-lab/FlaHMM-supplement).

### Optional: Preparing custom transposon annotations

To run the FlaHMM tool you need to prepare the transposon libraries (using EDTA) and calculate genome-wide transposon content. Detailed instructions can be found in:
* examples/run_EDTA - Scripts/instructions to generate EDTA transposon libraries.
* examples/make_bins - Scripts to calculate LTR/Gypsy content genome-wide per 2.5 kb, 5 kb, or 10 kb bins.

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

