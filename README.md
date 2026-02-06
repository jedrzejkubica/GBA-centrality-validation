# GBA-centrality-validation

This repository contains scripts for the validation of **[GBA centrality](https://github.com/jedrzejkubica/GBA-centrality)** as described in the submitted manuscript. We performed a leave-one-out (LOO) cross-validation and a tissue enrichment validation to compare the performance of GBA centrality and Random Walk with Restart (as implemented in MutliXrank[^1]).


## Python environment

Install MultiXrank and Python packages for the analyses.

```
python -m venv --system-site-packages ~/pyenvs/pyEnv_GBA-centrality
source ~/pyenvs/pyEnv_GBA-centrality/bin/activate
pip install --upgrade pip
pip install numpy networkx matplotlib scipy
pip install multixrank
```


## Part 1. Leave-one-out cross-validation

We created custom scripts (`run_GBA_centrality/leave_one_out.py` and `run_MultiXrank/run_leave_one_out.py`) for leave-one-out validation.


### Run GBA centrality

We assume that GBA-centrality is installed as described in [GBA centrality](https://github.com/jedrzejkubica/GBA-centrality) and that input data (interactome, uniprot, known causal genes) is prepared as described in [GBA centrality Interactome](https://github.com/jedrzejkubica/GBA-centrality/tree/main/Interactome).

Run [leave_one_out.py](run_GBA_centrality/leave_one_out.py) to calculate scores and ranks for known causal genes as if they were not known to be associated with the disease:

```
python leave_one_out.py \
    --network ~/GBA-input/interactome_human.sif \
    --seeds ~/GBA-input/causal_proteins.txt \
    --scores ~/GBA-output/scores_LOO.tsv \
    --ranks ~/GBA-output/ranks_LOO.tsv \
    2> ~/GBA-output/log_LOO.txt
```


### Run MultiXrank

We will be using MMAF as an example. If needed, adapt the command below for another phenotype:
```
PHENO="MMAF"
```

Prepare input data as follows:

```
cd run_MultiXrank
mkdir $PHENO
cp default/config.yml $PHENO/.
```

Create an interactome TSV file (with 2 columns: protein1 protein2) using the interactome file with which GBA centrality was run:

```
mkdir -p $PHENO/multiplex/1
cp ~/GBA-input/interactome_human.sif $PHENO/multiplex/1/.
awk -v OFS='\t' '{print $1, $3}' $PHENO/multiplex/1/interactome_human.sif > $PHENO/multiplex/1/interactome_human.tsv
```

Create a file with seeds (i.e. known disease causal genes). The most convenient way is to use the output file from leave-one-out of GBA centrality:
```
ln -s ~/GBA-output/ranks_LOO.tsv $PHENO/.
awk '{print $1}' $PHENO/ranks_LOO.tsv > $PHENO/seeds.txt
sed -i '1d' $PHENO/seeds.txt
```

Modify paths to interactome and seeds as needed in `$PHENO/config.yml`.

Run MultiXrank as follows:

```
python run_multixrank.py --pheno MMAF
```

The next script [run_leave_one_out.py](run_MultiXrank/run_leave_one_out.py) creates two files for each left-out gene: config_ENSGleftout.yml and seeds_ENSGleftout.txt, then it runs MultiXrank for each left-out gene and saves its rank to RWR_ranks_LOO.tsv.
```
python run_leave_one_out.py --pheno $PHENO
```


### Part 2. Perform the analyses

Use [validation.ipynb](validation.ipynb) to parse the results and perform the following analyses:
- compare the ratio of predicted causal genes enriched in the tissue of interest with the ratio of all genes enriched in that tissue

> [!NOTE]
> For the tissue enrichment validation we downloaded the Expression Atlas from Ensembl reference (v104):
> https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5214/resources/ExperimentDownloadSupplier.RnaSeqBaseline/tpms.tsv
> then we added manually a column (after the "ENSG" column) called "tissue_enrichment", which corresonds to the tissue enrichment of each gene
> (i.e. we divide tissue expression of each gene by the average expression of that gene in all tissues).

- compare empirical CDFs (cumulative distributions) for ranks of left-out genes from GBA centrality and RWR

WIP: The notebook [validation_parameters.ipynb](validation_parameters.ipynb) can be used to investigate the choice of the parameters on the methods' performance.


### Dependencies

For validation we used Python 3.9 with the following libraries:
- numpy==1.23.5
- networkx==3.2.1
- matplotlib==3.4.3
- scipy==1.9.3


### References

[^1]: Baptista, A., Gonzalez, A., & Baudot, A. (2022). Universal multilayer network exploration by random walk with restart. Communications Physics, 5(1), 1–9. https://doi.org/10.1038/s42005-022-00937-9