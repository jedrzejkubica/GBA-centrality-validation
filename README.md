# GBA-centrality-validation

This repository contains scripts for the validation of **[GBA centrality](https://github.com/jedrzejkubica/GBA-centrality)** as described in (x)[^1]. We performed a leave-one-out (LOO) cross-validation and a tissue enrichment validation to compare the performance of GBA centrality, Random Walk with Restart (as implemented in MutliXrank[^2]) and PageRank[^3] (as implemented in NetworkX[^4]). Furthermore, we investigated the choice of the parameters on the methods' performance.


## Part 1. Leave-one-out cross-validation

### Run GBA centrality 

Set up GBA-centrality as described here: https://github.com/jedrzejkubica/GBA-centrality

The script for leave-one-out validation [leave_one_out.py](run_GBA_centrality/leave_one_out.py) calculates scores and ranks for known causal genes as if they were not known to be associated with the disease. Example usage:

```
python leave_one_out.py \
    --interactome ~/GBA-input/interactome_human.sif \
    --causal ~/GBA-input/causal_genes_MMAF.txt \
    --uniprot ~/GBA-input/uniprot_parsed.tsv \
    --alpha 0.5 \
    --scores ~/GBA-output/MMAF/alpha05/scores_LOO.tsv \
    --ranks ~/GBA-output/MMAF/alpha05/ranks_LOO.tsv \
    2>~/GBA-output/MMAF/alpha05/log_LOO.txt
```


### Run MultiXrank

Set up a Python environment for MultiXrank using the following commands (as described here: https://github.com/anthbapt/multixrank):

```
python -m venv --system-site-packages ~/pyenvs/pyEnv_multixrank
source ~/pyenvs/pyEnv_multixrank/bin/activate
pip install --upgrade pip
pip install multixrank
```

Prepare input data as follows:
```
PHENO="MMAF"  # using MMAF phenotype as an example
```

```
cd run_MultiXrank

mkdir $PHENO
cp default/config.yml $PHENO/.  # modify as needed (parameters, paths to interactome and seeds generated in the next steps).
```

Create an interactome TSV file (with 2 columns: ENSG1 ENSG2) using the interactome SIF file from GBA centrality:
```
mkdir -p $PHENO/multiplex/1
cp ~/GBA-input/interactome_human.sif $PHENO/multiplex/1/.
awk -v OFS='\t' '{print $1, $3}' $PHENO/multiplex/1/interactome_human.sif > $PHENO/multiplex/1/interactome_human.tsv
```

Create a seeds file (i.e., disease-associated genes). The most convenient way is to use the output file from leave-one-out of GBA centrality:
```
ln -s ~/GBA-output/$PHENO/alpha05/ranks_LOO.tsv $PHENO/.
awk '{print $1}' $PHENO/ranks_LOO.tsv > $PHENO/seeds.txt
sed -i '1d' $PHENO/seeds.txt
```

The script to run MultiXrank [run_multixrank.py](run_MultiXrank/run_multixrank.py) can be used as follows:
```
python run_multixrank.py MMAF
```

The script for leave-one-out validation [run_leave_one_out.py](run_MultiXrank/run_leave_one_out.py) creates two files for each left-out gene: config_ENSGleftout.yml and seeds_ENSGleftout.txt, then it runs MultiXrank for each left-out gene and saves its rank to RWR_ranks_LOO.tsv. Example usage:
```
python run_leave_one_out.py $PHENO
```


### Part 2. Parse the results

We used the notebook [validation.ipynb](validation.ipynb) to parse the results and:
- compare the ratio of predicted causal genes enriched in the tissue of interest with the ratio of all genes enriched in that tissue

> [!NOTE]
> For the tissue enrichment validation we downloaded the Expression Atlas from Ensembl reference (v104):
> https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5214/resources/ExperimentDownloadSupplier.RnaSeqBaseline/tpms.tsv
> then we added manually a column (after the "ENSG" column) called "tissue_enrichment", which corresonds to the tissue enrichment of each gene
> (i.e., we divide tissue expression of each gene by the average expression of that gene in all tissues).

- compare empirical CDFs (cumulative distribution function) for ranks of left-out genes from GBA centrality, RWR and PageRank

The notebook [validation_parameters.ipynb](validation_parameters.ipynb) can be used to investigate the choice of the parameters on the methods' performance.


### Dependencies

For validation we used Python 3.9 and the following libraries:
- numpy==1.23.5
- networkx==3.2.1
- matplotlib==3.4.3
- scipy==1.9.3


### References

[^1]: TODO
[^2]: Baptista, A., Gonzalez, A., & Baudot, A. (2022). Universal multilayer network exploration by random walk with restart. Communications Physics, 5(1), 1–9. https://doi.org/10.1038/s42005-022-00937-9
[^3]: Page, L., Brin, S., Motwani, R., & Winograd, T. (1999). The PageRank Citation Ranking : Bringing Order to the Web. The Web Conference. https://api.semanticscholar.org/CorpusID:1508503
[^4]: https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.link_analysis.pagerank_alg.pagerank.html
