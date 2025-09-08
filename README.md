# GBA-centrality-validation

This repository contains scripts for the validation of **[GBA centrality](https://github.com/jedrzejkubica/GBA-centrality)** as described in (x)[^1]. We performed a leave-one-out (LOO) cross-validation and a tissue enrichment validation to compare the performance of GBA centrality, Random Walk with Restart (as implemented in MutliXrank[^2]) and PageRank[^3] (as implemented in NetworkX[^4]) in disease-gene prioritization. Furthermore, we investigated the choice of the parameters on the methods' performance.

### Part 1. Leave-one-out cross-validation

#### 🚀 Run GBA centrality 

Set-up and use GBA-centrality as described in: https://github.com/jedrzejkubica/GBA-centrality

All scripts are in [run_GBA_centrality/](run_GBA_centrality/):
- calculate scores for left-out genes using [leave_one_out_scores.py](run_GBA_centrality/leave_one_out_scores.py)
- calculate ranks for left-out genes using [leave_one_out_ranks.py](run_GBA_centrality/leave_one_out_ranks.py)
- if needed, use [run_all_experiments.sh](run_GBA_centrality/run_all_experiments.sh) to re-run GBA centrality for all six studied phenotypes

example usage:
```
source ~/pyEnv_GBA-centrality/bin/activate
```

```
python leave_one_out_scores.py \
    --interactome ~/GBA-input/interactome_human.sif \
    --causal ~/GBA-input/causal_genes_MMAF.txt \
    --uniprot ~/GBA-input/uniprot_parsed.tsv \
    --alpha 0.5 \
    1>~/GBA-output/scores_LOO.tsv \
    2>~/GBA-output/log_scores_LOO.txt
```

similarily for calculating ranks of left-out genes (instead of scores):
```
python leave_one_out_ranks.py \
    --interactome ~/GBA-input/interactome_human.sif \
    --causal ~/GBA-input/causal_genes_MMAF.txt \
    --uniprot ~/GBA-input/uniprot_parsed.tsv \
    --alpha 0.5 \
    1>~/GBA-output/ranks_LOO.tsv \
    2>~/GBA-output/log_ranks_LOO.txt
```


#### 🚀 Run MultiXrank

#### Set up the environment
```
python -m venv --system-site-packages ~/pyenvs/pyEnv_multixrank
source ~/pyenvs/pyEnv_multixrank/bin/activate
pip install --upgrade pip
pip install multixrank
```

#### How to prepare input data (for MMAF phenotype)
```
mkdir MMAF
cp default/config.yml MAAF/.  # make sure the parameters are OK!
```
We provide an example config.yml file [here](run_MultiXrank/default/config.yml)

Create an interactome TSV file using the interactome SIF file from GBA centrality:
```
mkdir -p MMAF/multiplex/1
cp ~/GBA-input/interactome_human.sif MMAF/multiplex/1/.
awk -v OFS='\t' '{print $1, $3}' interactome_human.sif > interactome_human.tsv  # make sure the interactome TSV file is OK!
```

Create a seeds TXT file with disease-associated genes. The most convenient way is to use the output file from GBA centrality:
```
ln -s ~/GBA-output/ranks_LOO.tsv MMAF/.
awk '{print $1}' MMAF/ranks_LOO.tsv > MMAF/seeds.txt
sed -i '1d' MMAF/seeds.txt
```

Set-up MultiXrank as described in: https://github.com/anthbapt/multixrank
```
# create an output folder
mkdir MMAF
```

All scripts are in [run_MultiXrank/](run_MultiXrank/):
- calculate scores for all genes using [run_multixrank.py](run_MultiXrank/run_multixrank.py)
```
python run_multixrank.py MMAF
```

- calculate ranks for left-out genes using [run_leave_one_out.py](run_MultiXrank/run_leave_one_out.py). It creates config_ENSGleftout.yml and seeds_ENSGleftout.txt for each left-out gene. Then it runs MultiXrank for each left-out gene and saves its rank to RWR_ranks_LOO.tsv
```
python run_leave_one_out.py MMAF
```

### TODO modify run_all_experiments.sh
- if needed, use [run_all_experiments.sh](run_MultiXrank/run_all_experiments.sh) to re-run MultiXrank for all studied phenotypes

### Part 2. Parse the results

The notebook [validation.ipynb](validation.ipynb) was be used to:
- compare the ratio of predicted causal genes enriched in the tissue of interest with the ratio of all genes enriched in that tissue

> [!NOTE]
> For the tissue enrichment validation we downloaded the Expression Atlas from Ensembl reference (v104):
> https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5214/resources/ExperimentDownloadSupplier.RnaSeqBaseline/tpms.tsv
> then we added manually a column (after the "ENSG" column) called "tissue_enrichment", which corresonds to the tissue enrichment of each gene
> (i.e., divide tissue expression of each gene by the average expression of that gene in all tissues).

- compare empirical CDFs (cumulative distribution function) for ranks of left-out genes from GBA centrality, RWR and PageRank

The notebook [validation_parameters.ipynb](validation_parameters.ipynb) can be used to investigate the choice of the parameters on the methods' performance.

### Dependencies

For validation we used Python 3.9 and the following Python libraries:
- numpy==1.23.5
- networkx==3.2.1
- matplotlib==3.4.3
- scipy==1.9.3

### References

[^1]: TODO
[^2]: Baptista, A., Gonzalez, A., & Baudot, A. (2022). Universal multilayer network exploration by random walk with restart. Communications Physics, 5(1), 1–9. https://doi.org/10.1038/s42005-022-00937-9
[^3]: Page, L., Brin, S., Motwani, R., & Winograd, T. (1999). The PageRank Citation Ranking : Bringing Order to the Web. The Web Conference. https://api.semanticscholar.org/CorpusID:1508503
[^4]: https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.link_analysis.pagerank_alg.pagerank.html
