# GBA-centrality-validation

This repository contains scripts for the validation of **[GBA centrality](https://github.com/jedrzejkubica/GBA-centrality)** as described in (x)[^1]. We performed a leave-one-out cross-validation and a tissue enrichment validation to compare the performance of GBA centrality, Random Walk with Restart (as implemented in MutliXrank[^2]) and PageRank[^3] (as implemented in NetworkX[^4]) in disease-gene prioritization. Furthermore, we investigated the choice of the parameters on the methods' performance.

### Part 1. Leave-one-out cross-validation

#### 🚀 Run GBA centrality 

Set-up and use GBA-centrality as described in: https://github.com/jedrzejkubica/GBA-centrality

All scripts are in [run_GBA_centrality/](run_GBA_centrality/):
- calculate scores for left-out genes using [leave_one_out_scores.py](run_GBA_centrality/leave_one_out_scores.py)
- calculate ranks for left-out genes using [leave_one_out_ranks.py](run_GBA_centrality/leave_one_out_ranks.py)
- if needed, use [run_all_experiments.sh](run_GBA_centrality/run_all_experiments.sh) to re-run GBA centrality for all studied phenotypes

example usage:
```
python leave_one_out_scores.py \
    -i interactome_human.sif \
    --causal_genes_file causal_genes_infertility.tsv \
    --uniprot_file uniprot_parsed.tsv \
    --patho MMAF \
    --alpha 0.5 \
    --d_max 10 \
    1> output/scores_leave_one_out.tsv \
    2> output/log_leave_one_out_scores.txt
```

#### 🚀 Run MultiXrank

#### How to prepare input data

Create an interactome TSV file using the interactome SIF file from GBA centrality by removing the "pp" column in interactome_human.sif.

Create a seeds TXT file with disease-associated genes. The most convenient way is to use the output file from GBA centrality:
```
ln -s run_GBA_centrality/output/MMAF/alpha05_d10/scores_leave_one_out.tsv .
awk '{print $1}' scores_leave_one_out.tsv > seeds.txt
sed -i '1d' seeds.txt  # remove header
```

Set-up MultiXrank as described in: https://github.com/anthbapt/multixrank
```
# create an output folder
mkdir MMAF
```

Create config.yml and modify if necessary. We provide an example file [here](run_MultiXrank/default/config.yml).

All scripts are in [run_MultiXrank/](run_MultiXrank/):
- calculate scores for all genes using [run_multixrank.py](run_MultiXrank/run_multixrank.py)
```
python run_multixrank.py MMAF
```

- calculate ranks for left-out genes using [run_leave_one_out.py](run_MultiXrank/run_leave_one_out.py). It creates config_ENSGleftout.yml and seeds_ENSGleftout.txt for each left-out gene. Then it runs MultiXrank for each left-out gene and saves its rank to RWR_ranks_leave_one_out.tsv
```
python run_leave_one_out.py MMAF
```

- if needed, use [run_all_experiments.sh](run_MultiXrank/run_all_experiments.sh) to re-run MultiXrank for all studied phenotypes

### Part 2. Parse the results

The notebook [validation.ipynb](validation.ipynb) can be used to:
- compare the ratio of predicted causal genes enriched in the tissue of interest with the ratio of all genes enriched in that tissue

> [!NOTE]
> For the tissue enrichment validation we downloaded the Expression Atlas from Ensembl reference (v104):
> https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5214/resources/ExperimentDownloadSupplier.RnaSeqBaseline/tpms.tsv
> then we added manually a column (after the "ENSG" column) called "tissue_enrichment", which corresonds to the tissue enrichment of each gene
> (i.e., divide tissue expression of each gene by the average expression of that gene in all tissues).

- compare empirical CDFs (cumulative distribution function) for ranks of left-out genes from GBA centrality, RWR and PageRank

The notebook [validation_parameters.ipynb](validation_parameters.ipynb) can be used to investigate the choice of the parameters on the methods' performance.

### References

[^1]: https://doi.org/
[^2]: Baptista, A., Gonzalez, A., & Baudot, A. (2022). Universal multilayer network exploration by random walk with restart. Communications Physics, 5(1), 1–9. https://doi.org/10.1038/s42005-022-00937-9
[^3]: Page, L., Brin, S., Motwani, R., & Winograd, T. (1999). The PageRank Citation Ranking : Bringing Order to the Web. The Web Conference. https://api.semanticscholar.org/CorpusID:1508503
[^4]: https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.link_analysis.pagerank_alg.pagerank.html
