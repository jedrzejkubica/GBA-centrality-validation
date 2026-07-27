# GBA-centrality-validation

This repository contains scripts for the validation of **[GBA centrality](https://github.com/jedrzejkubica/GBA-centrality)** as described in the submitted manuscript. We performed leave-one-out cross-validation (LOO CV) and tissue enrichment validation to compare the performance of GBA centrality with Random Walk with Restart (as implemented in MultiXrank[^1]) and NetCore[^2].


### Step 1. Run LOO CV for GBA centrality

We assume that GBA-centrality is installed as described in [GBA centrality](https://github.com/jedrzejkubica/GBA-centrality) and that input data (interactome, seeds) is prepared in `~/GBA-input/` as described in [GBA centrality Interactome](https://github.com/jedrzejkubica/GBA-centrality/tree/main/Interactome). Scores (`scores_LOO.tsv`) and ranks (`ranks_LOO.tsv`) for left-out genes will be saved in `~/GBA-output/` (`--out`). Logs are written to `log_LOO.txt`.

```
python run_GBA_centrality/leave_one_out.py \
    --network ~/GBA-input/interactome_human.sif \
    --seeds ~/GBA-input/causal_proteins.txt \
    --out ~/GBA-output/ \
    2> ~/GBA-output/log_LOO.txt
```


### Step 2. Run LOO CV for MultiXrank

Create a Python environment and install MultiXrank as described in [https://github.com/anthbapt/multixrank](https://github.com/anthbapt/multixrank).

A default config file is provided at [run_multixrank/default/config.yml](run_multixrank/default/config.yml). All output files will be saved in `~/multixrank-output/`.

```
mkdir ~/multixrank-output/
```

The MultiXrank scoring script takes the interactome SIF file and GBA centrality-derived left-out genes' ranks as input. The interactome SIF file from GBA centrality will be automatically converted to a TSV file (with two columns: node1, node2) as required by MultiXrank. The scoring script produces a file with scores for all genes in the network `multiplex_1.tsv` in `~/multixrank-output/`.

```
python run_multixrank/run_multixrank.py \
    --network ~/GBA-input/interactome_human.sif \
    --GBA_ranks ~/GBA-output/ranks_LOO.tsv \
    --config run_multixrank/default/config.yml \
    --out ~/multixrank-output/
```

Then, the leave-one-out script takes the same input files as the previous step. One subdirectory for each left-out gene is created for intermediate directories and files.

```
python run_multixrank/run_leave_one_out.py \
    --network ~/GBA-input/interactome_human.sif \
    --GBA_ranks ~/GBA-output/ranks_LOO.tsv \
    --config run_multixrank/default/config.yml \
    --out ~/multixrank-output/
```


### Step 3. Run LOO CV for NetCore

Create a Python environment and install NetCore in `run_netcore/` following instructions at [https://github.molgen.mpg.de/barel/NetCore](https://github.molgen.mpg.de/barel/NetCore) (note that NetCore requires Python 3.7).

The interactome TSV file created for MultiXrank will be reused. All output files will be saved in `~/netcore-output/`.

```
mkdir ~/netcore-output/
```

First, NetCore requires running edge permutations on the interactome. This script takes the interactome TSV file and produces a subdirectory `permutations/` in `~/netcore-output/`, required as input in the next step. Logs are written to `log_permut.txt`.

```
python run_netcore/run_permutations.py \
        --interactome ~/multixrank-output/interactome_human.tsv \
        --output-path ~/netcore-output/permutations/ \
        2> ~/netcore-output/log_permut.txt
```

Then, run the NetCore scoring script with an interactome TSV file (`-e`),the seeds file (`-c`), the permutations directory (`-pd`). It writes scores for all genes in the network to `random_walk_weights.txt` in `~/netcore-output/` (`-o`). Both stdout and stderr are captured in `log.txt`.

```
python run_netcore/NetCore/netcore/netcore.py \
        -e ~/multixrank-output/interactome_human.tsv \
        -s ~/GBA-input/causal_proteins.txt \
        -pd ~/netcore-output/permutations/interactome_human_edge_permutations/ \
        -o ~/netcore-output/ \
        1> ~/netcore-output/log.txt \
        2>&1
```

Then, run the leave-one-out bash script, which iterates over left-out genes and re-runs NetCore. The script follows this usage:

```
./run_leave_one_out.sh <interactome> <seeds_dir> <permutation_dir> <out_dir>
```

It takes the following arguments (in this order): interactome TSV file, the temporary files produced by MultiXrank (in `~/multixrank-output/tmp/`) and the permutations subdirectory, the output directory. One subdirectory is created for each left-out gene in `~/netcore-output/`, each containing `random_walk_weights.txt` with scores for all genes in the network, for example:

```
run_netcore/run_leave_one_out.sh \
    ~/multixrank-output/interactome_human.tsv \
    ~/multixrank-output/tmp/ \
    ~/netcore-output/permutations/interactome_human_edge_permutations/ \
    ~/netcore-output/
```


### Part 2. Perform the analyses

This part covers three analyses to compare GBA centrality, MultiXrank and NetCore.

[validation_CDF.py](validation_CDF.py) calculates and plots a cumulative distribution function (CDF) for left-out ranks. CDF curves show the proportion of left-out nodes recovered at or above rank x, for every rank x. The area under the curve (AUC) is calculated and shown in the legend.

```
python validation_CDF.py --help
```

[validation_TE.py](validation_TE.py) compares the ratio of the highest-scoring genes enriched in the tissue of interest with the ratio of all genes. It checks whether the highest-scoring genes are more enriched in the tissue than would be expected by chance.

> [!NOTE]
> For the tissue enrichment validation we downloaded the Expression Atlas from Ensembl reference (v104):
> https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5214/resources/ExperimentDownloadSupplier.RnaSeqBaseline/tpms.tsv
> then we manually added a column (after the "ENSG" column) called "tissue_enrichment", which corresponds to the tissue enrichment of each gene
> (i.e. tissue expression of each gene divided by the average expression of that gene in all tissues).

```
python validation_TE.py --help
```

[validation_ranksVsDeg.py](validation_ranksVsDeg.py) examines the relationship between the degrees of the left-out genes and the differences in their ranks between GBA centrality and MultiXrank, and between GBA centrality and NetCore.

```
python validation_ranksVsDeg.py --help
```


### Dependencies

For validation we used Python 3.9 with the following libraries:
- numpy 1.23
- networkx 3.2
- matplotlib 3.4
- scipy 1.13


### References

[^1]: Baptista, A., Gonzalez, A., & Baudot, A. (2022). Universal multilayer network exploration by random walk with restart. Communications Physics, 5(1), 1–9. https://doi.org/10.1038/s42005-022-00937-9
[^2]: Gal Barel, Ralf Herwig, NetCore: a network propagation approach using node coreness, Nucleic Acids Research, Volume 48, Issue 17, 25 September 2020, Page e98, https://doi.org/10.1093/nar/gkaa639