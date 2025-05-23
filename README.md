# GBA-centrality-validation

## This repository contains scripts for the validation of GBA centrality

Contents:

- **Part 1: Tissue-enrichment validation**

    We compare the ratio of predicted causal genes enriched in the tissue with the ratio of all genes enriched in the tissue. We compare the two ratios using Fisher's exact test to answer the question: "Are predicted causal genes significantly enriched in the tissue of interest?".

    Part 1.1. GBA centrality
    
    Part 1.2. Random Walk with Restart (RWR)

    Part 1.3. PageRank

    Part 1.4. Comparison of tissue enrichment between methods: GBA centrality vs RWR vs PageRank

- **Part 2: Leave-one-out validation**

    We compare empirical CDFs (cumulative distribution function) for ranks of left-out genes from GBA centrality, RWR and PageRank using the Wilcoxon signed-rank test.

[validation_GBA_centrality.ipynb](validation_GBA_centrality.ipynb)

Perform leave-one-out cross-validation and tissue enrichment validation; compare GBA centrality to Random Walk with Restart and PageRank.

> [!NOTE]
> For tissue enrichment validation download the Expression Atlas from Ensembl reference (v104):
> https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5214/resources/ExperimentDownloadSupplier.RnaSeqBaseline/tpms.tsv
>
> Then add manually a column (after the "ENSG" column) called "tissue_enrichment", which corresonds to the tissue enrichment
> of each gene (i.e., divide tissue expression of each gene by the average expression of that gene in all tissues).

- [validation_parameters.ipynb](validation_parameters.ipynb)

Investigate the choice of the parameters on the methods' performance (GBA centrality, RWR, PageRank)