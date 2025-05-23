#!/bin/bash
source ~/pythonvenvs/pyEnv_GBA/bin/activate

ALPHA=0.5
ALPHA_STR=`echo "${ALPHA//.}"`
DMAX=15

# MMAF
echo "Phenotype: MMAF"
python GBA_centrality.py -i Interactome/interactome_human.sif --uniprot_file Interactome/uniprot_parsed.tsv --causal_genes_file ~/workspace/data/input_GBA/causal_genes_infertility.tsv --alpha $ALPHA --d_max $DMAX --patho MMAF 1> output/MMAF/alpha${ALPHA_STR}_d${DMAX}/scores.tsv 2> output/MMAF/alpha${ALPHA_STR}_d${DMAX}/log.txt

python Validation/leave_one_out_scores.py -i Interactome/interactome_human.sif --uniprot_file Interactome/uniprot_parsed.tsv --causal_genes_file ~/workspace/data/input_GBA/causal_genes_infertility.tsv --alpha $ALPHA --d_max $DMAX --patho MMAF 1> output/MMAF/alpha${ALPHA_STR}_d${DMAX}/scores_leave_one_out.tsv 2> output/MMAF/alpha${ALPHA_STR}_d${DMAX}/log_leave_one_out_scores.txt
python Validation/leave_one_out_ranks.py -i Interactome/interactome_human.sif --uniprot_file Interactome/uniprot_parsed.tsv --causal_genes_file ~/workspace/data/input_GBA/causal_genes_infertility.tsv --alpha $ALPHA --d_max $DMAX --patho MMAF 1> output/MMAF/alpha${ALPHA_STR}_d${DMAX}/ranks_leave_one_out.tsv 2> output/MMAF/alpha${ALPHA_STR}_d${DMAX}/log_leave_one_out_ranks.txt

# NOA
echo "Phenotype: NOA"
python GBA_centrality.py -i Interactome/interactome_human.sif --uniprot_file Interactome/uniprot_parsed.tsv --causal_genes_file ~/workspace/data/input_GBA/causal_genes_infertility.tsv --alpha $ALPHA --d_max $DMAX --patho NOA 1> output/NOA/alpha${ALPHA_STR}_d${DMAX}/scores.tsv 2> output/NOA/alpha${ALPHA_STR}_d${DMAX}/log.txt

python Validation/leave_one_out_scores.py -i Interactome/interactome_human.sif --uniprot_file Interactome/uniprot_parsed.tsv --causal_genes_file ~/workspace/data/input_GBA/causal_genes_infertility.tsv --alpha $ALPHA --d_max $DMAX --patho NOA 1> output/NOA/alpha${ALPHA_STR}_d${DMAX}/scores_leave_one_out.tsv 2> output/NOA/alpha${ALPHA_STR}_d${DMAX}/log_leave_one_out_scores.txt
python Validation/leave_one_out_ranks.py -i Interactome/interactome_human.sif --uniprot_file Interactome/uniprot_parsed.tsv --causal_genes_file ~/workspace/data/input_GBA/causal_genes_infertility.tsv --alpha $ALPHA --d_max $DMAX --patho NOA 1> output/NOA/alpha${ALPHA_STR}_d${DMAX}/ranks_leave_one_out.tsv 2> output/NOA/alpha${ALPHA_STR}_d${DMAX}/log_leave_one_out_ranks.txt

# BC
echo "Phenotype: BC"
python GBA_centrality.py -i Interactome/interactome_human.sif --uniprot_file Interactome/uniprot_parsed.tsv --causal_genes_file ~/workspace/data/input_GBA/genes_for_HP_0003002 --alpha $ALPHA --d_max $DMAX --patho BC 1> output/BC/alpha${ALPHA_STR}_d${DMAX}/scores.tsv 2> output/BC/alpha${ALPHA_STR}_d${DMAX}/log.txt

python Validation/leave_one_out_scores.py -i Interactome/interactome_human.sif --uniprot_file Interactome/uniprot_parsed.tsv --causal_genes_file ~/workspace/data/input_GBA/genes_for_HP_0003002 --alpha $ALPHA --d_max $DMAX --patho BC 1> output/BC/alpha${ALPHA_STR}_d${DMAX}/scores_leave_one_out.tsv 2> output/BC/alpha${ALPHA_STR}_d${DMAX}/log_leave_one_out_scores.txt
python Validation/leave_one_out_ranks.py -i Interactome/interactome_human.sif --uniprot_file Interactome/uniprot_parsed.tsv --causal_genes_file ~/workspace/data/input_GBA/genes_for_HP_0003002 --alpha $ALPHA --d_max $DMAX --patho BC 1> output/BC/alpha${ALPHA_STR}_d${DMAX}/ranks_leave_one_out.tsv 2> output/BC/alpha${ALPHA_STR}_d${DMAX}/log_leave_one_out_ranks.txt

# CC
echo "Phenotype: CC"
python GBA_centrality.py -i Interactome/interactome_human.sif --uniprot_file Interactome/uniprot_parsed.tsv --causal_genes_file ~/workspace/data/input_GBA/genes_for_HP_0003003 --alpha $ALPHA --d_max $DMAX --patho CC 1> output/CC/alpha${ALPHA_STR}_d${DMAX}/scores.tsv 2> output/CC/alpha${ALPHA_STR}_d${DMAX}/log.txt

python Validation/leave_one_out_scores.py -i Interactome/interactome_human.sif --uniprot_file Interactome/uniprot_parsed.tsv --causal_genes_file ~/workspace/data/input_GBA/genes_for_HP_0003003 --alpha $ALPHA --d_max $DMAX --patho CC 1> output/CC/alpha${ALPHA_STR}_d${DMAX}/scores_leave_one_out.tsv 2> output/CC/alpha${ALPHA_STR}_d${DMAX}/log_leave_one_out_scores.txt
python Validation/leave_one_out_ranks.py -i Interactome/interactome_human.sif --uniprot_file Interactome/uniprot_parsed.tsv --causal_genes_file ~/workspace/data/input_GBA/genes_for_HP_0003003 --alpha $ALPHA --d_max $DMAX --patho CC 1> output/CC/alpha${ALPHA_STR}_d${DMAX}/ranks_leave_one_out.tsv 2> output/CC/alpha${ALPHA_STR}_d${DMAX}/log_leave_one_out_ranks.txt

# HYPCARD
echo "Phenotype: HYPCARD"
python GBA_centrality.py -i Interactome/interactome_human.sif --uniprot_file Interactome/uniprot_parsed.tsv --causal_genes_file ~/workspace/data/input_GBA/genes_for_HP_0001639 --alpha $ALPHA --d_max $DMAX --patho HYPCARD 1> output/HYPCARD/alpha${ALPHA_STR}_d${DMAX}/scores.tsv 2> output/HYPCARD/alpha${ALPHA_STR}_d${DMAX}/log.txt

python Validation/leave_one_out_scores.py -i Interactome/interactome_human.sif --uniprot_file Interactome/uniprot_parsed.tsv --causal_genes_file ~/workspace/data/input_GBA/genes_for_HP_0001639 --alpha $ALPHA --d_max $DMAX --patho HYPCARD 1> output/HYPCARD/alpha${ALPHA_STR}_d${DMAX}/scores_leave_one_out.tsv 2> output/HYPCARD/alpha${ALPHA_STR}_d${DMAX}/log_leave_one_out_scores.txt
python Validation/leave_one_out_ranks.py -i Interactome/interactome_human.sif --uniprot_file Interactome/uniprot_parsed.tsv --causal_genes_file ~/workspace/data/input_GBA/genes_for_HP_0001639 --alpha $ALPHA --d_max $DMAX --patho HYPCARD 1> output/HYPCARD/alpha${ALPHA_STR}_d${DMAX}/ranks_leave_one_out.tsv 2> output/HYPCARD/alpha${ALPHA_STR}_d${DMAX}/log_leave_one_out_ranks.txt

# DILCARD
echo "Phenotype: DILCARD"
python GBA_centrality.py -i Interactome/interactome_human.sif --uniprot_file Interactome/uniprot_parsed.tsv --causal_genes_file ~/workspace/data/input_GBA/genes_for_HP_0001644 --alpha $ALPHA --d_max $DMAX --patho DILCARD 1> output/DILCARD/alpha${ALPHA_STR}_d${DMAX}/scores.tsv 2> output/DILCARD/alpha${ALPHA_STR}_d${DMAX}/log.txt

python Validation/leave_one_out_scores.py -i Interactome/interactome_human.sif --uniprot_file Interactome/uniprot_parsed.tsv --causal_genes_file ~/workspace/data/input_GBA/genes_for_HP_0001644 --alpha $ALPHA --d_max $DMAX --patho DILCARD 1> output/DILCARD/alpha${ALPHA_STR}_d${DMAX}/scores_leave_one_out.tsv 2> output/DILCARD/alpha${ALPHA_STR}_d${DMAX}/log_leave_one_out_scores.txt
python Validation/leave_one_out_ranks.py -i Interactome/interactome_human.sif --uniprot_file Interactome/uniprot_parsed.tsv --causal_genes_file ~/workspace/data/input_GBA/genes_for_HP_0001644 --alpha $ALPHA --d_max $DMAX --patho DILCARD 1> output/DILCARD/alpha${ALPHA_STR}_d${DMAX}/ranks_leave_one_out.tsv 2> output/DILCARD/alpha${ALPHA_STR}_d${DMAX}/log_leave_one_out_ranks.txt

# OG
echo "Phenotype: OG"
python GBA_centrality.py -i Interactome/interactome_human.sif --uniprot_file Interactome/uniprot_parsed.tsv --causal_genes_file ~/workspace/data/input_GBA/causal_genes_infertility.tsv --alpha $ALPHA --d_max $DMAX --patho OG 1> output/OG/alpha${ALPHA_STR}_d${DMAX}/scores.tsv 2> output/OG/alpha${ALPHA_STR}_d${DMAX}/log.txt

python Validation/leave_one_out_scores.py -i Interactome/interactome_human.sif --uniprot_file Interactome/uniprot_parsed.tsv --causal_genes_file ~/workspace/data/input_GBA/causal_genes_infertility.tsv --alpha $ALPHA --d_max $DMAX --patho OG 1> output/OG/alpha${ALPHA_STR}_d${DMAX}/scores_leave_one_out.tsv 2> output/OG/alpha${ALPHA_STR}_d${DMAX}/log_leave_one_out_scores.txt
python Validation/leave_one_out_ranks.py -i Interactome/interactome_human.sif --uniprot_file Interactome/uniprot_parsed.tsv --causal_genes_file ~/workspace/data/input_GBA/causal_genes_infertility.tsv --alpha $ALPHA --d_max $DMAX --patho OG 1> output/OG/alpha${ALPHA_STR}_d${DMAX}/ranks_leave_one_out.tsv 2> output/OG/alpha${ALPHA_STR}_d${DMAX}/log_leave_one_out_ranks.txt


# PCD
echo "Phenotype: PCD"
python GBA_centrality.py -i Interactome/interactome_human.sif --uniprot_file Interactome/uniprot_parsed.tsv --causal_genes_file ~/workspace/data/input_GBA/causal_genes_infertility.tsv --alpha $ALPHA --d_max $DMAX --patho PCD 1> output/PCD/alpha${ALPHA_STR}_d${DMAX}/scores.tsv 2> output/PCD/alpha${ALPHA_STR}_d${DMAX}/log.txt

python Validation/leave_one_out_scores.py -i Interactome/interactome_human.sif --uniprot_file Interactome/uniprot_parsed.tsv --causal_genes_file ~/workspace/data/input_GBA/causal_genes_infertility.tsv --alpha $ALPHA --d_max $DMAX --patho PCD 1> output/PCD/alpha${ALPHA_STR}_d${DMAX}/scores_leave_one_out.tsv 2> output/PCD/alpha${ALPHA_STR}_d${DMAX}/log_leave_one_out_scores.txt
python Validation/leave_one_out_ranks.py -i Interactome/interactome_human.sif --uniprot_file Interactome/uniprot_parsed.tsv --causal_genes_file ~/workspace/data/input_GBA/causal_genes_infertility.tsv --alpha $ALPHA --d_max $DMAX --patho PCD 1> output/PCD/alpha${ALPHA_STR}_d${DMAX}/ranks_leave_one_out.tsv 2> output/PCD/alpha${ALPHA_STR}_d${DMAX}/log_leave_one_out_ranks.txt

# POF
echo "Phenotype: POF"
python GBA_centrality.py -i Interactome/interactome_human.sif --uniprot_file Interactome/uniprot_parsed.tsv --causal_genes_file ~/workspace/data/input_GBA/causal_genes_infertility.tsv --alpha $ALPHA --d_max $DMAX --patho POF 1> output/POF/alpha${ALPHA_STR}_d${DMAX}/scores.tsv 2> output/POF/alpha${ALPHA_STR}_d${DMAX}/log.txt

python Validation/leave_one_out_scores.py -i Interactome/interactome_human.sif --uniprot_file Interactome/uniprot_parsed.tsv --causal_genes_file ~/workspace/data/input_GBA/causal_genes_infertility.tsv --alpha $ALPHA --d_max $DMAX --patho POF 1> output/POF/alpha${ALPHA_STR}_d${DMAX}/scores_leave_one_out.tsv 2> output/POF/alpha${ALPHA_STR}_d${DMAX}/log_leave_one_out_scores.txt
python Validation/leave_one_out_ranks.py -i Interactome/interactome_human.sif --uniprot_file Interactome/uniprot_parsed.tsv --causal_genes_file ~/workspace/data/input_GBA/causal_genes_infertility.tsv --alpha $ALPHA --d_max $DMAX --patho POF 1> output/POF/alpha${ALPHA_STR}_d${DMAX}/ranks_leave_one_out.tsv 2> output/POF/alpha${ALPHA_STR}_d${DMAX}/log_leave_one_out_ranks.txt

echo "Done!"
