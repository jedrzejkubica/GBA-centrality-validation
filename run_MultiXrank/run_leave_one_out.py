import sys
import multixrank

PHENOTYPE = sys.argv[1]

CAUSAL_GENES_FILE = f"{PHENOTYPE}/seeds.txt"
CONFIG_FILE = f"{PHENOTYPE}/config.yml"
RANKS_LEAVE_ONE_OUT_FILE = f"{PHENOTYPE}/RWR_ranks_leave_one_out.tsv"


def save_seeds(gene_list, out_file):
    with open(out_file, "a+") as f:
        for gene in gene_list:
            f.write(gene + "\n")


def save_config(config_file, leftOut):
    with open(config_file, "r") as f_in:
        data = f_in.readlines()

    data[0] = f"seed: seeds_{leftOut}.txt\n"

    out_config_file = f"{PHENOTYPE}/config_{leftOut}.yml"

    with open(out_config_file, "w") as f_out:
        f_out.writelines(data)


# read causal genes
causalGenes = []
with open(CAUSAL_GENES_FILE, "r") as f:
    data = f.read().splitlines()
    for ENSG in data:
        causalGenes.append(ENSG)

print(f"Using {len(causalGenes)} {PHENOTYPE} genes")

for leftOut in causalGenes:
    print(f"Leaving out {leftOut}")

    # for each causal gene, create seeds_causalGene.txt file without the causal gene
    causalGenes_noLeftOut = causalGenes.copy()
    causalGenes_noLeftOut.remove(leftOut)
    save_seeds(causalGenes_noLeftOut, f"{PHENOTYPE}/seeds_{leftOut}.txt")

    # for each causal gene, create config_causalGene.yml file with seed: seeds_causalGene.txt
    save_config(CONFIG_FILE, leftOut)

    # run MultiXrank for each left-out causal gene
    multixrank_obj = multixrank.Multixrank(config=f"{PHENOTYPE}/config_{leftOut}.yml", wdir=f"{PHENOTYPE}")
    ranking_df = multixrank_obj.random_walk_rank()
    # sort ranking
    ranking_df.sort_values(by="score", ascending=False, inplace=True)
    ranking_df.reset_index(drop=True, inplace=True)
    multixrank_obj.write_ranking(ranking_df, path=f"{PHENOTYPE}/output_{PHENOTYPE}_{leftOut}")

    # save rank of left-out gene to ranks_leave_one_out.tsv
    rank = ranking_df.index[ranking_df["node"] == leftOut].tolist()[0] + 1
    with open(RANKS_LEAVE_ONE_OUT_FILE, "a+") as f_scores:
        f_scores.write(f"{leftOut}\t{rank}\n")

print("Done!")
