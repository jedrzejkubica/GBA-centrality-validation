import os
import sys
import logging
import argparse

import multixrank


def save_seeds(seeds, out_file):
    with open(out_file, "a+") as f:
        for gene in seeds:
            f.write(gene + "\n")


def save_config(phenotype, config, leftOut):
    with open(config, "r") as f_in:
        data = f_in.readlines()

    data[0] = f"seed: seeds_{leftOut}.txt\n"

    out_config_file = f"{phenotype}/config_{leftOut}.yml"

    with open(out_config_file, "w") as f_out:
        f_out.writelines(data)


def main(phenotype: str):
    causal = f"{phenotype}/seeds.txt"
    config = f"{phenotype}/config.yml"
    ranks = f"{phenotype}/RWR_ranks_LOO.tsv"

    logger.info("Running leave-one-out for MultiXrank")

    # read causal genes
    causal_genes = []
    with open(causal, "r") as f:
        data = f.read().splitlines()
        for ENSG in data:
            causal_genes.append(ENSG)
    
    logger.info(f"Found {len(causal_genes)} causal genes for {phenotype}")
	
    for leftOut in causal_genes:
        logger.info(f"Leaving out {leftOut}")

        # for each causal gene, create seeds_causalGene.txt file without the causal gene
        causal_genes_noLeftOut = causal_genes.copy()
        causal_genes_noLeftOut.remove(leftOut)
        save_seeds(causal_genes_noLeftOut, f"{phenotype}/seeds_{leftOut}.txt")

        # for each causal gene, create config_causalGene.yml file with seed: seeds_causalGene.txt
        save_config(config, leftOut)

        # run MultiXrank for each left-out causal gene
        multixrank_obj = multixrank.Multixrank(config=f"{phenotype}/config_{leftOut}.yml",
                                               wdir=f"{phenotype}")
        ranking_df = multixrank_obj.random_walk_rank()

        # sort ranking
        ranking_df.sort_values(by="score", ascending=False, inplace=True)
        ranking_df.reset_index(drop=True, inplace=True)
        multixrank_obj.write_ranking(ranking_df, path=f"{phenotype}/output_{phenotype}_{leftOut}")

        # save rank of left-out
        rank = ranking_df.index[ranking_df["node"] == leftOut].tolist()[0] + 1
    
        with open(ranks, "a+") as f_out:
            f_out.write(f"{leftOut}\t{rank}\n")

    logger.info("Done!")
    logger.info(f"Ranks saved to {ranks}")


if __name__ == "__main__":
	(pathToCode, script_name) = os.path.split(os.path.realpath(sys.argv[0]))
	# configure logging, sub-modules will inherit this config
	logging.basicConfig(format='%(asctime)s %(levelname)s %(name)s: %(message)s',
				datefmt='%Y-%m-%d %H:%M:%S',
				level=logging.DEBUG)
	# set up logger: we want script name rather than 'root'
	logger = logging.getLogger(script_name)

	parser = argparse.ArgumentParser(prog=script_name, description="""Run MultiXrank""")
	parser.add_argument("--pheno", type=str, required=True,
					 help="Phenotype of interest (str)")

	args = parser.parse_args()

	try:
		main(args.pheno)

	except Exception as e:
			# details on the issue should be in the exception name, print it to stderr and die
			sys.stderr.write("ERROR in " + script_name + " : " + repr(e) + "\n")
			sys.exit(1)