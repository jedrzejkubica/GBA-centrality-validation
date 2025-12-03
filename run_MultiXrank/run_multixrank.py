import sys
import multixrank

PHENOTYPE = sys.argv[1]

multixrank_obj = multixrank.Multixrank(config=f"{PHENOTYPE}/config.yml", wdir=f"{PHENOTYPE}")

ranking_df = multixrank_obj.random_walk_rank()

multixrank_obj.write_ranking(ranking_df, path=f"{PHENOTYPE}/output_{PHENOTYPE}")
