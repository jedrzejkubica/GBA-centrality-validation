import os
import sys
import logging
import argparse

import multixrank


def main(phenotype: str):
	logger.info("Running MultiXrank")
	multixrank_obj = multixrank.Multixrank(config=f"{phenotype}config.yml", wdir=f"{phenotype}")
	ranking_df = multixrank_obj.random_walk_rank()

	logger.info("Saving scores")
	multixrank_obj.write_ranking(ranking_df, path=f"{phenotype}/output_{phenotype}")

	logger.info("Done!")


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