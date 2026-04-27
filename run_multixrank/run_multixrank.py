import os
import pathlib
import sys
import logging
import argparse
import subprocess

import multixrank


def parse_interactome(interactome_file, out_dir):
    """
    remove the second column "pp" from interactome SIF and save to TSV

    awk -v OFS='\t' '{print $1, $3}' interactome_file > interactome_human.tsv
    """
    output_file = os.path.join(out_dir, "interactome_human.tsv")
    awk_command = f"awk -v OFS='\\t' '{{print $1, $3}}' {interactome_file} > {output_file}"
    subprocess.run(awk_command, shell=True, check=True)


def parse_seeds(ranks_file, out_dir):
    """
    copy seeds file to out_dir with name seeds.txt,
    reuse GBA LOO output to create seeds file
    """
    output_file = os.path.join(out_dir, "seeds.txt")
    os.makedirs(out_dir, exist_ok=True)

    with open(ranks_file, 'r') as f_in, open(output_file, 'w') as f_out:
        # skip the header
        next(f_in)
        for line in f_in:
            # take the first column (split on whitespace or tabs)
            first_col = line.strip().split()[0]
            f_out.write(first_col + "\n")


def main(network, GBA_ranks, out_dir, config_path):
    logger.info("Parsing interactome and seeds for MultiXrank")

    parse_interactome(network, out_dir)
    parse_seeds(GBA_ranks, out_dir)

    logger.info("Running MultiXrank")
    multixrank_obj = multixrank.Multixrank(config=f"{config_path}", wdir=f"{out_dir}")
    ranking_df = multixrank_obj.random_walk_rank()

    logger.info(f"Saving scores to {out_dir}")
    multixrank_obj.write_ranking(ranking_df, path=out_dir)

    logger.info("Done!")


if __name__ == "__main__":
    (pathToCode, script_name) = os.path.split(os.path.realpath(sys.argv[0]))

    parser = argparse.ArgumentParser(prog=script_name, description="Run MultiXrank")
    parser.add_argument('--network',
                        help=f'''filename (with path) with network in a SIF-like format
                                (3 tab-separated columns: node1 weight/interaction_type node2);''',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--GBA_ranks',
                        help='TSV file (with header) with NODE\tRANK per line',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--config',
                        help='config file',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--out',
                        help='directory where to write the output files',
                        type=pathlib.Path,
                        required=True)

    args = parser.parse_args()

    # create output directory if it doesn't exist
    out_dir = args.out
    os.makedirs(out_dir, exist_ok=True)

    # configure logging, sub-modules will inherit this config
    log_path = os.path.join(out_dir, "log.txt")
    logging.basicConfig(filename=log_path,
                        format='%(asctime)s %(levelname)s %(name)s: %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.DEBUG)
    # set up logger: we want script name rather than 'root'
    logger = logging.getLogger(script_name)
    try:
        main(args.network, args.GBA_ranks, out_dir, config_path=args.config)
    except Exception as e:
        # details on the issue should be in the exception name, print it to stderr and die
        sys.stderr.write("ERROR in " + script_name + " : " + repr(e) + "\n")
        sys.exit(1)
