import os
import sys
import pathlib
import shutil
import logging
import argparse
import subprocess

import multixrank


def parse_interactome(interactome_file, tmp_dir):
    """
    convert SIF-like interactome to two-column TSV as expected by multixrank
    """
    output_file = os.path.join(tmp_dir, "interactome_human.tsv")
    os.makedirs(tmp_dir, exist_ok=True)
    awk_command = f"awk -v OFS='\\t' '{{print $1, $3}}' {interactome_file} > {output_file}"
    subprocess.run(awk_command, shell=True, check=True)


def parse_seeds(ranks_file, tmp_dir):
    """
    extract seed list from a GBA ranks TSV file into tmp/seeds.txt
    """
    output_file = os.path.join(tmp_dir, "seeds.txt")
    os.makedirs(tmp_dir, exist_ok=True)

    with open(ranks_file, 'r') as f_in, open(output_file, 'w') as f_out:
        header = f_in.readline()
        if header.strip().split()[0].upper() != 'NODE':
            f_in.seek(0)
        for line in f_in:
            if not line.strip():
                continue
            first_col = line.strip().split()[0]
            f_out.write(first_col + "\n")


def save_seeds(seeds, out_file):
    with open(out_file, 'w') as f:
        for gene in seeds:
            f.write(gene + '\n')


def save_config(template_config, out_config, seed_path):
    with open(template_config, 'r') as f_in:
        data = f_in.readlines()

    replaced = False
    for i, line in enumerate(data):
        if line.strip().startswith('seed:'):
            data[i] = f"seed: {seed_path}\n"
            replaced = True
            break

    if not replaced:
        data.insert(0, f"seed: {seed_path}\n")

    with open(out_config, 'w') as f_out:
        f_out.writelines(data)


def main(network, GBA_ranks, out_dir, config_path):
    logger.info("Running leave-one-out for MultiXrank")

    out_dir = os.path.abspath(out_dir)
    tmp_dir = os.path.join(out_dir, 'tmp')
    os.makedirs(tmp_dir, exist_ok=True)

    parse_interactome(network, tmp_dir)
    parse_seeds(GBA_ranks, tmp_dir)

    base_config = os.path.join(out_dir, 'config.yml')
    shutil.copyfile(config_path, base_config)

    seeds_file = os.path.join(tmp_dir, 'seeds.txt')
    with open(seeds_file, 'r') as f:
        causal_genes = [line.strip() for line in f if line.strip()]

    logger.info(f"Found %d causal genes", len(causal_genes))

    ranks_file = os.path.join(out_dir, 'RWR_ranks_LOO.tsv')
    with open(ranks_file, 'w') as f_out:
        f_out.write(f'NODE\tRANK\n')

    for leftOut in causal_genes:
        logger.info(f"Leaving out {leftOut}")

        causal_genes_no_left_out = [g for g in causal_genes if g != leftOut]

        seeds_left_out = os.path.join(out_dir, f"seeds_{leftOut}.txt")
        save_seeds(causal_genes_no_left_out, seeds_left_out)

        config_left_out = os.path.join(out_dir, f"config_{leftOut}.yml")
        save_config(base_config, config_left_out, seeds_left_out)

        multixrank_obj = multixrank.Multixrank(config=config_left_out, wdir=out_dir)
        ranking_df = multixrank_obj.random_walk_rank()

        ranking_df.sort_values(by='score', ascending=False, inplace=True)
        ranking_df.reset_index(drop=True, inplace=True)

        output_path = os.path.join(out_dir, f"output_{leftOut}")
        multixrank_obj.write_ranking(ranking_df, path=output_path)

        matches = ranking_df.index[ranking_df['node'] == leftOut].tolist()
        rank = matches[0] + 1 if matches else 'NA'

        with open(ranks_file, 'a') as f_out:
            f_out.write(f"{leftOut}\t{rank}\n")

    logger.info("Done!")
    logger.info(f"Ranks saved to {ranks_file}")


if __name__ == '__main__':
    (pathToCode, script_name) = os.path.split(os.path.realpath(sys.argv[0]))

    parser = argparse.ArgumentParser(prog=script_name, description='Run leave-one-out for MultiXrank')
    parser.add_argument('--network',
                        help='filename (with path) with network in a SIF-like format',
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

    out_dir = args.out
    os.makedirs(out_dir, exist_ok=True)

    log_path = os.path.join(out_dir, 'log_LOO.txt')
    logging.basicConfig(filename=log_path,
                        format='%(asctime)s %(levelname)s %(name)s: %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.DEBUG)
    logger = logging.getLogger(script_name)

    try:
        main(str(args.network), str(args.GBA_ranks), str(out_dir), str(args.config))
    except Exception as e:
        sys.stderr.write('ERROR in ' + script_name + ' : ' + repr(e) + '\\n')
        sys.exit(1)
