############################################################################################
# Copyright (C) Jędrzej Kubica, Nicolas Thierry-Mieg, 2026
#
# This file was written by Jędrzej Kubica and Nicolas Thierry-Mieg
# (CNRS, France) Nicolas.Thierry-Mieg@univ-grenoble-alpes.fr
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program.
# If not, see <https://www.gnu.org/licenses/>.
############################################################################################

import os
import sys
import argparse
import pathlib
import logging

import statistics
import numpy
import networkx
import matplotlib.pyplot

sys.path.append("/home/kubicaj/Software/GBA-centrality")
import data_parser

# set up logger, using inherited config, in case we get called as a module
logger = logging.getLogger(__name__)

matplotlib.pyplot.set_loglevel('warning')


def parse_ranks(ranks_file):
    """
    Parses a TSV ranks file (with header) with two columns: NODE, RANK

    Returns:
    - node2rank: dict with key=node, value=rank
    """
    node2rank = {}
    
    with open(ranks_file, 'r') as f:
        header = f.readline()
        if not header.startswith("NODE\t"):
            raise Exception("Ranks file problem, wrong format")

        for line in f:
            split_line = line.rstrip().split("\t")
            if len(split_line) != 2:
                raise Exception(f"Invalid line in ranks file: {line}")

            (node, rank) = split_line
            try:
                rank = int(rank)
            except ValueError:
                raise Exception(f"Invalid rank value in ranks file: {rank} is not an integer")
            node2rank[node] = rank

    return(node2rank)


def generate_random_ranks(n_left_out, network_size):
    """
    Generate ranks with a random classifier

    arguments:
    - n_left_out: number of left-out genes
    - network_size: total number of nodes in the network

    returns:
    - random_ranks: list of average ranks for the left-out genes across n iterations
    """

    # use the integers method of a Generator instance instead of numpy.random.randint
    # https://numpy.org/doc/2.1/reference/random/generated/numpy.random.randint.html
    rng = numpy.random.default_rng()

    # 1D numpy array for left-out gene ranks
    # size=len(left-out genes)
    random_ranks = numpy.zeros(n_left_out, dtype=numpy.uint64)

    n = 10000  # 10,000 iterations for convergence
    for i in range(n):
        # 1D numpy array for left-out gene ranks
        # size=len(left-out genes)
        tmp = rng.integers(low=1, high=network_size, endpoint = True, size=random_ranks.shape, dtype=numpy.uint64)
        tmp = numpy.sort(tmp)
        random_ranks += tmp

    for j in range(len(random_ranks)):
        avg_rank = random_ranks[j] / n
        random_ranks[j] = avg_rank
    
    return(random_ranks)


def ranks_to_curve(ranks, network_size):
    """
    Computes a curve and AUC of rank vs cumulative distributions of left-out nodes

    arguments:
    - ranks: list of ranks for the left-out nodes
    - network_size: total number of nodes in the network

    returns:
    - curve: list of cumulative distributions of left-out nodes
    - AUC: normalized area under the curve
    """
    curve = []
    
    for i in range(1, network_size + 1):
        curve.append(sum(rank <= i for rank in ranks))

    AUC = 0.0
    for i in range(1, len(curve)):
        AUC += (curve[i] + curve[i-1]) / 2

    max_possible = (len(curve) - 1) * len(ranks)  # maximum area: all left-out ranked at the top (rank=1)
    AUC_norm = AUC / max_possible

    return(curve, AUC_norm)


def plot_CDF(GBA_curve, GBA_AUC, random_curve, network_size, out="CDF.png",
             multixrank_curve=None, multixrank_AUC=None, netcore_curve=None, netcore_AUC=None):
    x = range(network_size)
    matplotlib.pyplot.plot(x, GBA_curve, label="GBA centrality (AUC={:.3f})".format(GBA_AUC), color="#D81B60")
    if multixrank_curve:
        matplotlib.pyplot.plot(x, multixrank_curve, label="MultiXrank (AUC={:.3f})".format(multixrank_AUC), color="#FFC107")
    if netcore_curve:
        matplotlib.pyplot.plot(x, netcore_curve, label="NetCore (AUC={:.3f})".format(netcore_AUC), color="#1E88E5")
    matplotlib.pyplot.plot(x, random_curve, label="random classifier", color="#004D40")

    matplotlib.pyplot.xlabel("rank x", fontsize=12)
    matplotlib.pyplot.ylabel("Number of left-out genes where rank <= x", fontsize=12)
    matplotlib.pyplot.xticks(fontsize=11)
    matplotlib.pyplot.yticks(fontsize=11)
    matplotlib.pyplot.legend(loc='lower right', fontsize=12)

    matplotlib.pyplot.savefig(out, dpi=1000)

    logger.info(f"CDF curve saved to {out}")


def main(network_file, GBA_ranks_file, multixrank_ranks_file=None, netcore_ranks_file=None,
         cdf_path=None, weighted=False, directed=False):
    
    logger.info("Parsing network")
    (edge_list, node2idx, idx2node) = data_parser.parse_network(network_file, weighted, directed)

    # construct a networkx Graph from the edge list
    edge_list_no_weights = []
    for edge in edge_list:
        edge_list_no_weights.append((idx2node[edge[0]], idx2node[edge[1]]))
    interactome = networkx.from_edgelist(edge_list_no_weights)

    # calculate the mean and median degree of nodes in the network
    network_degrees = []
    for node in interactome.nodes():
        network_degrees.append(interactome.degree(node))
    logger.info(f"node degree mean: {round(statistics.mean(network_degrees))}, median: {round(statistics.median(network_degrees))}")

    logger.info("Parsing GBA centrality ranks")
    GBA_node2rank = parse_ranks(GBA_ranks_file)
    if multixrank_ranks_file:
        logger.info("Parsing MultiXrank ranks")
        multixrank_node2rank = parse_ranks(multixrank_ranks_file)
        assert len(GBA_node2rank) == len(multixrank_node2rank), "GBA and MultiXrank ranks files have different number of left-out genes"
    if netcore_ranks_file:
        logger.info("Parsing NetCore ranks")
        netcore_node2rank = parse_ranks(netcore_ranks_file)
        assert len(GBA_node2rank) == len(netcore_node2rank), "GBA and NetCore ranks files have different number of left-out genes"
    logger.info(f"Found {len(GBA_node2rank)} left-out genes")

    # calculate the mean and median degree of the left-out genes (seeds)
    seeds_degrees = []
    for protein in GBA_node2rank:
        seeds_degrees.append(interactome.degree(protein))
    logger.info(f"left-out degree mean: {round(statistics.mean(seeds_degrees))}, median: {round(statistics.median(seeds_degrees))}")

    logger.info("Calculating CDF curves and AUCs")
    (GBA_curve, GBA_AUC) = ranks_to_curve(GBA_node2rank.values(), len(node2idx))
    logger.info(f"AUC: {GBA_AUC:.3f} (GBA centrality)")

    optional_kwargs = {}  # dict for plotting
    if multixrank_ranks_file is not None:
        (multixrank_curve, multixrank_AUC) = ranks_to_curve(multixrank_node2rank.values(), len(node2idx))
        logger.info(f"AUC: {multixrank_AUC:.3f} (MultiXrank)")
        optional_kwargs.update({'multixrank_curve': multixrank_curve, 'multixrank_AUC': multixrank_AUC})

    if netcore_ranks_file is not None:
        (netcore_curve, netcore_AUC) = ranks_to_curve(netcore_node2rank.values(), len(node2idx))
        logger.info(f"AUC: {netcore_AUC:.3f} (NetCore)")
        optional_kwargs.update({'netcore_curve': netcore_curve, 'netcore_AUC': netcore_AUC})

    random_ranks = generate_random_ranks(len(GBA_node2rank), len(node2idx))
    (random_curve, random_AUC) = ranks_to_curve(random_ranks, len(node2idx))

    cdf_path.parent.mkdir(parents=True, exist_ok=True)  # Path.parent of a bare filename returns Path("."), and mkdir on "."
    plot_CDF(GBA_curve=GBA_curve,
             GBA_AUC=GBA_AUC,
             random_curve=random_curve,
             network_size=len(interactome.nodes()),
             out=cdf_path,
             **optional_kwargs)


if __name__ == "__main__":
    script_name = os.path.basename(sys.argv[0])
    # configure logging, sub-modules will inherit this config
    logging.basicConfig(format='%(asctime)s %(levelname)s %(name)s: %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.DEBUG)
    # set up logger: we want script name rather than 'root'
    logger = logging.getLogger(script_name)

    parser = argparse.ArgumentParser(
        prog=script_name,
        description="""
        TODO: add description
        """)
    parser.add_argument('--network', required=True)
    parser.add_argument('--GBA_ranks', required=True)
    parser.add_argument('--multixrank_ranks', required=False, default=None)
    parser.add_argument('--netcore_ranks', required=False, default=None)
    parser.add_argument('--cdf', type=pathlib.Path, required=False, default="CDF.png")
    parser.add_argument('--weighted', action='store_true', required=False)
    parser.add_argument('--directed', action='store_true', required=False)

    args = parser.parse_args()

    try:
        main(args.network, args.GBA_ranks,
             multixrank_ranks_file=args.multixrank_ranks,
             netcore_ranks_file=args.netcore_ranks,
             cdf_path=args.cdf,
             weighted=args.weighted,
             directed=args.directed)

    except Exception as e:
        # details on the issue should be in the exception name, print to stderr and die
        sys.stderr.write(f"ERROR in {script_name} : {repr(e)}\n")
        sys.exit(1)