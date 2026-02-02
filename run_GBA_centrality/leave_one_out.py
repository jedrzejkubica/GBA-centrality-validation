############################################################################################
# Copyright (C) Jędrzej Kubica, Nicolas Thierry-Mieg, 2024-2025
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
import logging
import pathlib

import argparse

PATH_TO_GBA="/home/kubicaj/Software/GBA-centrality"
sys.path.append(PATH_TO_GBA)
import GBA_centrality
import data_parser


def parse_uniprot(uniprot_file):
    '''
    Parses a tab-seperated Uniprot file produced by GBA_centrality/Interactome/uniprot_parser.py
    with 7 columns (one record per line):
    - Uniprot Primary Accession
    - Taxonomy Identifier
    - ENST (or a comma seperated list of ENSTs)
    - ENSG (or a comma seperated list of ENSGs)
    - Uniprot Secondary Accession (or a comma seperated list of Uniprot Secondary Accessions)
    - GeneID (or a comma seperated list of GeneIDs)
    - Gene Name (or a comma seperated list of Gene Names)

    Returns:
      - ENSG2gene: dict with key=ENSG, value=geneName
      - uniprot2ENSG: dict with key=Primary accession, value=ENSG

    Note: if more than one gene name is associated with a particular ENSG,
          then keeping the first gene name from the list
    '''
    ENSG2gene = {}
    uniprot2ENSG = {}

    try:
        f = open(uniprot_file, 'r')
    except Exception as e:
        logging.error("Opening provided uniprot file %s: %s", uniprot_file, e)
        raise Exception("cannot open provided Uniprot file")

    # skip header
    line = f.readline()
    if not line.startswith("Primary_AC\t"):
        logging.error("uniprot file %s is headerless? expecting headers but got %s",
                      uniprot_file, line)
        raise Exception("Uniprot file problem")

    for line in f:
        split_line = line.rstrip('\r\n').split('\t')

        # if some records are incomplete, die
        if len(split_line) != 7:
            logging.error("uniprot file %s line doesn't have 7 fields: %s",
                          uniprot_file, line)
            raise Exception("Uniprot file problem")

        (AC_primary, TaxID, ENSTs, ENSGs, AC_secondary, GeneIDs, geneNames) = split_line

        # make sure there is at least one ENSG and keep only the first one
        if ENSGs == "":
            continue
        ENSG = ENSGs.split(',')[0]

        # make sure there is at least one gene name and keep only the first one
        if geneNames == "":
            continue
        geneName = geneNames.split(',')[0]

        ENSG2gene[ENSG] = geneName
        uniprot2ENSG[AC_primary] = ENSG

    return(ENSG2gene, uniprot2ENSG)


def leave_one_out(network, node2idx, seeds, seeds_vector, alpha, PATH_TO_GBA, threads):
    '''
    arguments:
    - network: list of "edges", an edge is a tuple (source, dest, weight) where
        source and dest are ints, and weight is a float
    - node2idx: type=dict, key=node, value=unique identifier for the node, these are
        consecutive ints starting at 0
    - seeds: list of floats of length num_nodes, value=1 if node in seeds and 0 otherwise
    - alpha: attenuation coefficient (parameter set by user)
    - PATH_TO_GBA: path to GBA centrality
    - threads: number of threads to use, 0 to use all available cores

    returns:
    - scores_left_out: dict with scores for left-out nodes, key=node, value=score
    - ranks_left_out: dict with ranks for left-out nodes, key=node, value=rank in network
    '''
    # initialize dict to store left-out ranks
    scores_left_out = {}
    ranks_left_out = {}

    for node in seeds:
        logger.info("Leaving out %s", node)
        seeds_vector_copy = seeds_vector.copy()
        seeds_vector_copy[node2idx[node]] = 0
        scores = GBA_centrality.calculate_scores(network, node2idx, seeds_vector_copy, alpha, PATH_TO_GBA, threads)
        # save score
        scores_left_out[node] = scores[node2idx[node]]

        # save rank
        scores_sorted = sorted(scores, reverse=True)
        rank = scores_sorted.index(scores[node2idx[node]]) + 1  # + 1 because ranks start at 1 not 0
        ranks_left_out[node] = rank

    
    return(scores_left_out, ranks_left_out)


def scores_to_TSV(scores, uniprot2ENSG, ENSG2gene, scores_file):
    '''
    arguments:
    - scores: dict with scores for left-out causal genes, key=ENSG, value=score
    - uniprot2ENSG: dict with key=Primary AC, value=ENSG
    - ENSG2gene: dict with key=ENSG, value=gene_name
    - scores_file: path to output file with scores

    saves scores to a TSV file scores_file with 3 columns: ENSG gene_name score
    '''
    with open(scores_file, "w") as f:
        f.write("ENSG" + "\t" + "GENE" + "\t" + "SCORE" + "\n")
        for protein in scores:
            score = scores[protein]
            ENSG = uniprot2ENSG[protein]
            f.write(ENSG + "\t" + ENSG2gene[ENSG] + "\t" + "{:.3g}".format(score) + "\n")


def ranks_to_TSV(ranks, uniprot2ENSG, ENSG2gene, ranks_file):
    '''
    arguments:
    - ranks: dict with ranks for left-out causal genes, key=ENSG, value=rank
    - uniprot2ENSG: dict with key=Primary AC, value=ENSG
    - ENSG2gene: dict with key=ENSG, value=gene_name
    - ranks_file: path to output file with ranks

    saves ranks to a TSV file ranks_file with 3 columns: ENSG gene_name rank
    '''
    with open(ranks_file, "w") as f:
        f.write("ENSG" + "\t" + "GENE" + "\t" + "RANK" + "\n")
        for protein in ranks:
            rank = ranks[protein]
            ENSG = uniprot2ENSG[protein]
            f.write(ENSG + "\t" + ENSG2gene[ENSG] + "\t" + str(rank) + "\n")


def main(network_file, seeds_file, uniprot_file, alpha, weighted, directed,
         scores_file, ranks_file, PATH_TO_GBA, threads):

    logger.info("Parsing interactome")
    (network, node2idx) = data_parser.parse_network(network_file, weighted, directed)

    logger.info("Parsing protein-to-gene mapping")
    (ENSG2gene, uniprot2ENSG) = parse_uniprot(uniprot_file)

    logger.info("Parsing causal genes")
    (seeds, seeds_vector) = data_parser.parse_seeds(seeds_file, node2idx)

    logger.info("Calculating leave-one-out ranks")
    (scores, ranks) = leave_one_out(network, node2idx, seeds, seeds_vector, alpha, PATH_TO_GBA, threads)

    logger.info("Printing leave-one-out scores")
    scores_to_TSV(scores, uniprot2ENSG, ENSG2gene, scores_file)

    logger.info("Printing leave-one-out ranks")
    ranks_to_TSV(ranks, uniprot2ENSG, ENSG2gene, ranks_file)

    logger.info("Done!")


if __name__ == "__main__":
    (pathToCode, script_name) = os.path.split(os.path.realpath(sys.argv[0]))
    # configure logging, sub-modules will inherit this config
    logging.basicConfig(format='%(asctime)s %(levelname)s %(name)s: %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.DEBUG)
    # set up logger: we want script name rather than 'root'
    logger = logging.getLogger(script_name)

    parser = argparse.ArgumentParser(
        prog=script_name,
        description="""
        Leave-one-out validation for GBA centrality.

        For each causal gene, the method calculates its ranks when left-out from the causal gene list.
        """
    )

    parser.add_argument('--network',
                        help='''filename (with path) of network in a SIF-like format
                                (3 tab-separated columns: node1 weight/interaction_type node2), type=str
                                NOTE: second column is either weights (floats in ]0, 1])
                                or a single interaction type (eg "pp"); if weighted, use parameter --weighted''',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--seeds',
                        help='TXT file (without a header) with 1 seed per line',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--uniprot',
                        help='parsed Uniprot file from Interactome/uniprot_parser.py',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--alpha',
                        help='attenuation coefficient (0 < alpha < 1)',
                        default=0.5,
                        type=float)
    parser.add_argument('--weighted',
                        help='use if graph is weighted',
                        action='store_true')  # if present, set the value to True; otherwise False
    parser.add_argument('--directed',
                        help='use if graph is directed',
                        action='store_true')
    parser.add_argument('--scores',
                        help='output file to save scores',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--ranks',
                        help='output file to save ranks',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--threads',
                        help='number of parallel threads to run, default=0 to use all available cores',
                        default=0,
                        type=int)

    args = parser.parse_args()

    try:
        main(args.network, args.seeds, args.uniprot, args.alpha, args.weighted,
             args.directed, args.scores, args.ranks, PATH_TO_GBA, args.threads)
    except Exception as e:
        # details on the issue should be in the exception name, print it to stderr and die
        sys.stderr.write("ERROR in " + script_name + " : " + repr(e) + "\n")
        sys.exit(1)
