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


def leave_one_out(interactome, ENSG2idx, causal_genes, alpha, PATH_TO_GBA):
    '''
    arguments:
    - interactome: list of "edges", an edge is a tuple (source, dest, weight) where
      source and dest are ints, and weight is a float
    - ENSG2idx: type=dict, key=ENSG, value=unique identifier for the ENSG, these are
      consecutive ints starting at 0
    - causal_genes: list of floats, one per gene, value=1 if gene is causal and 0 otherwise
    - alpha: attenuation coefficient (parameter set by user)

    returns:
    - scores_left_out: dict with scores for left-out causal genes, key=ENSG, value=score
    - ranks_left_out: dict with ranks for left-out causal genes, key=ENSG, value=rank
    '''
    # initialize dict to store left-out ranks
    scores_left_out = {}
    ranks_left_out = {}

    for gene in ENSG2idx:
        if causal_genes[ENSG2idx[gene]] == 1:
            logger.info("Leaving out %s", gene)
            causal_genes_copy = causal_genes.copy()
            causal_genes_copy[ENSG2idx[gene]] = 0
            scores = GBA_centrality.calculate_scores(interactome, ENSG2idx, causal_genes_copy, alpha, PATH_TO_GBA)
            # save score
            scores_left_out[gene] = scores[ENSG2idx[gene]]

            # save rank
            scores_sorted = sorted(scores, reverse=True)
            rank = scores_sorted.index(scores[ENSG2idx[gene]]) + 1  # + 1 because ranks start at 1 not 0
            ranks_left_out[gene] = rank
        else:
            continue
    
    return(scores_left_out, ranks_left_out)


def scores_to_TSV(scores, ENSG2gene, scores_file):
    '''
    arguments:
    - scores: dict with scores for left-out causal genes, key=ENSG, value=score
    - ENSG2gene: dict with key=ENSG, value=geneName
    - scores_file: path to output file with scores

    saves scores to a TSV file scores_file with 3 columns: ENSG gene_name score
    '''
    with open(scores_file, "w") as f:
        f.write("ENSG" + "\t" + "GENE" + "\t" + "SCORE" + "\n")
        for gene in scores:
            score = scores[gene]
            f.write(gene + "\t" + ENSG2gene[gene] + "\t" + "{:.3g}".format(score) + "\n")


def ranks_to_TSV(ranks, ENSG2gene, ranks_file):
    '''
    arguments:
    - ranks: dict with ranks for left-out causal genes, key=ENSG, value=rank
    - ENSG2gene: dict with key=ENSG, value=geneName
    - ranks_file: path to output file with ranks

    saves ranks to a TSV file ranks_file with 3 columns: ENSG gene_name rank
    '''
    with open(ranks_file, "w") as f:
        f.write("ENSG" + "\t" + "GENE" + "\t" + "RANK" + "\n")
        for gene in ranks:
            rank = ranks[gene]
            f.write(gene + "\t" + ENSG2gene[gene] + "\t" + str(rank) + "\n")


def main(interactome_file, causal_genes_file, uniprot_file, alpha, weighted, directed,
         scores_file, ranks_file, PATH_TO_GBA):

    logger.info("Parsing interactome")
    (interactome, ENSG2idx, idx2ENSG) = data_parser.parse_interactome(interactome_file, weighted, directed)

    logger.info("Parsing gene-to-ENSG mapping")
    (ENSG2gene, gene2ENSG, uniprot2ENSG) = data_parser.parse_uniprot(uniprot_file)

    logger.info("Parsing causal genes")
    causal_genes = data_parser.parse_causal_genes(causal_genes_file, gene2ENSG, ENSG2idx)

    logger.info("Calculating leave-one-out ranks")
    (scores, ranks) = leave_one_out(interactome, ENSG2idx, causal_genes, alpha, PATH_TO_GBA)

    logger.info("Printing leave-one-out scores")
    scores_to_TSV(scores, ENSG2gene, scores_file)

    logger.info("Printing leave-one-out ranks")
    ranks_to_TSV(ranks, ENSG2gene, ranks_file)

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
        For each causal gene, the method calculates its ranks
        when left-out from the causal gene list.
        """
    )

    parser.add_argument('--interactome',
                        help='''filename (with path) of interactome in SIF format
                                (3 tab-separated columns: ENSG1 weight/interaction_type ENSG2), type=str
                                NOTE: second column is either weights (floats in [0, 1])
                                or one interaction type (eg "pp"); if weighted, use parameter --weighted''',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--causal',
                        help='TXT file (without a header) with 1 column: gene_name',
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

    args = parser.parse_args()

    try:
        main(args.interactome, args.causal, args.uniprot, args.alpha, args.weighted,
             args.directed, args.scores, args.ranks, PATH_TO_GBA)
    except Exception as e:
        # details on the issue should be in the exception name, print it to stderr and die
        sys.stderr.write("ERROR in " + script_name + " : " + repr(e) + "\n")
        sys.exit(1)
