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
    '''
    # initialize dict to store left-out scores
    scores_left_out = {}

    for gene in ENSG2idx:
        if causal_genes[ENSG2idx[gene]] == 1:
            logger.info("Leaving out %s", gene)
            causal_genes_copy = causal_genes.copy()
            causal_genes_copy[ENSG2idx[gene]] = 0
            scores = GBA_centrality.calculate_scores(interactome, ENSG2idx, causal_genes_copy, alpha, PATH_TO_GBA)
            scores_left_out[gene] = scores[ENSG2idx[gene]]
        else:
            continue
    
    return scores_left_out

def scores_to_TSV(scores, ENSG2gene):
    '''
    arguments:
    - scores: dict with scores for left-out causal genes, key=ENSG, value=score

    prints to STDOUT a TSV with 3 columns: ENSG gene_name score
    '''
    print("ENSG" + "\t" + "GENE" + "\t" + "SCORE")
    for gene in scores:
        score = scores[gene]
        print(gene + "\t" + ENSG2gene[gene] + "\t" + "{:.3g}".format(score))


def main(interactome_file, causal_genes_file, uniprot_file, alpha, weighted, directed, PATH_TO_GBA):

    logger.info("Parsing interactome")
    (interactome, ENSG2idx) = data_parser.parse_interactome(interactome_file, weighted, directed)

    logger.info("Parsing gene-to-ENSG mapping")
    (ENSG2gene, gene2ENSG) = data_parser.parse_uniprot(uniprot_file)

    logger.info("Parsing causal genes")
    causal_genes = data_parser.parse_causal_genes(causal_genes_file, gene2ENSG, ENSG2idx)

    logger.info("Calculating leave-one-out scores")
    scores = leave_one_out(interactome, ENSG2idx, causal_genes, alpha, PATH_TO_GBA)

    logger.info("Printing leave-one-out scores")
    scores_to_TSV(scores, ENSG2gene)

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
        For each causal gene, the method calculates its scores
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

    args = parser.parse_args()

    try:
        main(args.interactome, args.causal, args.uniprot, args.alpha, args.weighted,
             args.directed, PATH_TO_GBA)
    except Exception as e:
        # details on the issue should be in the exception name, print it to stderr and die
        sys.stderr.write("ERROR in " + script_name + " : " + repr(e) + "\n")
        sys.exit(1)
