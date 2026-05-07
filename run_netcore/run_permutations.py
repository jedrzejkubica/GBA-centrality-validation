############################################################################################
# Copyright (C) Jędrzej Kubica, Nicolas Thierry-Mieg, 2024-2026
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

import argparse
import logging
import os
import sys

# allow NetCore internal modules to find each other
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "NetCore", "netcore"))
from NetCore.netcore import permutations_test


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run NetCore network permutations.")
    parser.add_argument(
        "--interactome",
        required=True,
        help="path to network file")
    parser.add_argument(
        "--output-path",
        required=True,
        help="directory to save permutation output files")
    parser.add_argument(
        "--net-name",
        default="interactome_human",
        help="name used for the network in output files (optional)")
    parser.add_argument(
        "--num-perm",
        type=int,
        default=1,
        help="number of permutations to run (optional, default: 1)")
    parser.add_argument(
        "--swap-factor",
        type=int,
        default=1,
        help="swap factor used by permutation algorithm (optional, default: 1)")
    parser.add_argument(
        "--num-cores",
        type=int,
        default=8,
        help="number of CPU cores to use (optional, default: 8)")

    args = parser.parse_args()
    return args


def run_permutations(args: argparse.Namespace):
    logger = logging.getLogger(__name__)
    os.makedirs(args.output_path, exist_ok=True)

    logger.info("Running permutations")
    permutations_test.make_network_permutations(
        net_file=args.interactome,
        net_name=args.net_name,
        
        output_path=args.output_path,
        num_perm=args.num_perm,
        swap_factor=args.swap_factor,
        num_cores=args.num_cores)
    
    logger.info("Permutations saved to: %s", args.output_path)
    logger.info("Done!")


def main():
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    args = parse_arguments()
    run_permutations(args)


if __name__ == "__main__":
    main()
