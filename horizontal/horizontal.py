from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter
from pprint import pprint
from matplotlib import pyplot as plt
import numpy as np
import argparse
import math
from utils import to_json

GUIDE_LENGTH = 22

def include_guide(guide):
    """
    Want to filter out some subset of guides? Add your filters to this function.
    """
    if "AAAA" in guide.upper():
        return False
    return True

def horizontal_pool(input_fname, output_root):
    print("Reading input alignment")
    sequences = []
    for seq_record in SeqIO.parse(input_fname, "fasta"):
        sequences.append(str(seq_record.seq).upper())
    alignment_length = len(sequences[0])

    print("Getting guides")
    guides = []
    for i in range(alignment_length - GUIDE_LENGTH + 1):
        # Get all GUIDE_LENGTH sequences in this window
        all_candidates = [seq[i : i + GUIDE_LENGTH] for seq in sequences]
        # Remove gaps from each sequence
        all_candidates = [
            candidate.replace("-", "").upper() for candidate in all_candidates
        ]
        # Filter out bad guide sequences
        # Filter 1: all characters must be in ACTGU
        # Filter 2: guide must be GUIDE_LENGTH long
        all_candidates = [
            candidate
            for candidate in all_candidates
            if len(candidate) == GUIDE_LENGTH
            if all([char in ["A", "C", "T", "G", "U"] for char in candidate.upper()])
        ]

        if len(all_candidates) == 0:
            continue
        most_common_value, num_occurences = Counter(all_candidates).most_common(1)[0]
        coverage = num_occurences / len(sequences)
        if coverage < 0.5:
            continue

        guide = most_common_value
        if not include_guide(guide):
            continue

        guides.append({"guide": guide, "pos": i, "coverage": coverage})

    to_json("{}/horizontal_pool.json".format(output_root), guides)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_fname")
    parser.add_argument("-o", "--output_root", default=".")

    args = parser.parse_args()
    input_fname = args.input_fname
    output_root = args.output_root
    output_root = output_root.rstrip("/")

    horizontal_pool(input_fname, output_root)


if __name__ == "__main__":
    main()
