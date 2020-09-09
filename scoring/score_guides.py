from Bio import SeqIO
from glob import glob
from Bio.Seq import Seq
from pprint import pprint
import os
import numpy as np
import argparse
import json


def process_tool_output(fname):
    guide_to_scores = {}
    with open(fname, "r") as ropen:
        for idx, line in enumerate(ropen):
            if idx == 0:
                continue
            if len(line.strip()) == 0:
                continue
            split_line = line.strip().split(",")
            score = float(split_line[-2])

            seq = split_line[1]
            seq = str(Seq(seq).reverse_complement())
            for i in range(len(seq) - 22 + 1):
                guide = seq[i : i + 22]
                if guide not in guide_to_scores:
                    guide_to_scores[guide] = []
                guide_to_scores[guide].append(score)

    guide_to_score = {
        guide: np.mean(np.array(scores)) for guide, scores in guide_to_scores.items()
    }

    return guide_to_score


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--tool_output")
    parser.add_argument("-o", "--output_root")

    args = parser.parse_args()
    tool_output = args.tool_output
    output_root = args.output_root

    # Open design tool output files and get scores
    guide_to_score = process_tool_output(tool_output)
    with open("{}/guide_to_score.json".format(output_root), "w") as wopen:
        json.dump(guide_to_score, wopen, indent=4)


if __name__ == "__main__":
    main()
