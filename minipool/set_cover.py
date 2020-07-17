from Bio import SeqIO
from collections import defaultdict
from pprint import pprint
from utils import from_json, to_json
from paths import paths
import os
import argparse


def num_triplets(guide):
    count = 0
    for i in range(len(guide) - 3 + 1):
        triplet = guide[i: i + 3]
        if triplet.lower() in [token * 3 for token in ["a", "c", "t", "g", "u"]]:
            count += 1
    return count


def diverse_get_guide(guide_to_seqnames):
    max_covered, _ = max(
        [(len(seqnames), guide)
         for guide, seqnames in guide_to_seqnames.items()],
        key=lambda a: a[0],
    )
    top_guides = [
        guide
        for guide, seqnames in guide_to_seqnames.items()
        if len(seqnames) == max_covered
    ]
    _, top_guide = min(
        [(num_triplets(guide), guide) for guide in top_guides], key=lambda a: a[0]
    )
    return top_guide


def remove_covered_from_guides(guide_to_seqnames, covered):
    new_guides = {}
    for guide in guide_to_seqnames:
        hits = guide_to_seqnames[guide]
        new_hits = hits.difference(covered)
        if len(new_hits) == 0:
            continue
        new_guides[guide] = new_hits
    return new_guides


def set_cover(output_root, seeds):
    print("Reading input file")
    guide_to_seqnames = from_json(
        paths["guide_to_seqnames"].format(output_root))
    guides_to_keep = from_json(paths["guides_to_keep"].format(output_root))

    print("Converting guide2seqnames to use sets")
    new_guide_to_seqnames = {}
    for guide in guides_to_keep:
        new_guide_to_seqnames[guide] = set(guide_to_seqnames[guide])
    guide_to_seqnames = new_guide_to_seqnames

    print("Getting all_sequences")
    all_sequences = set()
    for _, seqnames in guide_to_seqnames.items():
        all_sequences.update(seqnames)

    print("Creating minimal guide set")
    covered = set()
    data = []
    for guide in seeds:
        covered.update(guide_to_seqnames[guide])
        data.append(
            {"guide": guide, "cumulative_sequences_covered": len(covered)})
        print(guide)
        print("Num covered: {} out of {}".format(
            len(covered), len(all_sequences)))
    while len(covered) < len(all_sequences):
        # Make sure that guide_to_seqnames is up to date
        guide_to_seqnames = remove_covered_from_guides(
            guide_to_seqnames, covered)

        # Add one more guide to the minipool
        guide = diverse_get_guide(guide_to_seqnames)
        covered.update(guide_to_seqnames[guide])
        data.append(
            {"guide": guide, "cumulative_sequences_covered": len(covered)})

        print(guide)
        print("Num covered: {} out of {}".format(
            len(covered), len(all_sequences)))

    to_json(paths["final_guides"].format(output_root), data)

    guide_to_seqnames = from_json(
        paths["guide_to_seqnames"].format(output_root))
    guide_to_seqnames = {dat["guide"]                         : guide_to_seqnames[dat["guide"]] for dat in data}
    to_json(paths["final_guide_to_seqnames"].format(
        output_root), guide_to_seqnames)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output_root", default=".")
    parser.add_argument("-s", "--seeds", nargs="*", default=[])
    args = parser.parse_args()

    seeds = args.seeds
    output_root = args.output_root
    output_root = output_root.rstrip("/")
    log_root = "{}/log".format(output_root)

    if not os.path.exists(log_root):
        os.makedirs(log_root)

    set_cover(output_root, seeds)


if __name__ == "__main__":
    main()
