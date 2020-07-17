from Bio import SeqIO
from Bio.Seq import Seq
import argparse
from collections import defaultdict
import os
from glob import glob
from utils import to_json
from paths import paths


def include_guide(guide):
    """
    Want to filter out some subset of guides? Add your filters to this function.
    """
    if "aaaa" in guide.lower():
        return False
    return True


def prepare_guide_candidates(input_fnames, output_root, reverse_complement):
    guide_to_seqnames_fname = paths["guide_to_seqnames"].format(output_root)
    all_guides_fname = paths["all_guides"].format(output_root)
    guides_to_keep_fname = paths["guides_to_keep"].format(output_root)

    fnames = []
    for input_sequence in input_fnames:
        fnames.extend(glob(input_sequence))

    print("Reading input sequences")
    name_to_seqs = defaultdict(list)
    for fname in fnames:
        for seq_record in SeqIO.parse(fname, "fasta"):
            sequence = str(seq_record.seq).replace("-", "")
            name_to_seqs[seq_record.name].append(sequence)

    print("Extracting guides")
    guide_to_seqnames = defaultdict(set)
    for seqname, sequences in name_to_seqs.items():
        for sequence in sequences:
            for i in range(len(sequence) - 22 + 1):
                guide = sequence[i: i + 22]
                guide = guide.upper()

                if any([char not in ["A", "C", "T", "G", "U"] for char in guide.upper()]):
                    continue

                if reverse_complement:
                    guide = str(Seq(guide).reverse_complement())

                if not include_guide(guide):
                    continue

                guide_to_seqnames[guide].add(seqname)

    to_json(guide_to_seqnames_fname, guide_to_seqnames)

    all_guides = guide_to_seqnames.keys()
    to_json(all_guides_fname, all_guides)
    to_json(guides_to_keep_fname, all_guides)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_sequences",
        nargs="+",
        help="One or more fasta files. Aligned sequences are accepted (gap characters such as '-' are automatically removed). File paths are expanded using glob()",
    )
    parser.add_argument("-r", "--reverse_complement", action="store_true")
    parser.add_argument("-o", "--output_root", default=".")
    args = parser.parse_args()

    input_sequences = args.input_sequences
    reverse_complement = args.reverse_complement
    output_root = args.output_root
    output_root = output_root.rstrip("/")
    log_root = "{}/log".format(output_root)

    if not os.path.exists(log_root):
        os.makedirs(log_root)

    prepare_guide_candidates(input_sequences, output_root, reverse_complement)


if __name__ == "__main__":
    main()
