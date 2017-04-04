from collections import defaultdict
import argparse

parser = argparse.ArgumentParser(description='Extract aligned sequences intp FASTA file')
parser.add_argument('filename')
parser.add_argument("--output-file")

def read_raw_data(filename):
    sequence_parts = defaultdict(list)
    first_name = None
    with open(args.filename) as f:
        for line in f:
            line_parts = line.split()
            if len(line_parts) > 1:
                name = line_parts[0]
                if "*" not in name or ":" not in name:
                    # assuming MHC names at least have gene*01:01 format
                    continue
                if not first_name:
                    first_name = name
                sequence_parts[name].extend(line_parts[1:])
    complete_sequences = {
        name: "".join(parts)
        for (name, parts)
        in sequence_parts.items()
    }
    return first_name, complete_sequences

def fill_missing_characters(first_allele_name, sequences):
    first_allele_sequence = sequences[first_allele_name]
    for name, seq in sequences.items():
        if len(seq) > len(first_allele_sequence):
            break
    position_to_default_amino_acid = {
        i: aa for (i, aa) in enumerate(first_allele_sequence)
    }
    results = {}
    for name, sequence in sequences.items():
        sequence = "".join(
            position_to_default_amino_acid[i] if aa in {"*", "-"} else aa
            for i, aa in enumerate(sequence))
        if sequence[-1] == "X":
            sequence = sequence[:-1]
        results[name] = sequence
    return results

if __name__ == "__main__":
    args = parser.parse_args()
    first_allele_name, sequences = read_raw_data(args.filename)
    sequences = fill_missing_characters(first_allele_name, sequences)
    print("Read %d alleles" % len(sequences))
    with open(args.output_file, "w") as f:
        for name, seq in sequences.items():
            f.write(">%s\n%s\n" % (name, seq))
