import sys

# IMGT CDR ranges
CDR_RANGES = {
    "CDR1": list(range(27, 39)),
    "CDR2": list(range(56, 66)),
    "CDR3": list(range(105, 118)),
}


def read_fasta(fasta_file, chain):
    if chain == "KL":
        chain = "L"

    sequences = {}
    cur_id, cur_seq = None, []

    with open(fasta_file) as f:
        for line in f:
            line = line.strip()

            if line.startswith(">"):
                if cur_id:
                    sequences[cur_id] = "".join(cur_seq)
                cur_id = line[1:]
                cur_seq = []
            else:
                cur_seq.append(line)

        if cur_id:
            sequences[cur_id] = "".join(cur_seq)

    for h, s in sequences.items():
        if chain in h:
            return s

    if len(sequences) == 1:
        return list(sequences.values())[0]

    raise ValueError(f"Chain {chain} not found")


def parse_and_extract(anarci_name, fasta_file, chain):
    seq = read_fasta(fasta_file, chain)

    anarci_file = f'{anarci_name}_{chain}.csv'

    with open(anarci_file) as f:
        header = [x.strip() for x in f.readline().split(",")]

        residues = []

        for line in f:
            cols = [x.strip() for x in line.split(",")]

            if not cols or len(cols) < 10:
                continue

            pos_map = {}

            # build IMGT → AA mapping
            for i in range(13, len(cols)):
                if i >= len(header):
                    break

                try:
                    imgt = int(header[i])
                except:
                    continue

                aa = cols[i]
                if aa != "-":
                    pos_map[imgt] = aa

            # extract each CDR
            for cdr_positions in CDR_RANGES.values():

                # build CDR sequence
                cdr_seq = "".join(pos_map[p] for p in cdr_positions if p in pos_map)

                if not cdr_seq:
                    continue

                start_idx = seq.find(cdr_seq)

                if start_idx == -1:
                    continue

                start = start_idx + 1
                end = start_idx + len(cdr_seq)

                residues.extend(range(start, end + 1))

    return sorted(set(residues))


if __name__ == "__main__":
    anarci_file = sys.argv[1]
    fasta_file = sys.argv[2]

    H = parse_and_extract(anarci_file, fasta_file, "H")
    L = parse_and_extract(anarci_file, fasta_file, "KL")
    #add 1000 to L to avoid overlap with H
    L = [x + 1000 for x in L]

    # merge + sort + unique
    final = sorted(set(H + L))
    print(" ".join(str(x) for x in final))
