from collections import defaultdict

def build_index(sequence, l):
    index = defaultdict(list)
    for i in range(len(sequence) - l + 1):
        kmer = sequence[i:i+l]
        index[kmer].append(i)
    return index

