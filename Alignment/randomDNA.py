import random
import sys

def generate_seq(lenseq):
    return ''.join(random.choices(['A', 'C', 'G', 'T'], k=lenseq))

def main():
    numseq = int(sys.argv[1])
    lenseq = int(sys.argv[2])

    sequences = [generate_seq(lenseq) for _ in range(numseq)]
    for seq in sequences:
        print(seq)

    combined = ''.join(sequences)
    freqs = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    for base in combined:
        freqs[base] += 1
    total = sum(freqs.values())
    for base in 'ACGT':
        print(f'{base}: {freqs[base] / total}')

if __name__ == "__main__":
    main()