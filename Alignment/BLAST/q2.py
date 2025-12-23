import subprocess
import re
import matplotlib.pyplot as plt

numpairs = 50
lenseq = [500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500]
match = 1
mismatch = -30
indel = [0, -20]
lengths = []

print("seq_length, match, mismatch, indel, mean_alignment_length")

for i in indel:
    results = []
    for l in lenseq:
        lengths = []
        for _ in range(numpairs):
            seq1 = subprocess.check_output(['python', 'randomDNA.py', '1', str(l)], text=True).splitlines()[0]
            seq2 = subprocess.check_output(['python', 'randomDNA.py', '1', str(l)], text=True).splitlines()[0]

            with open("query", "w") as f:
                f.write(seq1)
            with open("db", "w") as f:
                f.write(seq2)

            result = subprocess.run(
                ["python", "locAL.py", "-q", "query", "-d", "db",
                 "-m", str(match), "-s", str(mismatch), "-i", str(i)],
                capture_output=True,
                text=True
            )

            numbers = re.findall(r'\d+', result.stdout)
            seqlen = int(numbers[-1])
            lengths.append(seqlen)

        mean_length = sum(lengths) / len(lengths)

        print(l, match, mismatch, i, mean_length)