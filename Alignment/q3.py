import subprocess
import re
import matplotlib.pyplot as plt

numpairs = 50
lenseq = 1000
match = 1
mis_indel = [-2, -1.5, -1.4, -1.3, -1.2, -1.1, -1]
lengths = []

print("seq_length, match, mismatch, indel, mean_alignment_length")

for i in mis_indel:
    lengths = []
    for _ in range(numpairs):
        seq1 = subprocess.check_output(['python', 'randomDNA.py', '1', str(lenseq)], text=True).splitlines()[0]
        seq2 = subprocess.check_output(['python', 'randomDNA.py', '1', str(lenseq)], text=True).splitlines()[0]

        with open("query", "w") as f:
            f.write(seq1)
        with open("db", "w") as f:
            f.write(seq2)
        
        result = subprocess.run(
            ["python", "locAL.py", "-q", "query", "-d", "db",
             "-m", str(match), "-s", str(i), "-i", str(i)],
            capture_output=True,
            text=True
        )
        numbers = re.findall(r'-?\d+', result.stdout)
        seqlen = float(numbers[-1])
        lengths.append(seqlen)

    mean_length = sum(lengths) / len(lengths)
    print(lenseq, match, i, i, mean_length)