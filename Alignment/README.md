Tusha Karnani
April 26, 2025

locAL.py
Usage: python locAL.py -q [query file] -d [database file] -m [match score] -s [mismatch score] -i [indel score] -a 
Output: for the best local alignment, prints the beginning index of the query sequence, end index of the query sequence, beginning index of the database sequence, end index of the database sequence, score of the best local alignment, and the length of the locally aligned sequence
Notes: the optional -a provides the alignment in a BLAST-like format

locALoutput.txt
contains the output of running locAL.py on p1-query.txt and p1-db.txt with the -a option

randomDNA.py
Usage: python randomDNA.py [number of sequences] [length of sequences]
Output: randomly generated strings of nucleotides of the given length, constructed with each nucleotide having equal probability of occurring, followed by actual frequencies of all the bases

q2.py
Usage: python q2.py
Output: uses randomDNA.py to generate pairs of sequences, locAL.py to locally align them and provides the length of the random sequences, match score, mismatch score, indel score, and mean length of the locally aligned sequences with those specific parameters

q3.py
Usage: python q3.py
Output: uses randomDNA.py to generate pairs of sequences, locAL.py to locally align them and provides the length of the random sequences, match score, mismatch score, indel score, and mean length of the locally aligned sequences with those specific parameters

iBLAST.py
Usage: python iBLAST.py -q [query file] -d [database file] -m [match score] -s [mismatch score] -i [indel score] -t [score threshold] -p
Notes: the optional -p provides a plot of the number of hits for a score ≥ x as a function for the score threshold x ≤ T

A2 data - Sheet1
contains my data from the experiments I ran for q2 and q3
