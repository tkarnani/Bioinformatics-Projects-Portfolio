Tusha Karnani   
May 2025

locAL.py     
Usage: python locAL.py -q [query file] -d [database file] -m [match score] -s [mismatch score] -i [indel score] -a      
Output: for the best local alignment, prints the beginning index of the query sequence, end index of the query sequence, beginning index of the database sequence, end index of the database sequence, score of the best local alignment, and the length of the locally aligned sequence     
Note: the optional -a provides the alignment in a BLAST-like format

locALoutput.txt     
contains the output of running locAL.py on p1-query.txt and p1-db.txt with the -a option

seed_and_extend.py    
Usage: python seed_and_extend.py     
Output: a table containing the time taken to run, number of matches, number of clusters, and number of hits with a score over 50 for alignments with l-mers of lengths 50, 20 and 10

q3a.txt     
contains the output of the seed_and_extend.py function which indexes the database for l-mers of lengths 50, 20 and 10

q3b.txt     
contains the output of the seed_and_extend.py function which identifies l-mer matches between the query and the dataset based on the indexed dictionary created

q3b.txt     
contains the output of the seed_and_extend.py function which identifies clusters of proximal matches

q3d.txt     
contains the output of the seed_and_extend.py function which filters the 4-tuples created based on the score i.e. outputs ones with scores above 50
