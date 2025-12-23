Welcome to my bioinformatics projects site. This is where I collect coursework, side projects, and experiments in computational biology as an undergraduate at UC San Diego.​

I’m interested in genomics, data visualization, and using R and Python to explore biological data.

This site is a living portfolio where I document what I’m learning, from basic sequence analysis to more advanced methods, including a Python package to identify enriched motifs.

## Projects

1. ***BLAST-like Alignment Tools*** - Developed a Python-based local sequence alignment pipeline inspired by Smith–Waterman and BLAST. Implemented a configurable local alignment algorithm with BLAST-like output, a random DNA sequence generator, and simulation scripts to study how scoring parameters and sequence length affect expected alignment length between random sequences. Built an iterative BLAST-style search with score thresholds and visualization of hit frequency, and analyzed experimental results across parameter sweeps. This project strengthened my understanding of dynamic programming for sequence alignment, statistical behavior of random alignments, and practical tradeoffs in bioinformatics search algorithms. Find the project [here](https://github.com/tkarnani/Bioinformatics-Projects-Portfolio/tree/main/Alignment).

2. ***Local Alignment and Seed-and-Extend heuristics*** - Implemented a local sequence alignment tool with configurable scoring and BLAST-like output, and extended it with a seed-and-extend heuristic designed for larger, more realistic biological databases. Built an l-mer indexing and matching pipeline to identify seed matches, cluster proximal hits, extend alignments, and filter results by score thresholds, while benchmarking performance across different seed lengths. This project deepened my understanding of local alignment algorithms, BLAST-style heuristics, and the tradeoffs between sensitivity, runtime, and seed size in large-scale sequence search. Find the project [here](https://github.com/tkarnani/Bioinformatics-Projects-Portfolio/tree/main/Alignment/Seed-and-Extend).

3. ***MotifSeeker*** - Developed MotifSeeker, a Python command-line tool for identifying enriched DNA motifs from genomic peak data. The tool takes BED peak files and a reference genome FASTA to detect significantly enriched motifs, emulating core functionality of the HOMER motif analysis pipeline. Implemented statistical enrichment testing with user-defined p-value thresholds, sequence processing, and motif visualization, and packaged the tool for command-line use. This project strengthened my experience with genomic data formats (BED/FASTA), motif discovery workflows, statistical analysis of sequence data, and building reproducible bioinformatics software. Find the repository [here](https://github.com/tkarnani/MotifSeeker/tree/main).


## Contact Me

I’m always happy to connect about bioinformatics, computational biology, research, or collaboration opportunities.

- **Email:** tkarnani@ucsd.edu    
- **GitHub:** https://github.com/tkarnani     
- **LinkedIn:** https://www.linkedin.com/in/tusha-karnani  
