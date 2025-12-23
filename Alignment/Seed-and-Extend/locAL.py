import argparse

def read_fasta(file):
    with open(file, 'r') as f:
        lines = f.readlines()
        return ''.join([line.strip() for line in lines if not line.startswith(">")])

def align(query, db, match, mismatch, indel):
    m, n = len(query), len(db)
    scoring_matrix = [[0] * (n+1) for _ in range(m+1)]
    arrows = [[None] * (n+1) for _ in range(m+1)]

    max_score = 0
    max_pos = (0, 0)

    for i in range(1, m+1):
        for j in range(1, n+1):
            diag = scoring_matrix[i-1][j-1] + (match if query[i-1] == db[j-1] else mismatch)
            up = scoring_matrix[i-1][j] + indel
            left = scoring_matrix[i][j-1] + indel
            scoring_matrix[i][j] = max(0, diag, up, left)

            if scoring_matrix[i][j] == diag:
                arrows[i][j] = (i-1, j-1)
            elif scoring_matrix[i][j] == up:
                arrows[i][j] = (i-1, j)
            elif scoring_matrix[i][j] == left:
                arrows[i][j] = (i, j-1)

            if scoring_matrix[i][j] > max_score:
                max_score = scoring_matrix[i][j]
                max_pos = (i, j)

    aligned_q, aligned_d = "", ""
    i, j = max_pos
    end_query = i - 1
    end_db = j - 1

    while scoring_matrix[i][j] != 0:
        prev_i, prev_j = arrows[i][j]
        if prev_i == i - 1 and prev_j == j - 1:
            aligned_q = query[i-1] + aligned_q
            aligned_d = db[j-1] + aligned_d
        elif prev_i == i - 1:
            aligned_q = query[i-1] + aligned_q
            aligned_d = '-' + aligned_d
        else:
            aligned_q = '-' + aligned_q
            aligned_d = db[j-1] + aligned_d
        i, j = prev_i, prev_j

    begin_query = i
    begin_db = j
    length = len(aligned_q)

    return begin_query, end_query, begin_db, end_db, max_score, length, aligned_q, aligned_d

def print_alignment(aligned_q, aligned_d, line_width=60):
    match_line = ""

    for q_char, d_char in zip(aligned_q, aligned_d):
        if q_char == d_char:
            match_line += q_char
        else:
            match_line += ' ' 
    for i in range(0, len(aligned_q), line_width):
        q_segment = aligned_q[i:i+line_width]
        m_segment = match_line[i:i+line_width]
        d_segment = aligned_d[i:i+line_width]

        print(f"{q_segment}")
        print(f"{m_segment}")
        print(f"{d_segment}")
        print()

def main():
    parser = argparse.ArgumentParser(description="Local alignment using Smith-Waterman.")
    parser.add_argument("-q", "--query", required=True, help="Query FASTA file")
    parser.add_argument("-d", "--db", required=True, help="Database FASTA file")
    parser.add_argument("-m", "--match", type=float, required=True, help="Match score")
    parser.add_argument("-s", "--mismatch", type=float, required=True, help="Mismatch penalty")
    parser.add_argument("-i", "--indel", type=float, required=True, help="Indel penalty")
    parser.add_argument("-a", "--alignment", action="store_true", help="Output alignment")

    args = parser.parse_args()

    query_seq = read_fasta(args.query)
    db_seq = read_fasta(args.db)

    q_start, q_end, db_start, db_end, score, length, aligned_q, aligned_d = align(query_seq, db_seq, args.match, args.mismatch, args.indel)

    print(f"{q_start}, {q_end}, {db_start}, {db_end}, {score}, {length}")

    if args.alignment:
        print_alignment(aligned_q, aligned_d)

if __name__ == "__main__":
    main()