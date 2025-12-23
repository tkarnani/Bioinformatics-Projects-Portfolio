# import argparse
# import matplotlib.pyplot as plt
# import numpy as np

# def read_fasta(file):
#     with open(file, 'r') as f:
#         lines = f.readlines()
#         return ''.join([line.strip() for line in lines if not line.startswith(">")])

# def get_alignments(query, db, match, mismatch, indel, score_cutoff):
#     m, n = len(query), len(db)
#     prev = [0] * (n + 1)
#     curr = [0] * (n + 1)
#     hits = []

#     for i in range(1, m + 1):
#         for j in range(1, n + 1):
#             diag = prev[j - 1] + (match if query[i - 1] == db[j - 1] else mismatch)
#             up = prev[j] + indel
#             left = curr[j - 1] + indel
#             curr[j] = max(0, diag, up, left)

#             if curr[j] >= score_cutoff:
#                 hits.append((i, j, curr[j]))
#         prev, curr = curr, [0] * (n + 1)

#     return hits

# def traceback(query, db, hit_i, hit_j, match, mismatch, indel):
#     m, n = hit_i, hit_j
#     prev = np.zeros(n + 1, dtype=int)
#     curr = np.zeros(n + 1, dtype=int)
#     prev_arrows = [None] * (n + 1)
#     curr_arrows = [None] * (n + 1)

#     for i in range(1, m + 1):
#         curr[0] = 0
#         curr_arrows[0] = None
#         for j in range(1, n + 1):
#             diag = prev[j-1] + (match if query[i-1] == db[j-1] else mismatch)
#             up = prev[j] + indel
#             left = curr[j-1] + indel
#             best = max(0, diag, up, left)
#             curr[j] = best

#             if best == 0:
#                 curr_arrows[j] = None
#             elif best == diag:
#                 curr_arrows[j] = (i-1, j-1)
#             elif best == up:
#                 curr_arrows[j] = (i-1, j)
#             else:
#                 curr_arrows[j] = (i, j-1)

#         prev, curr = curr, prev
#         prev_arrows, curr_arrows = curr_arrows, prev_arrows

#     aligned_q, aligned_d = [], []
#     i, j = hit_i, hit_j
#     score = prev[j]

#     while i > 0 and j > 0 and prev_arrows[j] is not None:
#         prev_i, prev_j = prev_arrows[j]
#         if prev_i == i-1 and prev_j == j-1:
#             aligned_q.append(query[i-1])
#             aligned_d.append(db[j-1])
#         elif prev_i == i-1:
#             aligned_q.append(query[i-1])
#             aligned_d.append('-')
#         else:
#             aligned_q.append('-')
#             aligned_d.append(db[j-1])
#         i, j = prev_i, prev_j
    
#     return {
#         "db_start": j,
#         "db_end": hit_j-1,
#         "score": score,
#         "aligned_q": ''.join(reversed(aligned_q)),
#         "aligned_d": ''.join(reversed(aligned_d))
#     }

# def remove_overlaps(alignments):
#     alignments = sorted(alignments, key=lambda x: x["score"], reverse=True)
#     final_hits = []
#     for a in alignments:
#         keep = True
#         for b in final_hits:
#             overlap_start = max(a["db_start"], b["db_start"])
#             overlap_end = min(a["db_end"], b["db_end"])
#             db_overlap = (overlap_end - overlap_start + 1) if overlap_start <= overlap_end else 0
#             a_len = a["db_end"] - a["db_start"] + 1
#             if (db_overlap / a_len) >= 0.5:
#                 keep = False
#                 break
#         if keep:
#             final_hits.append(a)
#     print(final_hits)
#     return final_hits

# def plot_threshold_distribution(final_hits):
#     scores = sorted([a['score'] for a in final_hits], reverse=True)
#     thresholds = sorted(set(scores))
#     hits_above = []
#     for t in thresholds:
#         count = 0
#         for s in scores:
#             if s >= t:
#                 count += 1
#         hits_above.append(count)

#     plt.figure(figsize=(10, 6))
#     plt.plot(thresholds, hits_above, marker='o')
#     plt.title("Number of Hits vs. Score Threshold")
#     plt.xlabel("Score Threshold (x)")
#     plt.ylabel("Number of Hits with Score greater than or eaqual to x")
#     plt.grid(True)
#     plt.tight_layout()
#     plt.show()

# def main():
#     parser = argparse.ArgumentParser()
#     parser.add_argument("-q", "--query", required=True)
#     parser.add_argument("-d", "--db", required=True)
#     parser.add_argument("-m", "--match", type=float, default=1)
#     parser.add_argument("-s", "--mismatch", type=float, default=-2)
#     parser.add_argument("-i", "--indel", type=float, default=-2)
#     parser.add_argument("-t", "--threshold", type=int, default=30)
#     parser.add_argument("-p", "--plot", action="store_true")
#     args = parser.parse_args()

#     query = read_fasta(args.query)
#     db = read_fasta(args.db)
#     db = db[0:10001]

#     hits = get_alignments(query, db, args.match, args.mismatch, args.indel, args.threshold)
#     alignments = [traceback(query, db, i, j, args.match, args.mismatch, args.indel)
#                   for (i, j, _) in hits]
#     final_hits = remove_overlaps(alignments)

#     print(f"\nTotal distinct alignments with score ≥ {args.threshold}: {len(final_hits)}\n")
#     for i, a in enumerate(sorted(final_hits, key=lambda x: x['score'], reverse=True)[:15]):
#         print(f"Alignment {i + 1}: Score={a['score']}")
#         print(a['aligned_q'])
#         print(a['aligned_d'])
#         print()

#     if args.plot:
#         plot_threshold_distribution(final_hits)

# if __name__ == "__main__":
#     main()

################################################

import argparse
import matplotlib.pyplot as plt
import numpy as np
import time

def read_fasta(file):
    with open(file, 'r') as f:
        lines = f.readlines()
        return ''.join([line.strip() for line in lines if not line.startswith(">")])

def get_alignments(query, db, match, mismatch, indel, score_cutoff):
    m, n = len(query), len(db)
    score_matrix = np.zeros((m + 1, n + 1), dtype=int)
    traceback_matrix = np.full((m + 1, n + 1), None)

    hits = []

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            diag = score_matrix[i - 1][j - 1] + (match if query[i - 1] == db[j - 1] else mismatch)
            up = score_matrix[i - 1][j] + indel
            left = score_matrix[i][j - 1] + indel
            score = max(0, diag, up, left)

            score_matrix[i][j] = score

            if score == 0:
                traceback_matrix[i][j] = None
            elif score == diag:
                traceback_matrix[i][j] = (i - 1, j - 1)
            elif score == up:
                traceback_matrix[i][j] = (i - 1, j)
            else:
                traceback_matrix[i][j] = (i, j - 1)

            if score >= score_cutoff:
                hits.append((i, j, score))
    print(traceback_matrix)
    return hits, traceback_matrix

def traceback(query, db, i, j, score, traceback_matrix):
    aligned_q, aligned_d = [], []
    end_i, end_j = i, j

    while i > 0 and j > 0 and traceback_matrix[i][j] is not None:
        prev_i, prev_j = traceback_matrix[i][j]
        if prev_i == i - 1 and prev_j == j - 1:
            aligned_q.append(query[i - 1])
            aligned_d.append(db[j - 1])
        elif prev_i == i - 1:
            aligned_q.append(query[i - 1])
            aligned_d.append('-')
        else:
            aligned_q.append('-')
            aligned_d.append(db[j - 1])
        i, j = prev_i, prev_j

    return {
        "db_start": j,
        "db_end": end_j - 1,
        "score": score,
        "aligned_q": ''.join(reversed(aligned_q)),
        "aligned_d": ''.join(reversed(aligned_d))
    }

def remove_overlaps(alignments):
    alignments = sorted(alignments, key=lambda x: x["score"], reverse=True)
    final_hits = []
    for a in alignments:
        keep = True
        for b in final_hits:
            overlap_start = max(a["db_start"], b["db_start"])
            overlap_end = min(a["db_end"], b["db_end"])
            db_overlap = (overlap_end - overlap_start + 1) if overlap_start <= overlap_end else 0
            a_len = a["db_end"] - a["db_start"] + 1
            if (db_overlap / a_len) >= 0.5:
                keep = False
                break
        if keep:
            final_hits.append(a)
    return final_hits

def plot_threshold_distribution(final_hits):
    scores = sorted([a['score'] for a in final_hits], reverse=True)
    thresholds = sorted(set(scores))
    hits_above = []
    for t in thresholds:
        count = sum(s >= t for s in scores)
        hits_above.append(count)

    plt.figure(figsize=(10, 6))
    plt.plot(thresholds, hits_above, marker='o')
    plt.title("Number of Hits vs. Score Threshold")
    plt.xlabel("Score Threshold (x)")
    plt.ylabel("Number of Hits with Score greater than or equal to x")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-q", "--query", required=True)
    parser.add_argument("-d", "--db", required=True)
    parser.add_argument("-m", "--match", type=float, default=1)
    parser.add_argument("-s", "--mismatch", type=float, default=-2)
    parser.add_argument("-i", "--indel", type=float, default=-2)
    parser.add_argument("-t", "--threshold", type=int, default=30)
    parser.add_argument("-p", "--plot", action="store_true")
    args = parser.parse_args()

    start = time.time()

    query = read_fasta(args.query)
    db = read_fasta(args.db)
    db = db[0:10001]

    hits, traceback_matrix = get_alignments(query, db, args.match, args.mismatch, args.indel, args.threshold)
    alignments = [traceback(query, db, i, j, score, traceback_matrix) for (i, j, score) in hits]
    final_hits = remove_overlaps(alignments)

    end = time.time()
    print(end-start)

    print(f"\nTotal distinct alignments with score ≥ {args.threshold}: {len(final_hits)}\n")
    for i, a in enumerate(sorted(final_hits, key=lambda x: x['score'], reverse=True)[:15]):
        print(f"Alignment {i + 1}: Score={a['score']}")
        print(a['aligned_q'])
        print(a['aligned_d'])
        print()

    if args.plot:
        plot_threshold_distribution(final_hits)
    

if __name__ == "__main__":
    main()