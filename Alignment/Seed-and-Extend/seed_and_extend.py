from locAL import read_fasta, align
from collections import defaultdict
from collections import deque
import time

def build_index(sequence, l):
    index = defaultdict(list)
    sequence=sequence[0:10001]
    # print (f"For l = {l}, the indexed substrings are:")
    for i in range(len(sequence) - l + 1):
            kmer = sequence[i:i + l]
            index[kmer].append(i)
    # print(index)
    return index

def find_lmer_matches(query, db_index, l):
    matches = []
    # print (f"For l = {l}, the matches are present at the (query, db) indices:")
    for i in range(len(query) - l + 1):
        kmer = query[i:i + l]
        if kmer in db_index:
            for j in db_index[kmer]:
                matches.append((i, j))
    # print(matches)
    return matches

def cluster_matches_graph(matches, l):

    # print (f"For l = {l}, the 4-tuples for clusters are:")

    n = len(matches)
    visited = [False] * n
    clusters = []

    for index in range(n):
        if visited[index]:
            continue

        cluster = []
        queue = deque([index])
        visited[index] = True

        while queue:
            current = queue.popleft()
            cluster.append(matches[current])
            i1, j1 = matches[current]

            for neighbor in range(n):
                if not visited[neighbor]:
                    i2, j2 = matches[neighbor]
                    if abs(i2 - i1) <= l + 10 and abs(j2 - j1) <= l + 10:
                        visited[neighbor] = True
                        queue.append(neighbor)

        clusters.append(cluster)

    hits = [(min(i for i, _ in cluster), max(i for i, _ in cluster), 
             min(j for _, j in cluster), max(j for _, j in cluster)) 
             for cluster in clusters if cluster]

    # print(hits)
    return hits


def filter_hits_with_local(l, query, db, hits, match=1, mismatch=-2, indel=-2, threshold=50):
    passed = []
    # print (f"For l = {l}, the hits with score above 50 are:")
    for ia, ib, ja, jb in hits:
        q_sub = query[ia:ib+l+1]
        d_sub = db[ja:jb+l+1]
        result = align(q_sub, d_sub, match, mismatch, indel)
        score = result[4]
        if score >= threshold:
            passed.append((ia, ib, ja, jb, score))
            # print(ia, ib, ja, jb)
    return passed

def seed_and_extend(query_path, db_path, l_values, match=1, mismatch=-2, indel=-2, threshold=50):
    query = read_fasta(query_path)
    # print(len(query))
    db = read_fasta(db_path)
    # print(len(db))
    results = []

    for l in l_values:
        start = time.time()
        db_index = build_index(db, l)
        matches = find_lmer_matches(query, db_index, l)
        hits = cluster_matches_graph(matches, l)
        filtered = filter_hits_with_local(l, query, db, hits, match, mismatch, indel, threshold)
        end = time.time()
        results.append({
            'l': l,
            'time': round(end - start, 4),
            'matches': len(matches),
            'hits': len(hits),
            'filtered_hits': len(filtered)
        })
    return results

if __name__ == "__main__":
    query_file = "p4-query.fa"
    query = read_fasta(query_file)
    db_file = "p4-database.fa"
    lmers = [50, 20, 10]
    # lmers = [50, 20, 10, 5]
    db = read_fasta(db_file)

    for l in lmers:
        indices = build_index(db, l)
        find_lmer_matches(query, indices, l)
    results = seed_and_extend(query_file, db_file, lmers)

    print("l-mer\tTime(s)\tMatches\tClusters\tFiltered Hits (Score ≥ 50)")
    for row in results:
        print(f"{row['l']}\t{row['time']}\t{row['matches']}\t{row['hits']}\t\t{row['filtered_hits']}")