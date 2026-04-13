def greedy_clustering(msa, matrix):

    clusters = {}
    for i in range(len(msa)):
        added = False
        for cluster in clusters.values():
            avg_distance = sum(matrix[i][j] for j in cluster) / len(cluster)
            if avg_distance < 0.05:  # Threshold for clustering
                cluster.append(i)
                added = True
                break
        if not added:
            clusters[len(clusters)] = [i]

    return clusters
    